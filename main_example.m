%% ========================================================================
%  Project: dAoA-SLAM (Distributed Angle-of-Arrival SLAM)
%  Script:  Main Simulation Script
%  Description: 
%     This script performs Monte Carlo simulations for the proposed dAoA-SLAM 
%     algorithm using both temporal and spatial correlations.
%
%  Configuration:
%     - Mode: Proposed Method (Spatial + Temporal enabled)
%     - Dataset: scenarioCleanM2_new.mat
% =========================================================================

clear variables; 
close all; 
clc;

% --- Add Required Paths ---
addpath(genpath('.\functions'));

%% ========================================================================
%  1. Global Configuration & Data Loading
% =========================================================================

% Simulation settings
num_runs = 50;              % Number of Monte Carlo runs
save_filename = 'SimulationResults_Proposed.mat';

% Load Scenario Data
fprintf('Loading scenario data...\n');
load('scenarioCleanM2_new.mat', 'trueTrajectory', 'dataVA');

% Trajectory & Time Setup
parameters.maxSteps = 500;
numSteps = parameters.maxSteps;

% --- Pre-processing Trajectory & Visibility ---
% Adjust trajectory coordinate frame
trueTrajectory = (trueTrajectory - [-0.5; 5]) .* 1.8 + [-2; 1.0];
trueTrajectory = trueTrajectory(:, 1:parameters.maxSteps);

% Initialize sensor visibility (assume full visibility initially)
[numSensors, ~] = size(dataVA);
for sensor = 1:numSensors
    dataVA{sensor}.visibility = ones(size(dataVA{sensor}.visibility, 1), length(trueTrajectory));
end

% Manual correction for specific sensor position (based on scenario requirements)
dataVA{2}.positions(1, 5) = 6.7040;

%% ========================================================================
%  2. Parameter Definitions
% =========================================================================

% --- 2.1 Physics & Measurement Parameters ---
de2rad = pi/180; 
parameters.lengthStep = 0.03;       % Step length (m)
parameters.scanTime = 1;            % Scan time (s)
v_max = parameters.lengthStep / parameters.scanTime;

% Driving noise settings
parameters.drivingNoiseVariance = (v_max / 3 / parameters.scanTime)^2;

% Measurement noise settings
AOA_err_deg = 1;                    % AOA Error in degrees
parameters.measurementVariance = (AOA_err_deg * de2rad)^2;
parameters.measurementVarianceLHF = (AOA_err_deg * de2rad)^2;

% --- 2.2 Clutter & Detection Parameters ---
parameters.detectionProbability = 0.95;
parameters.regionOfInterestSize = 15;
parameters.regionOfInterestRangeX = [0, 6];
parameters.regionOfInterestRangeY = [0, 9];

% Clutter density calculations
roi_area = (2 * parameters.regionOfInterestSize)^2;
parameters.meanNumberOfClutter = 0.2;
parameters.clutterIntensity = parameters.meanNumberOfClutter / (2 * pi);

% Birth process parameters
parameters.meanNumberOfBirth = 10^(-4);
parameters.birthIntensity = parameters.meanNumberOfBirth / roi_area;
parameters.meanNumberOfUndetectedAnchors = 5;
parameters.undetectedAnchorsIntensity = parameters.meanNumberOfUndetectedAnchors / roi_area;
parameters.upSamplingFactor = 1;

% --- 2.3 Algorithm Core Parameters ---
parameters.numParticles = 15000;
parameters.known_track = 0;         % 0 = Unknown track, 1 = Known track
parameters.Isplotrt = 0;            % Real-time plotting (0 = Off)
parameters.detectionThreshold = 0.1;
parameters.survivalProbability = 0.999;

% --- 2.4 Anchor & Prior Parameters ---
parameters.priorKnownAnchors{1} = 1;
parameters.priorKnownAnchors{2} = 1;
parameters.priorCovarianceAnchor = 0.1^2 * eye(2);
parameters.anchorRegularNoiseVariance = 1e-4^2;
parameters.UniformRadius_pos = 0.1;
parameters.UniformRadius_vel = 0.1;
parameters.known_agent_inipos = 1;

% --- 2.5 Student-t Distribution Parameters ---
parameters.useTDis = false;
parameters.sigma_a_sq = 0.03;
parameters.beta = 0.2;
parameters.nu = 6;
parameters.unreliabilityThreshold = 1e-4;

% --- 2.6 Method ---
% enabling both Temporal and Spatial dAoA features
method_name = 'Proposed (Full)';

%% ========================================================================
%  3. Monte Carlo Simulation Loop
% =========================================================================

% --- Initialize Result Storage ---
% storage.ospa: [numSensors x numSteps x num_runs]
% storage.rmse: [numSteps x num_runs]
storage.ospa = zeros(numSensors, numSteps, num_runs);
storage.rmse = zeros(numSteps, num_runs);
storage.exec_time = zeros(num_runs, 1);

fprintf('\nStarting Monte Carlo Simulation (%d runs)...\n', num_runs);
fprintf('Method: %s\n', method_name);
fprintf('-----------------------------------------------------------------\n');


for run = 1:num_runs
    fprintf('--> Run %02d/%d | ', run, num_runs);
    
    % 3.1 Data Generation
    rng(run); 
    
    % Initialize priors based on trajectory start
    parameters.priorMean = [trueTrajectory(1:2, 1); 0; 0];
    parameters.trueTraj = trueTrajectory(1:2, :);
    
    % Generate Noisy and Cluttered Measurements
    measurements_clean = generateAngularMeasurements(trueTrajectory, dataVA, parameters);
    clutteredMeasurements = generateClutteredAngMeasurements(measurements_clean, parameters);
    
    % 3.2 Algorithm Execution
    t_start = tic;
    
    [estimatedTrajectory, estimatedAnchors, ~, numEstimatedAnchors, ~] = ...
        BPbasedSLAM_dAoA_v2(dataVA, clutteredMeasurements, parameters, trueTrajectory);
    
    run_time = toc(t_start);
    storage.exec_time(run) = run_time;
    
    % 3.3 Error Calculation (Metric Evaluation)
    positionErrorAgent = zeros(1, numSteps);
    dist_ospa_map = zeros(numSensors, numSteps);
    
    for step = 1:numSteps
        % -- 3.3.1 Agent Localization Error (RMSE) --
        % Only calculate once (using the first sensor loop context for convenience)
        positionErrorAgent(step) = calcDistance_(trueTrajectory(1:2, step), estimatedTrajectory(1:2, step));
        
        % -- 3.3.2 Map Estimation Error (OSPA) --
        for s_idx = 1:numSensors
            trueAnchorPositions = dataVA{s_idx}.positions;
            
            % Extract valid estimated anchor positions based on existence probability
            est_pos = [];
            current_anchors = estimatedAnchors{s_idx, step};
            num_est = numEstimatedAnchors(s_idx, step);
            
            for anchor = 1:num_est
                if current_anchors{anchor}.posteriorExistence >= parameters.detectionThreshold
                    est_pos = [est_pos, current_anchors{anchor}.x]; %#ok<AGROW>
                end
            end
            
            % Calculate OSPA for this sensor at this step
            dist_ospa_map(s_idx, step) = ospa_dist(trueAnchorPositions, est_pos, 10, 1);
        end
    end
    
    % 3.4 Store Results
    storage.ospa(:, :, run) = dist_ospa_map;
    storage.rmse(:, run) = positionErrorAgent';
    
    fprintf('Time: %.2fs | Final RMSE: %.4f m\n', run_time, positionErrorAgent(end));
end

%% ========================================================================
%  4. Statistical Analysis & Output
% =========================================================================

fprintf('\n=================================================================\n');
fprintf('                     SIMULATION RESULTS                          \n');
fprintf('=================================================================\n');

% --- Calculate Statistics ---
% RMSE Analysis
mean_rmse_curve = sqrt(mean(storage.rmse.^2, 2)); % RMSE vs Time
avg_rmse_total  = mean(mean_rmse_curve);          % Global Average
final_rmse      = mean_rmse_curve(end);           % Final Step RMSE

% OSPA Analysis (Average over runs)
mean_ospa_map = mean(storage.ospa, 3);
final_ospa_avg = mean(mean_ospa_map(:, end));     % Average OSPA at final step across sensors

% Timing Analysis
avg_time = mean(storage.exec_time);

% --- Display Summary ---
fprintf('%-20s | %-12s | %-12s | %-10s\n', 'Method', 'Avg RMSE (m)', 'Final OSPA', 'Time/Run(s)');
fprintf('-----------------------------------------------------------------\n');
fprintf('%-20s | %.4f       | %.4f       | %.2f\n', ...
        method_name, avg_rmse_total, final_ospa_avg, avg_time);
fprintf('-----------------------------------------------------------------\n');

% --- Save Results ---
results.parameters = parameters;
results.stats.rmse_curve = mean_rmse_curve;
results.stats.ospa_map = mean_ospa_map;
results.raw_data = storage;

save(save_filename, 'results');
fprintf('Detailed results saved to "%s".\n', save_filename);