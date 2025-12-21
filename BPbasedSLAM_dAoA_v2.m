function [estimatedTrajectory, estimatedAnchors, posteriorParticlesAnchorsstorage, numEstimatedAnchors, posteriorParticlesAgentstorage] = BPbasedSLAM_dAoA_v2(dataVA, clutteredMeasurements, parameters, trueTrajectory)
% BPbasedSLAM_dAoA_v2 Executes SLAM using Belief Propagation (BP) and Particle Filtering.
% Fuses Angle-of-Arrival (AoA) and Differential AoA (dAoA) information.
%
% Inputs:
%   dataVA                - Virtual anchor data structure
%   clutteredMeasurements - Cell array of measurements {step, sensor}
%   parameters            - Algorithm configuration structure
%   trueTrajectory        - Ground truth (for error calc)
%
% Outputs:
%   estimatedTrajectory   - Estimated agent states [4 x numSteps]
%   estimatedAnchors      - Estimated anchor states {sensor, step}{anchor}
%   posteriorParticlesAnchorsstorage - Particle storage for specific steps
%   numEstimatedAnchors   - Count of estimated anchors per step
%   posteriorParticlesAgentstorage   - Agent particle storage (optional)

% Author: YiJin
% ICASSP


    %% 1. Initialization and Allocation
    global_anchor_id_counter = 500; 
    [numSteps, numSensors] = size(clutteredMeasurements);
    numSteps = min(numSteps, parameters.maxSteps);
    numParticles = parameters.numParticles;
    
    % Probability parameters
    detectionProbability = parameters.detectionProbability;       
    survivalProbability = parameters.survivalProbability;         
    birthIntensity = parameters.birthIntensity;                   
    clutterIntensity = parameters.clutterIntensity;               
    undetectedAnchorsIntensity = parameters.undetectedAnchorsIntensity * ones(numSensors, 1);
    unreliabilityThreshold = parameters.unreliabilityThreshold; 
    
    % WLS / Spatial Correction parameters
    parameters.irls_history_size = 100;                 
    parameters.irls_trigger_interval = 5;               
    parameters.irls_min_history_for_correction = 10;  
    parameters.irls_max_solution_jump = 5.0;            % 
    parameters.irls_min_direction_consistency_ratio = 0.7; %
    
    parameters.enable_keyframe_selection = true;
    parameters.keyframe_min_distance = 1.0;     % m
    parameters.keyframe_min_angle = 5 * (pi/180); % rad

    % Memory Allocation
    estimatedTrajectory = zeros(4, numSteps);
    numEstimatedAnchors = zeros(2, numSteps);
    execTimePerStep = zeros(numSteps, 1);
    
    storing_idx = 50:50:numSteps; % Steps to save full particle data
    posteriorParticlesAnchorsstorage = cell(1, length(storing_idx));
    posteriorParticlesAgentstorage = cell(1, length(storing_idx));

    %% 2. State Initialization
    % --- Initialize Agent Particles ---
    if(parameters.known_track)
        posteriorParticlesAgent = repmat([trueTrajectory(:,1); 0; 0], 1, numParticles);
    else 
        posteriorParticlesAgent = zeros(4, numParticles);
        % Pos: Uniform circle; Vel: Uniform distribution
        posteriorParticlesAgent(1:2, :) = drawSamplesUniformlyCirc(parameters.priorMean(1:2), parameters.UniformRadius_pos, numParticles);
        posteriorParticlesAgent(3:4, :) = repmat(parameters.priorMean(3:4), 1, numParticles) + ...
            2 * parameters.UniformRadius_vel * rand(2, numParticles) - parameters.UniformRadius_vel;
    end
    
    estimatedTrajectory(:, 1) = mean(posteriorParticlesAgent, 2);

    % --- Initialize Anchors (Prior Knowledge) ---
    if parameters.priorKnownAnchors{1}
        [estimatedAnchors, posteriorParticlesAnchors] = initAnchors(parameters, dataVA, numSteps, numSensors);
        
        % Initialize history buffers for known anchors
        for sensor = 1:numSensors
            numKnownAnchors = size(estimatedAnchors{sensor, 1}, 2);
            for anchor = 1:numKnownAnchors
                posteriorParticlesAnchors{sensor}{anchor}.history = cell(1, parameters.irls_history_size);
                posteriorParticlesAnchors{sensor}{anchor}.history_idx = 0;
                posteriorParticlesAnchors{sensor}{anchor}.correction_count = 0;
            end
        end
    end
    
    for sensor = 1:numSensors
        numEstimatedAnchors(sensor, 1) = size(estimatedAnchors{sensor, 1}, 2);
    end

    %% 3. Main Time Loop
    for step = 2:numSteps 
        tic; 
        
        % --- 3.1 Agent Prediction ---
        if(parameters.known_track)
            predictedParticlesAgent = repmat([trueTrajectory(:,step); 0; 0], 1, numParticles);
        else
            predictedParticlesAgent = performPrediction(posteriorParticlesAgent, parameters); 
        end
        
        weightsSensors = nan(numParticles, numSensors); 
        
        % --- 3.2 Sensor Loop ---
        for sensor = 1:numSensors
            estimatedAnchors{sensor, step} = estimatedAnchors{sensor, step-1}; 
            measurements = clutteredMeasurements{step, sensor};
            numMeasurements = size(measurements, 2);
            
            % Update intensity of undetected anchors
            undetectedAnchorsIntensity(sensor) = undetectedAnchorsIntensity(sensor) * survivalProbability + birthIntensity;
            
            % --- 3.2.1 Anchor Prediction ---
            [predictedParticlesAnchors, weightsAnchor] = predictAnchors(posteriorParticlesAnchors{sensor}, parameters);
            
            % Generate new anchor hypotheses from measurements
            [newParticlesAnchors, newInputBP, parameters] = generateNewAnchorsAngular(...
                measurements, undetectedAnchorsIntensity(sensor), predictedParticlesAgent, parameters);
            
            % --- 3.2.2 Data Association & Measurement Prediction ---
            [predictedMeasurements, predictedUncertainties, predictedAngle] = predictMeasurementsAngular_t(...
                predictedParticlesAgent, predictedParticlesAnchors, weightsAnchor);
            
            % Calculate BP messages / Association probabilities
            [associationProbabilities, ~, messagelhfRatios, messagesNew] = calculateAssociationProbabilitiesGAAngular(...
                measurements, predictedMeasurements, predictedUncertainties, weightsAnchor, newInputBP, parameters);
            
            agent_pose_cov_current = cov(predictedParticlesAgent(1:2, :)');
            
            % --- 3.2.3 Legacy Anchor Update ---
            numAnchors = size(predictedParticlesAnchors, 3);
            weights = zeros(numParticles, numAnchors);
            
            for anchor = 1:numAnchors
                % Initialize weights with miss-detection probability
                weights(:, anchor) = repmat((1 - detectionProbability), numParticles, 1);
                
                % Update with measurements
                for m = 1:numMeasurements
                    measVar = measurements(2, m);
                    factor = 1/sqrt(2*pi*measVar) * detectionProbability / clutterIntensity;
                    lik = exp(-1/(2*measVar) * (measurements(1, m) - predictedAngle(:, anchor)).^2);
                    weights(:, anchor) = weights(:, anchor) + factor * messagelhfRatios(m, anchor) * lik;
                end
                
                % --- 3.2.4 dAoA Kalman Filter Update (Anchor State) ---
                % Track angular velocity of the anchor relative to the agent
                
                % KF Prediction
                kf_state_prev = posteriorParticlesAnchors{sensor}{anchor}.kf_state;
                kf_cov_prev   = posteriorParticlesAnchors{sensor}{anchor}.kf_covariance;
                
                dt = parameters.scanTime;
                F_kf = [1 dt; 0 1];
                Q_kf = [dt^4/4, dt^3/2; dt^3/2, dt^2] * parameters.sigma_a_sq;
                kf_state_pred = F_kf * kf_state_prev;
                kf_cov_pred   = F_kf * kf_cov_prev * F_kf' + Q_kf;
                
                % KF Update Decision (Based on Association Probability)
                [max_assoc_prob, m_idx] = max(associationProbabilities(2:end, anchor));
                is_associated = ~isempty(max_assoc_prob) && ...
                                (max_assoc_prob > associationProbabilities(1, anchor)) && ...
                                (max_assoc_prob > 0.1);
                
                if is_associated 
                    H_kf = [1 0];
                    R_kf = parameters.measurementVariance;
                    y_res = measurements(1, m_idx) - H_kf * kf_state_pred;
                    S_kf = H_kf * kf_cov_pred * H_kf' + R_kf;
                    K_kf = kf_cov_pred * H_kf' / S_kf;
                    
                    kf_state_updated = kf_state_pred + K_kf * y_res;
                    kf_cov_updated = (eye(2) - K_kf * H_kf) * kf_cov_pred;
                    
                    % Record history for WLS
                    if isfield(posteriorParticlesAnchors{sensor}{anchor}, 'history')
                        h_idx = mod(posteriorParticlesAnchors{sensor}{anchor}.history_idx, parameters.irls_history_size) + 1;
                        obs_data.aoa = measurements(1, m_idx);
                        obs_data.aoa_var = measurements(2, m_idx);
                        obs_data.agent_pose = mean(predictedParticlesAgent, 2);
                        obs_data.agent_cov = agent_pose_cov_current;
                        
                        posteriorParticlesAnchors{sensor}{anchor}.history{h_idx} = obs_data;
                        posteriorParticlesAnchors{sensor}{anchor}.history_idx = h_idx;
                    end
                else
                    kf_state_updated = kf_state_pred;
                    kf_cov_updated = kf_cov_pred;
                end
                
                posteriorParticlesAnchors{sensor}{anchor}.kf_state = kf_state_updated;
                posteriorParticlesAnchors{sensor}{anchor}.kf_covariance = kf_cov_updated;

                % --- 3.2.5 Weight Normalization & Resampling ---
                predExistence = sum(weightsAnchor(:, anchor)); 
                aliveUpdate = sum(predExistence * 1/numParticles * weights(:, anchor));
                deadUpdate = 1 - predExistence;
                posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence = aliveUpdate / (aliveUpdate + deadUpdate);
                
                % Resample anchor particles
                idx_res = resampleSystematic(weights(:, anchor) / sum(weights(:, anchor)), numParticles);
                posteriorParticlesAnchors{sensor}{anchor}.x = predictedParticlesAnchors(:, idx_res(1:numParticles), anchor);
                posteriorParticlesAnchors{sensor}{anchor}.w = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence / numParticles * ones(numParticles, 1);
                
                % --- 3.2.6 Spatial dAoA: WLS Correction ---
                if isfield(posteriorParticlesAnchors{sensor}{anchor}, 'history') && ...
                   mod(step, parameters.irls_trigger_interval) == 0 && ...
                   posteriorParticlesAnchors{sensor}{anchor}.history_idx >= parameters.irls_min_history_for_correction && ...
                   strcmp(estimatedAnchors{sensor, step}{anchor}.type, 'virtual_nlos')
                    
                    [corrected_pos, corrected_cov] = correctAnchorPosition_WLS(posteriorParticlesAnchors{sensor}{anchor}, parameters);
                    
                    % Particle Rejuvenation if correction is valid
                    if ~any(isnan(corrected_pos)) && norm(mean(posteriorParticlesAnchors{sensor}{anchor}.x, 2) - corrected_pos) < 2.0
                        try
                            num_rejuvenate = round(parameters.beta * numParticles);
                            num_survive = numParticles - num_rejuvenate;
                            
                            % Keep best existing particles
                            prior_particles = posteriorParticlesAnchors{sensor}{anchor}.x;
                            weights_fusion = mvnpdf(prior_particles', corrected_pos', corrected_cov);
                            [~, sorted_idx] = sort(weights_fusion, 'descend');
                            surviving_p = prior_particles(:, sorted_idx(1:num_survive));
                            
                            % Sample new particles from WLS result
                            L_chol = chol(corrected_cov, 'lower');
                            standard_samples = trnd(parameters.nu, [2, num_rejuvenate]);
                            rejuvenated_p = corrected_pos + L_chol * standard_samples;
                            
                            posteriorParticlesAnchors{sensor}{anchor}.x = [surviving_p, rejuvenated_p];
                            posteriorParticlesAnchors{sensor}{anchor}.correction_count = posteriorParticlesAnchors{sensor}{anchor}.correction_count + 1;
                        catch 
                            % Skip fusion if matrix issues occur
                        end
                    end
                end
                
                % Update estimates
                estimatedAnchors{sensor, step}{anchor}.x = mean(posteriorParticlesAnchors{sensor}{anchor}.x, 2);
                estimatedAnchors{sensor, step}{anchor}.posteriorExistence = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence;
                
                % Accumulate weights for Agent Update
                weights(:, anchor) = predExistence * weights(:, anchor) + deadUpdate;
                weights(:, anchor) = log(weights(:, anchor));
                weights(:, anchor) = weights(:, anchor) - max(weights(:, anchor));    
            end 
        
            numEstimatedAnchors(sensor, step) = size(estimatedAnchors{sensor, step}, 2);
            
            % Sum weights across anchors for this sensor
            weightsSensors(:, sensor) = sum(weights, 2);
            weightsSensors(:, sensor) = weightsSensors(:, sensor) - max(weightsSensors(:, sensor));
            
            undetectedAnchorsIntensity(sensor) = undetectedAnchorsIntensity(sensor) * (1 - parameters.detectionProbability);
        
            % --- 3.2.7 Instantiate New Anchors ---
            for m = 1:numMeasurements
                % Calc existence probability
                post_exist = messagesNew(m) * newParticlesAnchors(m).constant / ...
                             (messagesNew(m) * newParticlesAnchors(m).constant + 1);
                
                new_id = numAnchors + m;
                posteriorParticlesAnchors{sensor}{new_id}.posteriorExistence = post_exist;
                posteriorParticlesAnchors{sensor}{new_id}.x = newParticlesAnchors(m).x;
                posteriorParticlesAnchors{sensor}{new_id}.w = post_exist / numParticles;
                
                % Init KF and History
                posteriorParticlesAnchors{sensor}{new_id}.kf_state = [measurements(1, m); 0];
                posteriorParticlesAnchors{sensor}{new_id}.kf_covariance = diag([parameters.measurementVariance, 1]); 
                posteriorParticlesAnchors{sensor}{new_id}.history = cell(1, parameters.irls_history_size);
                posteriorParticlesAnchors{sensor}{new_id}.history_idx = 0;
                posteriorParticlesAnchors{sensor}{new_id}.correction_count = 0;
          
                estimatedAnchors{sensor, step}{new_id}.x = mean(newParticlesAnchors(m).x, 2);
                estimatedAnchors{sensor, step}{new_id}.posteriorExistence = post_exist;
                estimatedAnchors{sensor, step}{new_id}.generatedAt = step;
                estimatedAnchors{sensor, step}{new_id}.uniqueID = global_anchor_id_counter;
                estimatedAnchors{sensor, step}{new_id}.type = 'virtual_nlos'; 
                global_anchor_id_counter = global_anchor_id_counter + 1;
            end
        
            % --- 3.2.8 Pruning ---
            [estimatedAnchors{sensor, step}, posteriorParticlesAnchors{sensor}, ~] = ...
                deleteUnreliableVA_v2(estimatedAnchors{sensor, step}, posteriorParticlesAnchors{sensor}, unreliabilityThreshold);
            numEstimatedAnchors(sensor, step) = size(estimatedAnchors{sensor, step}, 2);
        end 

        %% 4. Agent Update
        
        % --- 4.1 Temporal dAoA Likelihood ---
        log_weights_d_aoa = zeros(parameters.numParticles, 1);
        
        for sensor = 1:numSensors
            numAnchors_sensor = size(estimatedAnchors{sensor, step}, 2);
            for anchor = 1:numAnchors_sensor
                % Use only high-confidence anchors
                if estimatedAnchors{sensor, step}{anchor}.posteriorExistence > 0.5
        
                    d_aoa_meas = posteriorParticlesAnchors{sensor}{anchor}.kf_state(2);
                    d_aoa_var = posteriorParticlesAnchors{sensor}{anchor}.kf_covariance(2, 2);
        
                    % Skip if variance is too high (not converged)
                    if d_aoa_var > 1.0; continue; end
        
                    % Predicted dAoA (based on current agent particles)
                    anchor_pos = estimatedAnchors{sensor, step}{anchor}.x;
                    delta_pos = anchor_pos - predictedParticlesAgent(1:2, :);
                    agent_velocities = predictedParticlesAgent(3:4, :);
        
                    r_sq = sum(delta_pos.^2, 1);
                    numerator = delta_pos(1, :) .* agent_velocities(2, :) - delta_pos(2, :) .* agent_velocities(1, :);
                    predicted_d_aoa = numerator ./ r_sq;
        
                    % Log-likelihood update
                    log_lik = -0.5 * ((d_aoa_meas - predicted_d_aoa').^2 / d_aoa_var) - 0.5 * log(2*pi*d_aoa_var);
                    log_weights_d_aoa = log_weights_d_aoa + log_lik;
                end
            end
        end
        
        % --- 4.2 Fuse Weights & Resample Agent ---
        total_log_weights_aoa = sum(weightsSensors, 2);
        total_log_weights_combined = total_log_weights_aoa + log_weights_d_aoa;
        
        % Normalize
        total_log_weights_combined = total_log_weights_combined - max(total_log_weights_combined);
        weights_final = exp(total_log_weights_combined);
        weights_final = weights_final / sum(weights_final);
        
        if(parameters.known_track)
            estimatedTrajectory(:, step) = mean(predictedParticlesAgent, 2);
            posteriorParticlesAgent = predictedParticlesAgent;
        else
            estimatedTrajectory(:, step) = predictedParticlesAgent * weights_final;
            posteriorParticlesAgent = predictedParticlesAgent(:, resampleSystematic(weights_final, parameters.numParticles));
        end
        
        % Store debug data
        if(any(storing_idx == step))
            posteriorParticlesAnchorsstorage{storing_idx == step} = posteriorParticlesAnchors;
        end
        execTimePerStep(step) = toc;
        
        error_agent = calcDistance_(trueTrajectory(1:2, step), estimatedTrajectory(1:2, step));

        % % print for test
        if mod(step, 100) == 0
            if step == 100
                fprintf('\nEpoch %d ', step);
            else
                fprintf('Epoch %d ', step);
            end
            fprintf('NumAncs S1: %d ', numEstimatedAnchors(1, step));
            fprintf('NumAncs S2: %d ', numEstimatedAnchors(2, step));
            fprintf('PosErr agent: %.4f ', error_agent);
            fprintf('Time: %.4f \n', execTimePerStep(step));
        end
        
        
    end 
end