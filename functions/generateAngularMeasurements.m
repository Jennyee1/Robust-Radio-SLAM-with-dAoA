% YiJin 2025/05/15

function [measurementsCell, dataVA] = generateAngularMeasurements(targetTrajectory, dataVA, parameters)
% GENERATEANGULARMEASUREMENTS Generates noisy bearing-only measurements.
%
% This function calculates the relative angle (bearing) between the target 
% and visible anchors (virtual agents), adding Gaussian noise.
%
% Inputs:
%   targetTrajectory : [2 x numSteps] matrix of target positions.
%   dataVA           : Cell array containing sensor/anchor data.
%   parameters       : Struct containing system parameters (noise, etc.).
%
% Outputs:
%   measurementsCell : {numSteps x numSensors} cell array. Each cell contains
%                      [2 x k] matrix: row 1 is angle, row 2 is variance.
%   dataVA           : Updated sensor data (structure preserved).

    % Extract parameters
    measurementVarianceAngular = parameters.measurementVariance;
    measurementVarianceLHF     = parameters.measurementVarianceLHF;
    
    [~, numSteps] = size(targetTrajectory);
    numSensors = length(dataVA);
    
    % Initialize output cell array
    measurementsCell = cell(numSteps, numSensors);

    for sensor = 1:numSensors
        positions = dataVA{sensor}.positions;   % Anchor positions
        visibility = dataVA{sensor}.visibility; % Visibility matrix
        
        [~, numAnchors] = size(positions);

        for step = 1:numSteps
            k = 0; % Counter for visible anchors
            
            % Pre-allocate temporary measurement matrix for this step
            % Format: [Measurement Value; Measurement Variance]
            measurements = zeros(2, numAnchors); 

            for anchor = 1:numAnchors
                % Check if the specific anchor is visible at the current step
                if visibility(anchor, step)
                    k = k + 1;
                    
                    % 1. Calculate True Angle (atan2)
                    % Angle from Target to Anchor
                    dx = positions(1, anchor) - targetTrajectory(1, step);
                    dy = positions(2, anchor) - targetTrajectory(2, step);
                    trueAngle = atan2(dy, dx);
                    
                    % 2. Add Gaussian Noise
                    noiseStd = sqrt(measurementVarianceAngular);
                    noisyAngle = trueAngle + noiseStd * randn;
                    
                    % 3. Store Data
                    measurements(1, k) = noisyAngle;
                    measurements(2, k) = measurementVarianceLHF; 
                end
            end
            
            % Store only valid measurements (trim unused columns)
            if k > 0
                measurementsCell{step, sensor} = measurements(:, 1:k);
            else
                measurementsCell{step, sensor} = [];
            end
        end
    end
end
