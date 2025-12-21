
% 
function [clutteredMeasurements] = generateClutteredAngMeasurements(trueMeasurementsCell, parameters)
% GENERATECLUTTEREDANGMEASUREMENTS Adds clutter and simulates detection probability.
%
% This function takes clean measurements and simulates a realistic sensor by:
% 1. Applying a detection probability (Pd) to true measurements.
% 2. Adding false alarms (clutter) based on a Poisson distribution.
%
% Inputs:
%   trueMeasurementsCell : {numSteps x numSensors} cell array of clean measurements.
%   parameters           : Struct containing Pd, clutter intensity, etc.
%
% Outputs:
%   clutteredMeasurements: Cell array structured like input, but with clutter.

    % range->angle(rad)
measurementVarianceAngular = parameters.measurementVariance;
detectionProbability = parameters.detectionProbability;
meanNumberOfClutter = parameters.meanNumberOfClutter;
maxRange = parameters.regionOfInterestSize;
maxAngle = 2*pi;
[numSteps, numSensors] = size(trueMeasurementsCell);
clutteredMeasurements = cell(numSteps,1);

for sensor = 1:numSensors
    for step = 1:numSteps
        trueMeasurements = trueMeasurementsCell{step,sensor};
        [~, numAnchors] = size(trueMeasurements);
        detectedMeasurements = trueMeasurements;
        
        numFalseAlarms = poissrnd(meanNumberOfClutter);
        falseAlarms = zeros(2,numFalseAlarms);
        if(~isempty(falseAlarms))
            falseAlarms(1,:) = maxAngle*rand(numFalseAlarms,1)-pi;
            falseAlarms(2,:) = measurementVarianceAngular;
        end
        if step>1
            detectionIndicator = (rand(numAnchors,1) < detectionProbability);
            detectionIndicator(1) = 1;
            detectedMeasurements = squeeze(trueMeasurements(:,detectionIndicator));
            clutteredMeasurement = [detectedMeasurements, falseAlarms];
            clutteredMeasurement = clutteredMeasurement(:,randperm(numFalseAlarms+sum(detectionIndicator)));
        else
            clutteredMeasurement = [detectedMeasurements];
            % clutteredMeasurement = clutteredMeasurement(:,randperm(sum(detectionIndicator)));
        end
        
        % clutteredMeasurement = [falseAlarms, detectedMeasurements];
        clutteredMeasurements{step,sensor} = clutteredMeasurement;
    end
end
end