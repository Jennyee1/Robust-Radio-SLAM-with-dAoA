% based on generateNewAnchors.m from Florian Meyer, Erik Leitinger, 20/05/17.

% Generates new potential anchors based on current measurements.
% These new anchors are candidates for newly detected features.

function [newParticlesAnchors,inputBP, parameters] = generateNewAnchorsAngular(newMeasurements,undetectedTargetsIntensity,predictedParticlesAgent,parameters)
clutterIntensity = parameters.clutterIntensity;  % Get clutter intensity from parameters
numParticles = parameters.numParticles;
numMeasurements = length(newMeasurements(1,:));
detectionProbability = parameters.detectionProbability;

inputBP = zeros(numMeasurements,1);  % Initialize input messages for BP for each measurement
newParticlesAnchors = [];  % Initialize structure for new anchor particles

if(numMeasurements)
    newParticlesAnchors = struct('x',zeros(2,numParticles),'w',zeros(numParticles,1),'posteriorExistence',0);
end

% f(z|x,a) likelihood function of measurements
constants = calculateConstantsUniformAngular(predictedParticlesAgent, newMeasurements, parameters);  % range->angle(rad)


% Loop through each measurement to generate a potential new anchor
for measurement = 1:numMeasurements

    %inputBP(measurement) = clutterIntensity + constants(measurement) * %undetectedTargetsIntensity * detectionProbability;  % "Williams style"
    inputBP(measurement) = 1 + (constants(measurement) * undetectedTargetsIntensity * detectionProbability)/clutterIntensity;  % "BP style"
    
    measurementToAnchor = newMeasurements(1,measurement);  % range->angle(rad)
    measurementVariance = newMeasurements(2,measurement);


    newParticlesAnchors(measurement).x = sampleFromLikelihoodAngular(measurementToAnchor, measurementVariance, ...
    predictedParticlesAgent, numParticles, parameters.regionOfInterestSize);

    % Store a constant factor for the new anchor's existence probability calculation
    newParticlesAnchors(measurement).constant = constants(measurement) * undetectedTargetsIntensity *...
        detectionProbability /clutterIntensity;   % detectionProbability: Pd, clutterIntensity: f_{FA}
    % Initialize particle weights for the new anchor (uniform initially)
    newParticlesAnchors(measurement).w = ones(numParticles,1)/numParticles;
    newParticlesAnchors(measurement).kf_state = [newMeasurements(1,measurement); 0]; 
    newParticlesAnchors(measurement).kf_covariance = diag([parameters.measurementVariance, 1]);

end

end