function [constants] = calculateConstantsUniformAngular(predictedParticlesAgent, newMeasurements, parameters)
% CALCULATECONSTANTSUNIFORMANGULAR Evaluates measurement likelihoods against a uniform spatial prior.
%
% This is used to determine if a measurement is likely to be a "New Anchor" 
% by comparing it against a background model (uniform distribution of potential anchors).

    % Define Region of Interest (ROI) for uniform particle generation
    regionOfInterestRangeX = parameters.regionOfInterestRangeX;
    regionOfInterestRangex = (regionOfInterestRangeX(2) - regionOfInterestRangeX(1)) * 3;
    regionOfInterestRangeY = parameters.regionOfInterestRangeY;
    regionOfInterestRangey = (regionOfInterestRangeY(2) - regionOfInterestRangeY(1)) * 3;
    
    regionOfInterest = regionOfInterestRangex * regionOfInterestRangey;
    
    % Upsample particles for better Monte Carlo integration density
    upSamplingFactor = parameters.upSamplingFactor;
    numParticles = parameters.numParticles * upSamplingFactor;
    predictedParticlesAgent = repmat(predictedParticlesAgent, [1, upSamplingFactor]);
    
    % Generate random particles uniformly distributed in the extended ROI
    particles = 3 * [regionOfInterestRangex; regionOfInterestRangey] .* rand(2, numParticles) - ...
                [regionOfInterestRangex; regionOfInterestRangey];
            
    constantWeight = 1 / regionOfInterest;
    numMeasurements = length(newMeasurements(1, :));
    
    % Calculate predicted angles from Agent to random particles
    predictedAngle = atan2(particles(2, :) - predictedParticlesAgent(2, :), ...
                           particles(1, :) - predictedParticlesAgent(1, :));

    % --- Gaussian Likelihood Calculation ---
    constants = zeros(numMeasurements, 1);
    for measurement = 1:numMeasurements
        % Calculate minimum angle difference
        angle_diff = normalize_angle_diff(repmat(newMeasurements(1, measurement), 1, numParticles), predictedAngle);
        
        constantLikelihood = 1 / sqrt(2 * pi * newMeasurements(2, measurement));
        likelihoods = exp((-1/2) * (angle_diff.^2) / newMeasurements(2, measurement));
        
        % Average likelihood over all random particles (Monte Carlo integration)
        constants(measurement) = sum(1/numParticles * constantLikelihood * likelihoods);
    end   
    
    % Normalize by the spatial density weight
    constants = constants / constantWeight;
end