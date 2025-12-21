function [associationProbabilities, associationProbabilitiesNew, messagelegacy, messagesNew] = calculateAssociationProbabilitiesGAAngular(...
    measurements, predictedMeasurements, predictedUncertaintys, weightsAnchor, newInputBP, parameters)
% CALCULATEASSOCIATIONPROBABILITIESGAANGULAR Computes data association probabilities via Belief Propagation.
%
% This function determines the probability that a measurement originates from:
% 1. An existing anchor (legacy),
% 2. A new anchor, or
% 3. Clutter (false alarm).

    detectionProbability = parameters.detectionProbability;
    clutterIntensity = parameters.clutterIntensity;
    
    [~, numMeasurements] = size(measurements);
    numAnchors = length(predictedMeasurements);
    
    % Sum of weights represents the predicted existence probability of each anchor
    predictedExistence = sum(weightsAnchor, 1);
    
    % Initialize BP input messages matrix (Equation 25 in reference paper)
    % Rows: Measurements (Row 1 = Missed Detection)
    % Cols: Anchors
    inputBP = zeros(numMeasurements + 1, numAnchors);
       
    % Row 1: Hypothesis that anchor k is NOT detected
    inputBP(1, :) = (1 - detectionProbability);
    
    for anchor = 1:numAnchors
        for measurement = 1:numMeasurements
            % Combined uncertainty (Innovation covariance)
            predictedUncertaintyTmp = predictedUncertaintys(anchor) + measurements(2, measurement);
                % --- Gaussian Likelihood ---
                factor = 1 / sqrt(2 * pi * predictedUncertaintyTmp) * detectionProbability / clutterIntensity;
                diff_ang = normalize_angle_diff(measurements(1, measurement), predictedMeasurements(anchor));
                inputBP(measurement + 1, anchor) = factor * exp(-1 / (2 * predictedUncertaintyTmp) * (diff_ang).^2);
        end
    
        % Sanity check for complex numbers
        if ~isreal(inputBP) && any(imag(inputBP(:)) ~= 0)
            warning('Complex values detected in inputBP.');
        end
        
        % Combine existence probability with association likelihoods (Equation 25)
        inputBP(:, anchor) = getInputBP(predictedExistence(anchor), inputBP(:, anchor));
        
        if any(isnan(inputBP(:)))
            warning('NaN values detected in inputBP.');
        end
    end
    
    % Execute Loopy Belief Propagation
    % associationProbabilities: Probabilities for legacy anchors (cols sum to 1)
    % associationProbabilitiesNew: Probabilities that measurements are new features
    [associationProbabilities, associationProbabilitiesNew, messagelegacy, messagesNew] = ...
        dataAssociationBP(inputBP, newInputBP, 30, 10^(-6), 10^6);

    if ~isreal(messagesNew) && any(imag(messagesNew(:)) ~= 0)
        warning('Complex values detected in messagesNew.');
    end
    if any(isnan(messagesNew))
        warning('NaN values detected in messagesNew.');
    end
end