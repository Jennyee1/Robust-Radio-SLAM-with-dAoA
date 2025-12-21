% YiJin 2025/7/15

function [predictedMeans, predictedUncertainties, predictedAngle] = predictMeasurementsAngular_t(predictedParticles, anchorPositions, weightsAnchor)
[~, numParticles, numAnchors] = size(anchorPositions);


predictedMeans = zeros(numAnchors,1);
predictedUncertainties = zeros(numAnchors,1);
predictedAngle = zeros(numParticles,numAnchors);
for anchor = 1:numAnchors    

    current_angles = atan2( anchorPositions(2,:,anchor)-predictedParticles(2,:), anchorPositions(1,:,anchor)-predictedParticles(1,:) )';
    predictedAngle(:,anchor) = current_angles;  
    current_weights = weightsAnchor(:,anchor);
    sumOfWeights = sum(current_weights);
    if sumOfWeights < 1e-15
        sumOfWeights = 1e-15;
    end
    
    mean_angle = (current_angles' * current_weights) / sumOfWeights;
    predictedMeans(anchor) = mean_angle;
    mean_of_squares = ((current_angles.^2)' * current_weights) / sumOfWeights;
        
    % E[X^2] - (E[X])^2
    predictedUncertainties(anchor) = mean_of_squares - mean_angle^2;
    if predictedUncertainties(anchor)< 0
        flag = 1;
    end
    if any(isnan(predictedMeans(:))) || any(isnan(predictedUncertainties(:)))
        flag_nan = 1;
    end
end


end

