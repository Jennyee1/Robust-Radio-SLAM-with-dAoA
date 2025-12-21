% Florian Meyer, 08/01/16.


function [predictedMeans, predictedUncertainties, predictedAngle] = predictMeasurementsAngular(predictedParticles, anchorPositions, weightsAnchor)
[~, numParticles, numAnchors] = size(anchorPositions);

predictedMeans = zeros(numAnchors,1);
predictedUncertainties = zeros(numAnchors,1);
predictedAngle = zeros(numParticles,1);
for anchor = 1:numAnchors    
    % predictedAngle(:,anchor) = atan2( predictedParticles(2,:) - anchorPositions(2,:,anchor), predictedParticles(1,:) - anchorPositions(1,:,anchor) );
    predictedAngle(:,anchor) = atan2( anchorPositions(2,:,anchor)-predictedParticles(2,:), anchorPositions(1,:,anchor)-predictedParticles(1,:) );
    % sqrt((predictedParticles(1,:) - anchorPositions(1,:,anchor)).^2+(predictedParticles(2,:) - anchorPositions(2,:,anchor)).^2)';    
    predictedMeans(anchor) = predictedAngle(:,anchor).'*weightsAnchor(:,anchor)/sum(weightsAnchor(:,anchor),1);

    predictedUncertainties(anchor) = (((predictedAngle(:,anchor) - predictedMeans(anchor)).*(predictedAngle(:,anchor) - predictedMeans(anchor))).'*weightsAnchor(:,anchor))/sum(weightsAnchor(:,anchor),1);
end



end

