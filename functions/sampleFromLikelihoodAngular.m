

function [ samples ] = sampleFromLikelihoodAngular(measurementToAnchor, measurementVariance, agentPosition, numParticles, maxRange)
samples = zeros(2,numParticles);

% range->angle(rad)
% r = measurementToAnchor + sqrt(measurementVariance)*randn(1,numParticles);
% phi = 2*pi*rand(1,numParticles);
AoA = measurementToAnchor + sqrt(measurementVariance)*randn(1,numParticles);
range = maxRange*rand(1,numParticles);

AoA_agent_anchor = AoA;

% samples(1,:) = agentPosition(1,:) + r.*cos(phi);
% samples(2,:)= agentPosition(2,:) + r.*sin(phi);
samples(1,:) = agentPosition(1,:) + range.*cos(AoA_agent_anchor);
samples(2,:)= agentPosition(2,:) + range.*sin(AoA_agent_anchor);

end

