

% Initialization function for anchors (PAs and initial VAs).
% Sets up initial particle distributions and estimated states.

function [ estimatedAnchors, posteriorParticlesAnchors ] = initAnchors( parameters, dataVA, numSteps, numSensors )

numParticles = parameters.numParticles;
posteriorParticlesAnchors = cell(numSensors,1);
estimatedAnchors = cell(numSensors,numSteps);
for sensor = 1:numSensors
    anchorPositions = dataVA{sensor}.positions(:,parameters.priorKnownAnchors{sensor});
    [~, numAnchors] = size(anchorPositions);
    posteriorParticlesAnchors{sensor} = cell(1,numAnchors);
    estimatedAnchors{sensor,1} = cell(1,numAnchors);
end 

priorCovarianceAnchor = parameters.priorCovarianceAnchor;  
agent_initial_pos = parameters.priorMean(1:2); %

for sensor = 1:numSensors
    anchorPositions = dataVA{sensor}.positions(:,parameters.priorKnownAnchors{sensor});
    [~, numAnchors] = size(anchorPositions); 
    for anchor = 1:numAnchors
       
        anchor_pos = dataVA{sensor}.positions(:, anchor);
        delta_y = anchor_pos(2) - agent_initial_pos(2);
        delta_x = anchor_pos(1) - agent_initial_pos(1);
        initial_aoa = atan2(delta_y, delta_x);

        posteriorParticlesAnchors{sensor}{anchor}.x = zeros(2,numParticles); 
        posteriorParticlesAnchors{sensor}{anchor}.w = zeros(numParticles,1);  
        posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence = 1;
        posteriorParticlesAnchors{sensor}{anchor}.w(:) = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence/numParticles*ones(numParticles,1);
        
        posteriorParticlesAnchors{sensor}{anchor}.kf_state = [initial_aoa; 0]; 
        posteriorParticlesAnchors{sensor}{anchor}.kf_covariance = diag([parameters.measurementVariance, 1]); 

        if(anchor == 1)
            posteriorParticlesAnchors{sensor}{anchor}.x =  mvnrnd(anchorPositions(:,anchor),priorCovarianceAnchor,numParticles)';

            estimatedAnchors{sensor,1}{anchor}.x = anchorPositions(:,anchor);
            estimatedAnchors{sensor,1}{anchor}.posteriorExistence = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence;
            estimatedAnchors{sensor,1}{anchor}.generatedAt = 1;
            estimatedAnchors{sensor,1}{anchor}.notObservationCount = 0;
            estimatedAnchors{sensor,1}{anchor}.uniqueID = sensor;
            estimatedAnchors{sensor,1}{anchor}.type = 'physical_los'; 
        else
            
            anchorPositions(:,anchor) = mvnrnd(anchorPositions(:,anchor),priorCovarianceAnchor,1);
            posteriorParticlesAnchors{sensor}{anchor}.x =  mvnrnd(anchorPositions(:,anchor),priorCovarianceAnchor,numParticles)';
            estimatedAnchors{sensor,1}{anchor}.x = anchorPositions(:,anchor);
            estimatedAnchors{sensor,1}{anchor}.posteriorExistence  = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence;
            estimatedAnchors{sensor,1}{anchor}.generatedAt = 1;
            estimatedAnchors{sensor,1}{anchor}.notObservationCount = 0;
            estimatedAnchors{sensor,1}{anchor}.uniqueID = sensor;
            estimatedAnchors{sensor,1}{anchor}.type = 'physical_los'; 
        end
    end
end

end

