function [estimatedAnchors, posteriorParticlesAnchors, unreliableAnchors_indices] = deleteUnreliableVA_v2(estimatedAnchors, posteriorParticlesAnchors, unreliabilityThreshold)
% DELETEUNRELIABLEVA_V2 Prunes anchors with low existence probability.
%
% Handles empty cells safely to ensure robustness against previously deleted or invalid entries.
%
% Inputs:
%   unreliabilityThreshold: Probability threshold below which an anchor is deleted.

    numAnchors = length(posteriorParticlesAnchors); 
    unreliableAnchors_indices = [];
    
    for anchor_idx = 1:numAnchors  
        
        % Skip physical LOS anchors (index 1 is typically the LOS anchor)
        current_anchor_particles = posteriorParticlesAnchors{anchor_idx};
        current_anchor_type = estimatedAnchors{anchor_idx}.type;
        
        if strcmp(current_anchor_type, 'physical_los')
            continue; 
        end
        
        % Check for empty cells or invalid existence probabilities
        is_empty = isempty(current_anchor_particles);
        
        if is_empty || isnan(current_anchor_particles.posteriorExistence)
            % Mark for deletion if empty or invalid
            unreliableAnchors_indices = [unreliableAnchors_indices, anchor_idx];
        else
            % Check existence probability
            priorExistence = current_anchor_particles.posteriorExistence;
            if(priorExistence < unreliabilityThreshold) 
                unreliableAnchors_indices = [unreliableAnchors_indices, anchor_idx];
            end
        end
    end
    
    % Perform deletion (maintaining array compactness)
    if ~isempty(unreliableAnchors_indices)
        posteriorParticlesAnchors(unreliableAnchors_indices) = [];
        estimatedAnchors(unreliableAnchors_indices) = [];
    end
end