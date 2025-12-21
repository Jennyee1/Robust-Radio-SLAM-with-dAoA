function [assocProbExisting, assocProbNew, messagelhfRatios, messagelhfRatiosNew] = dataAssociationBP(legacy, new, checkConvergence, threshold, numIterations)
% DATAASSOCIATIONBP Solves the data association problem using Belief Propagation (BP).
%
% Inputs:
%   legacy: [M+1 x N] Matrix of likelihoods for legacy anchors (Row 1 = Missed Detection).
%   new:    [M x 1] Likelihoods for measurements being new features.
%
% Outputs:
%   assocProbExisting: Marginal posterior probabilities for legacy anchors.
%   assocProbNew:      Marginal posterior probabilities for new features.
%   messagelhfRatios:  Messages passed back to legacy anchors (used for state update).
%   messagelhfRatiosNew: Messages passed back to new feature initialization.

    [m, n] = size(legacy);    
    m = m - 1;   % m = numMeasurements, n = numAnchors
    
    assocProbNew = ones(m, 1);
    assocProbExisting = ones(m + 1, n);
    messagelhfRatios = ones(m, n);  
    messagelhfRatiosNew = ones(m, 1); 
    
    if(n == 0 || m == 0)  
        return;
    end
    if(isempty(new))  
        new = 1;
    end
    
    om = ones(1, m);  
    on = ones(1, n);  
    
    % Initialize messages from measurement nodes (b) to anchor nodes (a)
    muba = ones(m, n); 
    
    % --- BP Iteration Loop ---
    for iteration = 1:numIterations
        mubaOld = muba;
      
        % 1. Compute messages a -> b (muab)
        prodfact = muba .* legacy(2:end, :);          % Product term
        sumprod = legacy(1, :) + sum(prodfact, 1);    % Sum over all measurements (col sum)
        normalization = (sumprod(om, :) - prodfact);  % Exclude current term
        
        normalization(normalization == 0) = eps;      % Avoid division by zero
        muab = legacy(2:end, :) ./ normalization;     % Update muab
      
        % 2. Compute messages b -> a (muba)
        summuab = new + sum(muab, 2);                 % Sum over all anchors + new feature hypothesis
        normalization = summuab(:, on) - muab;        % Exclude current term
        
        normalization(normalization == 0) = eps;
        muba = 1 ./ normalization;
        
        if any(isnan(muba(:))) || any(isnan(muab(:))) || any(isnan(summuab(:)))
            flag_nan = 1;
        end
      
        % Check for convergence
        if(mod(iteration, checkConvergence) == 0)
            distance = max(max(abs(log(muba ./ mubaOld))));
            if(distance < threshold)
                break
            end
        end
    end
    
    % --- Compute Marginals (Posteriors) ---
    assocProbExisting(1, :) = legacy(1, :);
    assocProbExisting(2:end, :) = legacy(2:end, :) .* muba;
    
    % Normalize existing anchor probabilities
    for target = 1:n 
        assocProbExisting(:, target) = assocProbExisting(:, target) / sum(assocProbExisting(:, target));
    end
    
    % Compute messages/probabilities for output
    messagelhfRatios = muba;  
    assocProbNew = new ./ summuab;  
    
    % Normalize New Feature messages
    messagelhfRatiosNew = [ones(m, 1), muab];  
    messagelhfRatiosNew = messagelhfRatiosNew ./ repmat(sum(messagelhfRatiosNew, 2), [1, n + 1]);
    messagelhfRatiosNew = messagelhfRatiosNew(:, 1);
    
end