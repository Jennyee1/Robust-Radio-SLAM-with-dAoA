% Author: YiJin

function [corrected_mean, corrected_covariance] = correctAnchorPosition_WLS(anchor_data, params)
% CORRECTANCHORPOSITION_WLS Solves anchor position using Weighted Least Squares (WLS).
%
% Includes keyframe selection logic to ensure geometric diversity of observations.
%
% Inputs:
%   anchor_data : Struct containing particle history and stats.
%   params      : Struct containing thresholds and configuration.

    % --- 0. Parameters ---
    max_solution_jump = params.irls_max_solution_jump;
    min_direction_consistency_ratio = params.irls_min_direction_consistency_ratio; 
    
    % --- 1. Extract History ---
    history_raw = anchor_data.history;
    history_raw = history_raw(~cellfun('isempty', history_raw));
    
    initial_mean = mean(anchor_data.x, 2);
    initial_covariance = cov(anchor_data.x');
    
    if length(history_raw) < params.irls_min_history_for_correction
        corrected_mean = initial_mean;
        corrected_covariance = initial_covariance;
        return;
    end
    
    % --- 2. Keyframe Selection ---
    if params.enable_keyframe_selection
        % Filter history for geometric diversity
        history = select_keyframes(history_raw, initial_mean, params);
    else
        history = history_raw; 
    end
    
    K = length(history);
    if K < 3 % Insufficient data
        corrected_mean = initial_mean;
        corrected_covariance = initial_covariance;
        return;
    end
    
    % --- 3. WLS Computation ---
    agent_poses_history = zeros(2, K);
    for i = 1:K, agent_poses_history(:, i) = history{i}.agent_pose(1:2); end
    
    A = zeros(K, 2);
    b = zeros(K, 1);
    weights = zeros(K, 1);
    
    for k = 1:K
        z_k = history{k}.aoa;
        x_hat_k = agent_poses_history(:, k);
        R_k = history{k}.aoa_var;
        P_x_k = history{k}.agent_cov(1:2, 1:2);
        
        sin_z = sin(z_k);     cos_z = cos(z_k);
        
        % Linearized measurement model: A*x = b
        A(k, :) = [sin_z, -cos_z];
        b(k) = x_hat_k(1) * sin_z - x_hat_k(2) * cos_z;
        
        % Uncertainty propagation (Agent Pose + Measurement Noise)
        J_pose = [sin_z, -cos_z];
        var_pose = J_pose * P_x_k * J_pose';
        
        J_aoa = x_hat_k(1) * cos_z + x_hat_k(2) * sin_z;
        var_aoa = J_aoa^2 * R_k;
        
        total_var =  var_pose + var_aoa; 
        weights(k) = 1 / max(total_var, 1e-9);
    end
    
    % Apply weights
    sqrt_w = sqrt(weights);
    A_w = A .* sqrt_w;
    b_w = b .* sqrt_w;
    
    At_W_A = A_w' * A_w;
    
    % Check conditioning
    if rcond(At_W_A) < 1e-9
        corrected_mean = initial_mean;
        corrected_covariance = initial_covariance;
        return;
    end
    
    wls_solution = At_W_A \ (A_w' * b_w);
    
    % --- 4. Validation ---
    % Check 1: Solution Jump
    if norm(wls_solution - initial_mean) > max_solution_jump
        corrected_mean = initial_mean;
        corrected_covariance = initial_covariance;
        return;
    end
    
    % Check 2: Direction Consistency (Cheirality check)
    direction_consistent_count = 0;
    for k = 1:K
        z_k = history{k}.aoa;
        x_hat_k = agent_poses_history(:, k);
        
        ray_direction_vector = [cos(z_k); sin(z_k)];
        anchor_direction_vector = wls_solution - x_hat_k;
        
        if (anchor_direction_vector' * ray_direction_vector) > 0
            direction_consistent_count = direction_consistent_count + 1;
        end
    end
    
    if (direction_consistent_count / K) < min_direction_consistency_ratio
        corrected_mean = initial_mean;
        corrected_covariance = initial_covariance;
        return;
    end
    
    % --- 5. Finalize Output ---
    corrected_mean = wls_solution;
    
    % Estimate posterior covariance using residuals
    residuals = A * corrected_mean - b;
    weighted_rss = (residuals.^2)' * weights;
    dof = K - 2;
    
    if dof > 0
        sigma2_hat = weighted_rss / dof;
        if rcond(At_W_A) > 1e-12
            corrected_covariance = sigma2_hat * inv(At_W_A);
        else
            corrected_covariance = initial_covariance;
        end
    else
        corrected_covariance = initial_covariance;
    end
end

function keyframes = select_keyframes(full_history, anchor_estimate, params)
% SELECT_KEYFRAMES Greedily selects geometrically diverse frames.
    if isempty(full_history)
        keyframes = {}; 
        return; 
    end
    
    % Always include first frame
    keyframes = {full_history{1}};
    last_keyframe_pose = full_history{1}.agent_pose(1:2);
    
    delta = anchor_estimate - last_keyframe_pose;
    last_keyframe_angle_to_anchor = atan2(delta(2), delta(1));
    
    % Iterate through history
    for i = 2:length(full_history)
        current_obs = full_history{i};
        current_pose = current_obs.agent_pose(1:2);
        
        % Check translation
        dist_from_last = norm(current_pose - last_keyframe_pose);
        
        % Check parallax angle
        delta = anchor_estimate - current_pose;
        current_angle_to_anchor = atan2(delta(2), delta(1));
        
        angle_diff = abs(current_angle_to_anchor - last_keyframe_angle_to_anchor);
        angle_diff = atan2(sin(angle_diff), cos(angle_diff)); % Normalize to [0, pi]
        
        % Selection criteria: sufficient translation AND sufficient parallax
        if dist_from_last > params.keyframe_min_distance && angle_diff > params.keyframe_min_angle
            keyframes{end+1} = current_obs; %#ok<AGROW>
            last_keyframe_pose = current_pose; 
            last_keyframe_angle_to_anchor = current_angle_to_anchor;
        end
    end
end