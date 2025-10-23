% 

function [ estimatedTrajectory, estimatedAnchors, posteriorParticlesAnchorsstorage, numEstimatedAnchors,posteriorParticlesAgentstorage ] =  BPbasedSLAM_dAoA( dataVA, clutteredMeasurements, parameters, trueTrajectory )
mycolors = [ ...
0.66,0.00,0.00; ...
0.00,0.30,0.70; ...
0.92,0.75,0.33; ...
0.00,0.00,0.00; ...
];
if parameters.Isplotrt
    h_paranc = cell(2,50);
    h_estimatedAnc = cell(2,50);
    for sensor = 1:3
        for anchor = 1:20
            h_paranc{sensor, anchor} = plot(nan, nan, ...  %
            'color', mycolors(sensor, :), ...
            'marker', 'o', ...
            'MarkerSize', 3, ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('传感器%d-锚点%d', sensor, anchor));
            h_estimatedAnc{sensor, anchor} = plot(nan, nan, 'color', 'k', ...
            'marker', '+', 'MarkerSize', 5,  'LineWidth', 2.5);
        end
    end
    % legend('Location', 'best'); 
end

global_anchor_id_counter = 500; % 
[numSteps,numSensors] = size(clutteredMeasurements);
numSteps = min(numSteps,parameters.maxSteps);
numParticles = parameters.numParticles;
detectionProbability = parameters.detectionProbability;  % Pd
priorMean = parameters.priorMean; 
survivalProbability = parameters.survivalProbability;    % Ps
undetectedAnchorsIntensity = parameters.undetectedAnchorsIntensity*ones(numSensors,1); % mu*f_{n,n}
birthIntensity = parameters.birthIntensity;   % mu*f_{b,n}
clutterIntensity = parameters.clutterIntensity; % mu*f_{FA}
unreliabilityThreshold = parameters.unreliabilityThreshold; % P_{prun}
execTimePerStep = zeros(numSteps,1);
known_track = parameters.known_track;




parameters.irls_history_size = 200; % 
parameters.irls_trigger_interval = 10; %
parameters.irls_min_history_for_correction = 10; 
parameters.irls_min_baseline_std = 0.5;
parameters.irls_max_solution_jump = 2.0; %
parameters.irls_min_direction_consistency_ratio = 0.7; %

parameters.enable_keyframe_selection = true; % 
parameters.keyframe_min_distance = 1.0;     % 
parameters.keyframe_min_angle = 5 * (pi/180); % 
beta = 0.25; % 
sigma_a_sq = 1.5; % 

% allocate memory
estimatedTrajectory = zeros(4,numSteps);
numEstimatedAnchors = zeros(2,numSteps);
% storing_idx = 30:30:numSteps;
storing_idx = 50:50:numSteps;
posteriorParticlesAnchorsstorage = cell(1,length(storing_idx));
posteriorParticlesAgentstorage = cell(1,length(storing_idx));

% figure(1);

% --- initial state vectors ---
if(known_track)
    posteriorParticlesAgent = repmat([trueTrajectory(:,1);0;0],1,numParticles);
else 
    posteriorParticlesAgent(1:2,:) = drawSamplesUniformlyCirc(priorMean(1:2), parameters.UniformRadius_pos ,parameters.numParticles);
    posteriorParticlesAgent(3:4,:) = repmat(priorMean(3:4),1,parameters.numParticles) + 2*parameters.UniformRadius_vel * rand( 2, parameters.numParticles ) - parameters.UniformRadius_vel;
end

% Estimate the initial trajectory as the mean of the initial agent particles
estimatedTrajectory(:,1) = mean(posteriorParticlesAgent,2);

% Initialize anchors (PAs and VAs) using the initAnchors function
if parameters.priorKnownAnchors{1}
    [ estimatedAnchors, posteriorParticlesAnchors ] =  initAnchors( parameters, dataVA, numSteps, numSensors );
    for sensor = 1:numSensors
        numKnownAnchors = size(estimatedAnchors{sensor, 1}, 2);
        for anchor = 1:numKnownAnchors
            posteriorParticlesAnchors{sensor}{anchor}.history = cell(1, parameters.irls_history_size);
            posteriorParticlesAnchors{sensor}{anchor}.history_idx = 0;
            posteriorParticlesAnchors{sensor}{anchor}.correction_count = 0;
        end
    end
end
anchor_los_pos = zeros(2, numSensors);
for sensor = 1:numSensors
numEstimatedAnchors(sensor, 1) = size(estimatedAnchors{sensor,1},2);
anchor_los_pos(:, sensor) = dataVA{sensor}.positions(:,parameters.priorKnownAnchors{sensor});
end

if parameters.Isplotrt
    figure(1);
end
%% main loop
for step = 2:numSteps 
    tic
    if(known_track)
        predictedParticlesAgent = repmat([trueTrajectory(:,step);0;0],1,numParticles);
    else
        predictedParticlesAgent = performPrediction( posteriorParticlesAgent, parameters ); % f(x_n|x_(n-1))
    end
    weightsSensors = nan(numParticles,numSensors);
    % 
    for sensor = 1:numSensors
    
        estimatedAnchors{sensor,step} = estimatedAnchors{sensor,step-1};
        measurements = clutteredMeasurements{step,sensor};
        numMeasurements = size(measurements,2);
        undetectedAnchorsIntensity(sensor) = undetectedAnchorsIntensity(sensor) * survivalProbability + birthIntensity;

        [predictedParticlesAnchors, weightsAnchor] = predictAnchors( posteriorParticlesAnchors{sensor}, parameters );
        
        % Create new anchors (potential new features) based on current measurements
        [newParticlesAnchors,newInputBP, parameters] = generateNewAnchorsAngular(measurements, undetectedAnchorsIntensity(sensor) , predictedParticlesAgent, parameters);
        
        % Predict measurements from the predicted agent particles to the predicted anchor particles
        [predictedMeasurements, predictedUncertainties, predictedAngle] = predictMeasurementsAngular_t(predictedParticlesAgent, predictedParticlesAnchors, weightsAnchor);
        
        % Compute association probabilities using Belief Propagation (BP) for Data Association (DA)
        % This is a core part of the BP-SLAM algorithm, handling measurement-to-feature association 
        [associationProbabilities, associationProbabilitiesNew, messagelhfRatios, messagesNew] = calculateAssociationProbabilitiesGAAngular(...
        measurements, predictedMeasurements, predictedUncertainties, weightsAnchor, newInputBP, parameters);
        agent_mean_pos = mean(predictedParticlesAgent(1:2,:), 2);
        agent_pose_cov_current = cov(predictedParticlesAgent(1:2,:)');

        % PF anchor state update
        numAnchors = size(predictedParticlesAnchors,3);
        weights = zeros(numParticles,numAnchors);
        for anchor = 1:numAnchors
            weights(:,anchor) = repmat((1-detectionProbability),numParticles,1);
            anchor_type = estimatedAnchors{sensor, step}{anchor}.type;
            constraint_weights = ones(numParticles, 1); 
            % % 
            % if ~strcmp(anchor_type, 'physical_los') && ~isempty(anchor_los_pos(:, sensor))
            %     % 1. 
            %     nlos_anchor_pos = mean(predictedParticlesAnchors(:,:,anchor), 2);
            %     margin = 0.4;
            %     epsilon = 1e-3;
            % 
            %     % 2. 
            %     %    predictedParticlesAgent(1:2,:) -> 2 x numParticles
            %     %    nlos_anchor_pos/anchor_los_pos -> 2 x 1
            %     dist_to_nlos_per_particle = vecnorm(predictedParticlesAnchors(:,:,anchor) - agent_mean_pos);
            %     dist_to_los_per_particle = vecnorm(anchor_los_pos(:, sensor) - agent_mean_pos)+0.1;
            % 
            %     % 3. 
            %     violation_dist_per_particle = dist_to_los_per_particle - dist_to_nlos_per_particle;
            %     dist_diff_per_particle = dist_to_nlos_per_particle - dist_to_los_per_particle;
            %     satisfied_idx = dist_diff_per_particle > margin;
            %     violated_idx = dist_diff_per_particle < -margin;
            %     transition_idx = ~satisfied_idx & ~violated_idx;
            %     constraint_weights(violated_idx) = epsilon;
            %     constraint_weights(transition_idx) = (dist_diff_per_particle(transition_idx) + margin) / (2 * margin);
            % 
            %     % 4. 
            %     % beta = 100 ./ dist_to_los_per_particle; % 
            %     % constraint_weights = 1 ./ (1 + exp(beta .* violation_dist_per_particle));
            %     % constraint_weights = constraint_weights'; % numParticles x 1 
            % end
            for measurement = 1:numMeasurements
                measurementVariance = measurements(2,measurement);
                
                if parameters.useTDis
                    likelihood = student_t_pdf(measurements(1,measurement), predictedAngle(:,anchor), measurementVariance, parameters.studentT_nu);
                    %  g()  P(z|a) * Pd / f_FA
                    weights(:,anchor) = weights(:,anchor) + likelihood * detectionProbability*messagelhfRatios(measurement,anchor) / clutterIntensity;
                else
                    %
                    %factor = 1/sqrt(2*pi*measurementVariance)*detectionProbability; % "Williams style"
                    factor = 1/sqrt(2*pi*measurementVariance)*detectionProbability/clutterIntensity; % "BP style"
                    % weights(:,anchor) = weights(:,anchor) + factor*messagelhfRatios(measurement,anchor)*...
                    % exp(-1/(2*measurementVariance)*(measurements(1,measurement)-predictedAngle(:,anchor)).^2);
                    %

                    % 
                    likelihood_term = exp(-1/(2*measurementVariance)*(measurements(1,measurement)-predictedAngle(:,anchor)).^2);
                    constrained_likelihood_term = likelihood_term .* constraint_weights;
                    % weights(:,anchor) = weights(:,anchor) + factor*messagelhfRatios(measurement,anchor)*...
                    %     exp(-1/(2*measurementVariance)*(measurements(1,measurement)-predictedAngle(:,anchor)).^2);
                    weights(:,anchor) = weights(:,anchor) + factor*messagelhfRatios(measurement,anchor) * constrained_likelihood_term;
                end
            end
            
            %%
            % --- 1.  ---
            kf_state_prev = posteriorParticlesAnchors{sensor}{anchor}.kf_state;
            kf_cov_prev = posteriorParticlesAnchors{sensor}{anchor}.kf_covariance;
            
            % (F, Q, H, R)
            dt = parameters.scanTime;
            F = [1 dt; 0 1];
            
            Q = [dt^4/4, dt^3/2; dt^3/2, dt^2] * sigma_a_sq;
            H = [1 0];
            R = parameters.measurementVariance; %

            kf_state_pred = F * kf_state_prev;
            kf_cov_pred = F * kf_cov_prev * F' + Q;
            
            % --- 2.  ---
            prob_not_detected = associationProbabilities(1, anchor);
            association_probs_for_measurements = associationProbabilities(2:end, anchor);
            [max_assoc_prob, measurement_idx] = max(association_probs_for_measurements);

            is_associated = ~isempty(max_assoc_prob) && (max_assoc_prob > prob_not_detected) && (max_assoc_prob > 0.3);
            if is_associated 
                z = measurements(1, measurement_idx); 
                y = z - H * kf_state_pred;
                S = H * kf_cov_pred * H' + R;
                K = kf_cov_pred * H' / S;
                kf_state_updated = kf_state_pred + K * y;
                kf_cov_updated = (eye(2) - K * H) * kf_cov_pred;
                if isfield(posteriorParticlesAnchors{sensor}{anchor}, 'history')
                    history_idx_current = mod(posteriorParticlesAnchors{sensor}{anchor}.history_idx, parameters.irls_history_size) + 1;
                    observation_data.aoa = z;
                    observation_data.aoa_var = measurements(2, measurement_idx);
                    observation_data.agent_pose = mean(predictedParticlesAgent, 2); 
                    observation_data.agent_cov = agent_pose_cov_current;
                    
                    posteriorParticlesAnchors{sensor}{anchor}.history{history_idx_current} = observation_data;
                    posteriorParticlesAnchors{sensor}{anchor}.history_idx = history_idx_current;
                end
            else
                kf_state_updated = kf_state_pred;
                kf_cov_updated = kf_cov_pred;
            end
            
            % --- 3.
            posteriorParticlesAnchors{sensor}{anchor}.kf_state = kf_state_updated;
            posteriorParticlesAnchors{sensor}{anchor}.kf_covariance = kf_cov_updated;

            % legacy PF 
            predictedExistence = sum(weightsAnchor(:,anchor)); 
            aliveUpdate = sum(predictedExistence*1/numParticles*weights(:,anchor));  % 
            deadUpdate = 1 - predictedExistence;
            posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence = aliveUpdate/(aliveUpdate+deadUpdate);
            
            % Resample particles based on weights to update anchor position distribution
            idx_resampling = resampleSystematic(weights(:,anchor)/sum(weights(:,anchor)),numParticles);
            posteriorParticlesAnchors{sensor}{anchor}.x = predictedParticlesAnchors(:,idx_resampling(1:numParticles),anchor);
            posteriorParticlesAnchors{sensor}{anchor}.w = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence/numParticles*ones(numParticles,1);
            
            % --- IRLS ---
            if isfield(posteriorParticlesAnchors{sensor}{anchor}, 'history') && ...
               mod(step, parameters.irls_trigger_interval) == 0 && ...
               posteriorParticlesAnchors{sensor}{anchor}.history_idx >= parameters.irls_min_history_for_correction && ...
               strcmp(estimatedAnchors{sensor, step}{anchor}.type, 'virtual_nlos')

                % [corrected_pos, corrected_cov] = correctAnchorPosition_WLS(posteriorParticlesAnchors{sensor}{anchor}, numParticles);
                [corrected_pos, corrected_cov] = correctAnchorPosition_WLS(posteriorParticlesAnchors{sensor}{anchor}, parameters);

                % --- (Particle Rejuvenation) ---
                if ~any(isnan(corrected_pos)) && norm(mean(posteriorParticlesAnchors{sensor}{anchor}.x, 2) - corrected_pos) < 2.0
                    
                    try
                        
                        num_rejuvenate = round(beta * numParticles);
                        num_survive = numParticles - num_rejuvenate;
                        
                        prior_particles = posteriorParticlesAnchors{sensor}{anchor}.x;
                        weights2 = mvnpdf(prior_particles', corrected_pos', corrected_cov);
            
                        [~, sorted_indices] = sort(weights2, 'descend');
                        survivor_indices = sorted_indices(1:num_survive);
                        surviving_particles = prior_particles(:, survivor_indices);
                        
                        L = chol(corrected_cov, 'lower');
                        rejuvenated_particles = corrected_pos + L * randn(2, num_rejuvenate);
                        
                        new_particles = [surviving_particles, rejuvenated_particles];
                        
                        posteriorParticlesAnchors{sensor}{anchor}.x = new_particles;
                        posteriorParticlesAnchors{sensor}{anchor}.correction_count = posteriorParticlesAnchors{sensor}{anchor}.correction_count + 1;
                        
                    catch ME
                        % disp(['Particle Rejuvenation failed: ', ME.message, '. Skipping fusion.']);
                    end
                end
            end

            % Update estimated anchor position (mean of particles) and existence probability
            estimatedAnchors{sensor,step}{anchor}.x = mean(posteriorParticlesAnchors{sensor}{anchor}.x,2);
            estimatedAnchors{sensor,step}{anchor}.posteriorExistence = posteriorParticlesAnchors{sensor}{anchor}.posteriorExistence;
            
            weights(:,anchor) = predictedExistence*weights(:,anchor) + deadUpdate;
            weights(:,anchor) = log(weights(:,anchor));
            weights(:,anchor) = weights(:,anchor) - max(weights(:,anchor));    
            if parameters.Isplotrt
                % h_paranc{sensor, anchor} = plot(posteriorParticlesAnchors{sensor}{anchor}.x(1,:), posteriorParticlesAnchors{sensor}{anchor}.x(2,:), 'color', mycolors(sensor,:),'marker','o', 'MarkerSize', 3,'LineWidth', 1.5);
                set( h_paranc{sensor, anchor}, 'XData', posteriorParticlesAnchors{sensor}{anchor}.x(1,:), 'YData', posteriorParticlesAnchors{sensor}{anchor}.x(2,:) );
                set( h_estimatedAnc{sensor, anchor}, 'XData', estimatedAnchors{sensor,step}{anchor}.x(1,:), 'YData', estimatedAnchors{sensor,step}{anchor}.x(2,:) );
            end
        end
    
        % legacy PF number
        numEstimatedAnchors(sensor, step) = size(estimatedAnchors{sensor,step},2);
        weightsSensors(:,sensor) = sum(weights,2);
        weightsSensors(:,sensor) = weightsSensors(:,sensor) - max(weightsSensors(:,sensor));
        
        undetectedAnchorsIntensity(sensor) = undetectedAnchorsIntensity(sensor) * (1-parameters.detectionProbability);
    
        % Update new anchors (newly generated anchors)
        for measurement = 1:numMeasurements
            % Calculate posterior existence probability
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.posteriorExistence = ...
            messagesNew(measurement)*newParticlesAnchors(measurement).constant/(messagesNew(measurement)*newParticlesAnchors(measurement).constant + 1);
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.x = newParticlesAnchors(measurement).x;
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.w = posteriorParticlesAnchors{sensor}{numAnchors+measurement}.posteriorExistence/numParticles;
            
            initial_aoa = measurements(1, measurement);
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.kf_state = [initial_aoa; 0];
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.kf_covariance = diag([parameters.measurementVariance, 1]); 
            
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.history = cell(1, parameters.irls_history_size);
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.history_idx = 0;
            posteriorParticlesAnchors{sensor}{numAnchors+measurement}.correction_count = 0;
      
            estimatedAnchors{sensor,step}{numAnchors+measurement}.x = mean(newParticlesAnchors(measurement).x,2);
            estimatedAnchors{sensor,step}{numAnchors+measurement}.posteriorExistence = posteriorParticlesAnchors{sensor}{numAnchors+measurement}.posteriorExistence;
            estimatedAnchors{sensor,step}{numAnchors+measurement}.generatedAt = step;
            estimatedAnchors{sensor,step}{numAnchors+measurement}.notObservationCount = 0;
            estimatedAnchors{sensor,step}{numAnchors+measurement}.uniqueID = global_anchor_id_counter;
            global_anchor_id_counter = global_anchor_id_counter + 1;
            estimatedAnchors{sensor,step}{numAnchors+measurement}.type = 'virtual_nlos'; 
            if parameters.Isplotrt
                % h_paranc{sensor, anchor} = plot(posteriorParticlesAnchors{sensor}{anchor}.x(1,:), posteriorParticlesAnchors{sensor}{anchor}.x(2,:), 'color', mycolors(sensor,:),'marker','o', 'MarkerSize', 3,'LineWidth', 1.5);
                set( h_paranc{sensor, numAnchors+measurement}, 'XData', posteriorParticlesAnchors{sensor}{numAnchors+measurement}.x(1,:), 'YData', posteriorParticlesAnchors{sensor}{numAnchors+measurement}.x(2,:) );
                set( h_estimatedAnc{sensor, numAnchors+measurement}, 'XData', estimatedAnchors{sensor,step}{numAnchors+measurement}.x(1,:), 'YData', estimatedAnchors{sensor,step}{numAnchors+measurement}.x(2,:) );
            end
        end
    
        % Delete 
        % [estimatedAnchors{sensor,step}, posteriorParticlesAnchors{sensor}] = deleteUnreliableVA( estimatedAnchors{sensor,step}, posteriorParticlesAnchors{sensor}, unreliabilityThreshold );
        [estimatedAnchors{sensor,step}, posteriorParticlesAnchors{sensor}, unreliableAnchors_indices] = deleteUnreliableVA_v2( estimatedAnchors{sensor,step}, posteriorParticlesAnchors{sensor}, unreliabilityThreshold );
        numEstimatedAnchors(sensor, step) = size(estimatedAnchors{sensor,step},2);
        if parameters.Isplotrt
            if ~isempty(unreliableAnchors_indices)
                for i =  1: length(unreliableAnchors_indices)
                    set( h_paranc{sensor, unreliableAnchors_indices(i)}, 'XData', [], 'YData', [] );
                    set( h_estimatedAnc{sensor, unreliableAnchors_indices(i)}, 'XData', [], 'YData', [] );
                end
            end
        end
    end

    %%
    log_weights_d_aoa = zeros(parameters.numParticles, 1);
    for sensor = 1:numSensors
        numAnchors_sensor = size(estimatedAnchors{sensor, step}, 2);
        for anchor = 1:numAnchors_sensor
            if estimatedAnchors{sensor,step}{anchor}.posteriorExistence > 0.5
    
                d_aoa_meas = posteriorParticlesAnchors{sensor}{anchor}.kf_state(2);
                d_aoa_var = posteriorParticlesAnchors{sensor}{anchor}.kf_covariance(2,2);
    
                if d_aoa_var > 1.0
                    continue;
                end
    
                anchor_pos = estimatedAnchors{sensor,step}{anchor}.x;
    
                delta_pos = anchor_pos - predictedParticlesAgent(1:2, :);
                agent_velocities = predictedParticlesAgent(3:4, :);
    
                r_sq = sum(delta_pos.^2, 1);
    
                % numerator: (x-ax)*vy - (y-ay)*vx
                numerator = delta_pos(1, :) .* agent_velocities(2, :) - delta_pos(2, :) .* agent_velocities(1, :);
                predicted_d_aoa = numerator ./ r_sq;
                log_likelihood_d_aoa = -0.5 * ((d_aoa_meas - predicted_d_aoa').^2 / d_aoa_var) - 0.5 * log(2*pi*d_aoa_var);
                log_weights_d_aoa = log_weights_d_aoa + log_likelihood_d_aoa;
            end
        end
    end
    
    total_log_weights_aoa = sum(weightsSensors, 2);
    
    total_log_weights_combined = total_log_weights_aoa + log_weights_d_aoa;
    
    total_log_weights_combined = total_log_weights_combined - max(total_log_weights_combined);
    weights_final = exp(total_log_weights_combined);
    weights_final = weights_final / sum(weights_final);


    % % MODIFIED CODE ADAPTIVE RESAMPLING
    % weightsSensors = sum(weightsSensors,2);
    % weightsSensors = weightsSensors - max(weightsSensors);
    % weightsSensors = exp(weightsSensors);
    % weightsSensors = weightsSensors/sum(weightsSensors);

    if(known_track)
        estimatedTrajectory(:,step) = mean(predictedParticlesAgent,2);
        posteriorParticlesAgent = predictedParticlesAgent;
    else
        estimatedTrajectory(:,step) = predictedParticlesAgent * weights_final;
        posteriorParticlesAgent = predictedParticlesAgent(:, resampleSystematic(weights_final, parameters.numParticles));
        % estimatedTrajectory(:,step) = predictedParticlesAgent*weightsSensors;
        % posteriorParticlesAgent = predictedParticlesAgent(:,resampleSystematic(weightsSensors,numParticles));
    end


    if(any(storing_idx == step))
        posteriorParticlesAnchorsstorage{storing_idx == step} = posteriorParticlesAnchors;
        % posteriorParticlesAgentstorage{storing_idx == step} = predictedParticlesAgent;
    end
    
    % error output
    execTimePerStep(step) = toc;
    error_agent = calcDistance_(trueTrajectory(1:2,step),estimatedTrajectory(1:2,step));
    if mod(step, 100)==0
        fprintf('Epoch %d ',step);
        fprintf('NumofAncs S1: %d ',numEstimatedAnchors(1, step));
        fprintf('NumofAncs S2: %d ',numEstimatedAnchors(2, step));
        % fprintf('NumofAncs S3: %d \n',numEstimatedAnchors(3, step));
        fprintf('PositionError agent: %d ',error_agent);
        fprintf('Time: %4.4f \n',execTimePerStep(step));
        % fprintf('--------------------------------------------------- \n\n')
    end

    %% Plot for test (real time)
    if parameters.Isplotrt
    
        % plot estimated anchors and agent parameters.trueTraj
        h_traj = plot(estimatedTrajectory(1,step),estimatedTrajectory(2,step),'color',[0 .5 0],'marker','+', 'MarkerSize', 5,'LineWidth', 1.5);
        h_true = plot(parameters.trueTraj(1,step),parameters.trueTraj(2,step),'color',[0.6 .5 0],'marker','+', 'MarkerSize', 3,'LineWidth', 1.5);
        
        drawnow;
    end

end

