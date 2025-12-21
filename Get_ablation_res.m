% Analyze_Ablation_Detailed
clear; clc;

filename = 'Res_ablation.mat'; 

if ~exist(filename, 'file')
    error('文件 %s 不存在，请检查路径。', filename);
end

load(filename); 

fprintf('==========================================================================================\n');
fprintf('Detailed Analysis for: %s\n', filename);
fprintf('==========================================================================================\n');
fprintf('%-20s | %-16s | %-16s | %-16s | %-16s | %-10s\n', ...
    'Method', 'RMSE(Avg)±Std', 'RMSE(Final)', 'OSPA(Avg)', 'OSPA(Final)', 'Time(s)');
fprintf('------------------------------------------------------------------------------------------\n');

num_scenarios = length(results);

for k = 1:num_scenarios
    rmse_matrix = results(k).rmse;
    val_rmse_avg_total = mean(rmse_matrix(:), 'omitnan');
    rmse_per_run = mean(rmse_matrix, 1, 'omitnan');
    val_rmse_std = std(rmse_per_run, 'omitnan');
    last_steps = 10; 
    val_rmse_final = mean(mean(rmse_matrix(end-last_steps+1:end, :), 1, 'omitnan'), 'omitnan');
    
    mean_curve = mean(rmse_matrix, 2, 'omitnan');
    val_rmse_max = max(mean_curve);

    ospa_matrix = results(k).ospa;
    
    ospa_combined = squeeze(mean(ospa_matrix, 1, 'omitnan'));
    val_ospa_avg = mean(ospa_combined(:), 'omitnan');
    val_ospa_final = mean(mean(ospa_combined(end-last_steps+1:end, :), 1, 'omitnan'), 'omitnan');
    val_time = mean(results(k).exec_time);

    str_rmse_full = sprintf('%.3f ± %.3f', val_rmse_avg_total, val_rmse_std);
    
    fprintf('%-20s | %-16s | %-16.4f | %-16.4f | %-16.4f | %-10.2f\n', ...
        results(k).name, str_rmse_full, val_rmse_final, val_ospa_avg, val_ospa_final, val_time);
end
fprintf('------------------------------------------------------------------------------------------\n');