clearvars; clc; close all;
addpath('./Dataset/');

% Sensitivity settings explored
delta_grid = [0.01 0.01;
              0.01 1;
              1    0.01;
              1    1];

NumExp = 10; % number of replicates per setting

% Preallocate result storage
Results_mean = zeros(size(delta_grid,1),6);
Results_se   = zeros(size(delta_grid,1),6);

for d = 1:size(delta_grid,1)
    delta1 = delta_grid(d,1);
    delta2 = delta_grid(d,2);

    all_reps = zeros(NumExp,6); % store [TPR FPR MCC AUC AUC20 Sparsity]

    for rep = 1:NumExp
        filename = sprintf('Output/Sensitivity_delta1_%g_delta2_%g_TPR_FPR_MCC_AUC_AUC20_sparsity_BMASTER_rep_%d.csv',...
                            delta1, delta2, rep);
        vals = csvread(filename);
        all_reps(rep,:) = vals(:)'; % ensure row vector
    end

    % Compute mean and standard error across replicates
    Results_mean(d,:) = mean(all_reps,1);
    Results_se(d,:)   = std(all_reps,0,1)/sqrt(NumExp);
end

% Display results in a table
VarNames = {'TPR','FPR','MCC','AUC','AUC20','Sparsity'};
RowNames = arrayfun(@(i) sprintf('delta1=%.2g, delta2=%.2g', ...
                 delta_grid(i,1), delta_grid(i,2)), 1:size(delta_grid,1), 'UniformOutput',false);

SummaryTable_mean = array2table(Results_mean, 'VariableNames', VarNames, 'RowNames', RowNames);
SummaryTable_se   = array2table(Results_se,   'VariableNames', VarNames, 'RowNames', RowNames);

disp('Mean performance across replicates:')
disp(SummaryTable_mean)

disp('Standard error across replicates:')
disp(SummaryTable_se)

% Optionally save to CSV
writetable(SummaryTable_mean, 'Output/SensitivitySummary_mean.csv','WriteRowNames',true);
writetable(SummaryTable_se,   'Output/SensitivitySummary_se.csv','WriteRowNames',true);
