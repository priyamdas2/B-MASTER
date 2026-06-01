clearvars; clc; close all;
addpath('./Dataset/');

delta_1 = 0.1;        % rate lambda1square (prior)
delta_2 = 0.1;        % rate lambda2 (prior), delta_2 decrese <=> lambda2 increase, .1 value

%% Tuning parameters

RunNowMCMC = 1;
DisplayExtraUpdates = 0;
overWriteOutput = 1;
NumIter = 2000;        % number of mcmc iterations
BurnIn = 100;          % burnin
conf = 0.90;
AUCgridsize = 0.0001;

%% Data Reading

X_raw = csvread("X.csv");
IndicatorNonzeroTrue = csvread("NonZeroLocation_BMASTER_estimated_NumIter_2000_burnin_100.csv");
AllCoeffsTrue = csvread("B_BMASTER_estimated_NumIter_2000_burnin_100.csv");
TopMany = 5;
[rowIdx, colIdx, topVals] = topKAbsEntries(AllCoeffsTrue, TopMany);

rep = 1;
RandSeed = rep;
tic;
filename = ['Y_rep_',num2str(rep),'.csv'];
Y_raw = csvread(filename);
% Transforming X and Y to have (mean = 0, sd = 1)
[X,Y] = XYraw2XY(X_raw,Y_raw);
N = size(X_raw,1);
P = size(X_raw,2);
Q = size(Y_raw,2);
C = ones(P,Q);

%% MCMC hyper-parameters
sum_C = sum(C,[1,2]);
r_1 = 1;           % shape lambda1square (prior), mean = shape/rate
r_2 = 1;           % shape lambda2 (prior), mean = shape/rate
shape_a1 = sum_C + r_1;
shape_a2 = 0.5*sum_C + r_2;
%% Hyperparameter values and variables' initial values
B_initial = 0.01*ones(P,Q);
tauSquareMat_initial = 1*ones(P,Q);
gammapSquareArray_initial = 1*ones(P,1);
sigmaqSquareArray_initial = 1*ones(Q,1);
lambda1square_initial = 10^(1);
lambda2square_initial = 10^(1);


%% Filename for cached raw MCMC output

filename = ['Output iter sensitivity/BMASTER_RAW_MCMC_ITER_Sensitivity_NumIter_', ...
    num2str(NumIter), ...
    '_burnin_', ...
    num2str(BurnIn), ...
    '.mat'];

%% Run B-MASTER only if cached file does not already exist

if ~isfile(filename)
    
    [B_values,lambda1square_values,lambda2square_values] = MCMC_BMASTER_v4(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
        B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
        lambda1square_initial, lambda2square_initial);
    
    % Save raw MCMC outputs
    
    save(filename, ...
        'B_values', ...
        'lambda1square_values', ...
        'lambda2square_values', ...
        'NumIter', ...
        'BurnIn', ...
        'RandSeed', ...
        '-v7.3');
    
    Time_taken = toc;
    
    filename = ['Output iter sensitivity/ITER_Sensitivity_compTime_NumIter_',num2str(NumIter),'.csv'];
    csvwrite(filename,Time_taken)
    
else
    
    % Load previously saved MCMC output
    
    load(filename)    
end

%% Convergence diagnostics
% Example call
stats = mcmc_post_analysis(B_values, rowIdx, colIdx, ...
                           lambda1square_values, lambda2square_values, BurnIn);

% Quick table of MCSE/SD% and Geweke z (for reporting)
names     = string({stats.name})';
mcse_pct  = [stats.mcse_pct]';
geweke_z  = [stats.geweke_z]';
T = table(names, mcse_pct, geweke_z)

% Example plots outside the function:
j = 1;  % e.g., lambda1^2

x = stats(j).trace;

figure;
plot(x, '-', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1);

xlabel('Iteration (post-burn-in)');
ylabel(stats(j).name);
title(sprintf('%s: trace plot', stats(j).name));

plot_mcmc_traces_ITER_sensitivity(stats, 'File', ...
    sprintf('Output iter sensitivity/ITER_mcmc_traces_NumIter_%d.png', NumIter));

%% Summary
    
    B_post_mean = mean(B_values(:,:,BurnIn:NumIter),3);
    lambda_1_post_mean = sqrt(mean(lambda1square_values(BurnIn:NumIter)));
    lambda_2_post_mean = sqrt(mean(lambda2square_values(BurnIn:NumIter)));
    [pvalMat,IndicatorNonzero_est] = DeriveIndicatorNonzero(B_values,BurnIn,NumIter,conf);
    B_est = B_post_mean.*IndicatorNonzero_est;
    TrueSparsity = 1 - sum(sum(IndicatorNonzeroTrue))/(P*Q);
    [TPR_FPR_MCC_AUC_AUC20_sparsity, Final_FPR_TPR] = Summary_BMASTER(B_values, AUCgridsize, BurnIn, NumIter, IndicatorNonzeroTrue, IndicatorNonzero_est, pvalMat, X, Y);
    
    % Y_est = X*B_est
    % norm(Y)
    % norm(Y_est)
    NumYinfluenced = sum(IndicatorNonzero_est,2);
    [NumYinfluencedOrdered, idx] = sort(NumYinfluenced, 'descend');
    MasterRanks = [idx,NumYinfluencedOrdered];
    
    
    
    filename = ['Output iter sensitivity/ITER_Sensitivity_NonZeroLocation_BMASTER_NumIter_',num2str(NumIter),'.csv'];
    csvwrite(filename,IndicatorNonzero_est)
    filename = ['Output iter sensitivity/ITER_Sensitivity_TPR_FPR_MCC_AUC_AUC20_sparsity_BMASTER_NumIter_',num2str(NumIter),'.csv'];
    csvwrite(filename,TPR_FPR_MCC_AUC_AUC20_sparsity)
    filename = ['Output iter sensitivity/ITER_Sensitivity_AUC_curve_BMASTER_NumIter_',num2str(NumIter),'.csv'];
    csvwrite(filename,Final_FPR_TPR)
    filename = ['Output iter sensitivity/ITER_Sensitivity_MasterRanks_NumIter_',num2str(NumIter),'.csv'];
    csvwrite(filename,MasterRanks)
    
