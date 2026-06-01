clear; clc;

addpath('./BMASTER supp files/');

%% Tuning parameters

RunNowMCMC = 1;
overWriteOutput = 1;
DisplayExtraUpdates = 0;
NumIter = 1000;        % number of mcmc iterations
BurnIn = 100;          % burnin

%% Data Reading

X_raw = csvread("Data/X.csv");

RandSeed = 1;
tic;
filename = ['Data/Y.csv'];
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
delta_1 = .1;        % rate lambda1square (prior)
r_2 = 1;           % shape lambda2 (prior), mean = shape/rate
delta_2 = .1;        % rate lambda2 (prior), delta_2 decrese <=> lambda2 increase, .1 value
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

filename = ['Output/BMASTER_RAW_MCMC_NumIter_', ...
    num2str(NumIter), ...
    '_burnin_', ...
    num2str(BurnIn), ...
    '.mat'];

%% Run B-MASTER only if cached file does not already exist

if ~isfile(filename)
    
    [B_values,lambda1square_values,lambda2values] = MCMC_BMASTER_v4(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
        B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
        lambda1square_initial, lambda2square_initial);
    
    % Save raw MCMC outputs
    
    save(filename, ...
        'B_values', ...
        'lambda1square_values', ...
        'lambda2values', ...
        'NumIter', ...
        'BurnIn', ...
        'RandSeed', ...
        '-v7.3');
    
    Time_taken = toc;
    filename = ['Output/Output_REAL_compTime.csv'];
    csvwrite(filename,Time_taken)
    
else
    % Load previously saved MCMC output
    load(filename)
end

%% Summary

B_post_mean = mean(B_values(:,:,BurnIn:NumIter),3);
lambda_1_post_mean = sqrt(mean(lambda1square_values(BurnIn:NumIter)));
lambda_2_mean = mean(lambda2values(BurnIn:NumIter));



conf_list = [0.90, 0.95, 0.99];

for cc = 1:length(conf_list)
    
    conf = conf_list(cc);
    
    [pvalMat, IndicatorNonzero_est] = DeriveIndicatorNonzero(B_values,BurnIn,NumIter,conf);
    B_est = B_post_mean .* IndicatorNonzero_est;
    
    NumYinfluenced = sum(IndicatorNonzero_est,2);
    [NumYinfluencedOrdered, idx] = sort(NumYinfluenced, 'descend');
    MasterRanks = [idx, NumYinfluencedOrdered];
    
    filename = ['Output/Output_conf_',num2str(round(conf*100)),'_REAL_B_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
    csvwrite(filename,B_est)
    
    filename = ['Output/Output_conf_',num2str(round(conf*100)),'_REAL_NonZeroLocation_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
    csvwrite(filename,IndicatorNonzero_est)
    
    filename = ['Output/Output_conf_',num2str(round(conf*100)),'_REAL_MasterRanks_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
    csvwrite(filename,MasterRanks)
    
end

filename = ['Output/Output_REAL_B_NonSparse_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
csvwrite(filename,B_post_mean)

filename = ['Output/Output_REAL_pValues_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
csvwrite(filename,pvalMat)

