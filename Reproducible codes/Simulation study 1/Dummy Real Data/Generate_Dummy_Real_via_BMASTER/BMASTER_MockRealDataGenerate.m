clear all

%% Tuning parameters
overWriteOutput = 0;

NumIter = 2000;        % number of mcmc iterations
BurnIn = 100;          % burnin
RandSeed = 1;
RunNowMCMC = 1;
conf = 0.95;
DisplayExtraUpdates = 1;

%% Data Generation/Reading

X_raw = csvread("Data/X.csv");
Y_raw = csvread("Data/Y.csv");

%% Transforming X and Y to have (mean = 0, sd = 1)

N = size(X_raw,1);
P = size(X_raw,2);
Q = size(Y_raw,2);
C = ones(P,Q);
[X,Y] = XYraw2XY(X_raw,Y_raw);
%% MCMC

if(RunNowMCMC == 1)
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
    
    [B_values,lambda1square_values,lambda2square_values] = MCMC_BMASTER_v3(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
    B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
    lambda1square_initial, lambda2square_initial);
        
    %% Posterior summary
    
    
    lambda_1_post_mean = sqrt(mean(lambda1square_values(BurnIn:NumIter)));
    lambda_2_post_mean = sqrt(mean(lambda2square_values(BurnIn:NumIter)));
    if(overWriteOutput == 1)
        save('Data/B_REAL_posteriorSamples.mat', 'B_values', '-v7.3');
    end
else
    N = size(X,1);
    B_values_temp = load('Data/B_REAL_posteriorSamples.mat');
    B_values = cell2mat(struct2cell(B_values_temp));
    P = size(B_post_mean,1);
    Q = size(B_post_mean,2);
end

%% Summary


B_post_mean = mean(B_values(:,:,BurnIn:NumIter),3);
[pvalMat,IndicatorNonzero_est] = DeriveIndicatorNonzero(B_values,BurnIn,NumIter,conf);
B_est = B_post_mean.*IndicatorNonzero_est;
EstSparsity = 1 - sum(sum(IndicatorNonzero_est))/(P*Q);
EstSparsity
Y_est = X*B_est;
norm(Y_est)
norm(Y - Y_est)

NumYinfluenced = sum(IndicatorNonzero_est,2);
[NumYinfluencedOrdered, idx] = sort(NumYinfluenced, 'descend');
MasterRanks = [idx,NumYinfluencedOrdered];

%% Only uncomment the following lines to overwrite

filename = ['Data/B_BMASTER_estimated_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
csvwrite(filename,B_est)
filename = ['Data/NonZeroLocation_BMASTER_estimated_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
csvwrite(filename,IndicatorNonzero_est)
filename = ['Data/MasterRanks_BMASTER_estimated_NumIter_',num2str(NumIter),'_burnin_',num2str(BurnIn),'.csv'];
csvwrite(filename,MasterRanks)
