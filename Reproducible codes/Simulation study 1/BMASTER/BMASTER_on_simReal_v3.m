%% Tuning parameters

NumExp = 1;
RunNowMCMC = 1;
DisplayExtraUpdates = 0;
overWriteOutput = 1;
NumIter = 1000;        % number of mcmc iterations
BurnIn = 100;          % burnin
conf = 0.90;
AUCgridsize = 0.0001;

%% Data Reading

X_raw = csvread("X.csv");
IndicatorNonzeroTrue = csvread("NonZeroLocation_BMASTER_estimated_NumIter_2000_burnin_100.csv");

for rep = 1:NumExp
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
    
    Time_taken = toc;
    
    filename = ['Output_simReal_NonZeroLocation_BMASTER_rep_',num2str(rep),'.csv'];
    csvwrite(filename,IndicatorNonzero_est)
    filename = ['Output_simReal_TPR_FPR_MCC_AUC_AUC20_sparsity_BMASTER_rep_',num2str(rep),'.csv'];
    csvwrite(filename,TPR_FPR_MCC_AUC_AUC20_sparsity)
    filename = ['Output_simReal_AUC_curve_BMASTER_rep_',num2str(rep),'.csv'];
    csvwrite(filename,Final_FPR_TPR)
    filename = ['Output_MasterRanks_rep_',num2str(rep),'.csv'];
    csvwrite(filename,MasterRanks)
    filename = ['Output_compTime_rep_',num2str(rep),'.csv'];
    csvwrite(filename,Time_taken)
end