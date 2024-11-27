RandSeed = 1;
IsRhoNonZero = 0;
SampleMultFactor = 1;
P = 20;
N = P*SampleMultFactor;

%% Data generation parameters
MasterNonZeroProp = 0.20; % proportion of variables which are Master predictors
partialMasterNonzeroProp = 0.25; % after masters are eliminated
rlb = 0.4;
rub = 0.6;
Q = P;

%% Tuning parameters
RunNowMCMC = 1;
DisplayExtraUpdates = 1;
overWriteOutput = 1;
NumIter = 100;        % number of mcmc iterations
BurnIn = 20;          % burnin
conf = 0.90;
AUCgridsize = 0.001;

%% Data generation

if(IsRhoNonZero == 1)
    rho = 0.5;
else
    rho = 0;
end
% X generation
X_raw = generate_X_Raw(RandSeed,rho, N, P);
RangeOfNonZeroPropInPartialMaster = [rlb, rub];

% Nonzero-indicator generator
IndicatorNonzeroTrue = generate_NonZeroIndicatorMat(RandSeed,...
    P, Q, MasterNonZeroProp, partialMasterNonzeroProp,...
    RangeOfNonZeroPropInPartialMaster);

% B_True generator
B_True_nonsparse = generateBetaMatrix(RandSeed,P, Q);
B_True = B_True_nonsparse.*IndicatorNonzeroTrue;

TrueSparsity = 1 - sum(sum(IndicatorNonzeroTrue))/(P*Q);
if(IsRhoNonZero == 1)
    filename = ['Output_Scalability_rhoNonZero_TrueSparsity_P_',num2str(P),'_NumIter_',num2str(NumIter),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
    csvwrite(filename,round(TrueSparsity ,3));
else
    filename = ['Output_Scalability_rhoZero_TrueSparsity_P_',num2str(P),'_NumIter_',num2str(NumIter),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
    csvwrite(filename,round(TrueSparsity ,3));
end


Y_raw = X_raw*B_True + normrnd(0,0.1,[N,Q]);




% Transforming X and Y to have (mean = 0, sd = 1)
[X,Y] = XYraw2XY(X_raw,Y_raw);
%% MCMC hyper-parameters
C = ones(P,Q);
sum_C = sum(C,[1,2]);
r_1 = 1;           % shape lambda1square (prior), mean = shape/rate
delta_1 = 0.1;        % rate lambda1square (prior)
r_2 = 1;           % shape lambda2 (prior), mean = shape/rate
delta_2 = 0.1;        % rate lambda2 (prior), delta_2 decrese <=> lambda2 increase, .1 value
shape_a1 = sum_C + r_1;
shape_a2 = 0.5*sum_C + r_2;
%% Hyperparameter values and variables' initial values
B_initial = 0.01*ones(P,Q);
tauSquareMat_initial = 1*ones(P,Q);
gammapSquareArray_initial = 1*ones(P,1);
sigmaqSquareArray_initial = 1*ones(Q,1);
lambda1square_initial = 10^(-10);
lambda2square_initial = 10^(2);

[B_values,lambda1square_values,lambda2square_values] = MCMC_BMASTER_v4(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
    B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
    lambda1square_initial, lambda2square_initial);
Time_taken = toc;
%% Summary

B_post_mean = mean(B_values(:,:,BurnIn:NumIter),3);
lambda_1_post_mean = sqrt(mean(lambda1square_values(BurnIn:NumIter)));
lambda_2_post_mean = sqrt(mean(lambda2square_values(BurnIn:NumIter)));
[pvalMat,IndicatorNonzero_est] = DeriveIndicatorNonzero(B_values,BurnIn,NumIter,conf);
B_est = B_post_mean.*IndicatorNonzero_est;
TrueSparsity = 1 - sum(sum(IndicatorNonzeroTrue))/(P*Q);
[TPR_FPR_MCC_AUC_AUC20_sparsity, Final_FPR_TPR] = Summary_BMASTER(B_values, AUCgridsize, BurnIn, NumIter, IndicatorNonzeroTrue, IndicatorNonzero_est, pvalMat, X, Y);

NumYinfluenced = sum(IndicatorNonzero_est,2);
[NumYinfluencedOrdered, idx] = sort(NumYinfluenced, 'descend');
MasterRanks = [idx,NumYinfluencedOrdered];
TPR_FPR_MCC_AUC_AUC20_sparsity


if(NumIter >= 100)
    if(IsRhoNonZero == 1)
        %filename = ['Output_Scalability_rhoNonZero_NonZeroLocation_P_',num2str(P),'_NumIter_',num2str(NumIter),'.csv'];
        %csvwrite(filename,IndicatorNonzero_est)
        filename = ['Output_Scalability_rhoNonZero_summary_P_',num2str(P),'_NumIter_',num2str(NumIter),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
        csvwrite(filename,TPR_FPR_MCC_AUC_AUC20_sparsity)
        %filename = ['Output_Scalability_rhoNonZero_AUC_curve_P_',num2str(P),'_NumIter_',num2str(NumIter),'.csv'];
        %csvwrite(filename,Final_FPR_TPR)
        %filename = ['Output_Scalability_rhoNonZero_MasterRanks_P_',num2str(P),'_NumIter_',num2str(NumIter),'.csv'];
        %csvwrite(filename,MasterRanks)
        filename = ['Output_Scalability_rhoNonZero_compTime_P_',num2str(P),'_NumIter_',num2str(NumIter),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
        csvwrite(filename,Time_taken)
    else
        %filename = ['Output_Scalability_rhoZero_NonZeroLocation_P_',num2str(P),'_NumIter_',num2str(NumIter),'.csv'];
        %csvwrite(filename,IndicatorNonzero_est)
        filename = ['Output_Scalability_rhoZero_summary_P_',num2str(P),'_NumIter_',num2str(NumIter),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
        csvwrite(filename,TPR_FPR_MCC_AUC_AUC20_sparsity)
        %filename = ['Output_Scalability_rhoZero_AUC_curve_P_',num2str(P),'_NumIter_',num2str(NumIter),'.csv'];
        %csvwrite(filename,Final_FPR_TPR)
        %filename = ['Output_Scalability_rhoZero_MasterRanks_P_',num2str(P),'_NumIter_',num2str(NumIter),'.csv'];
        %csvwrite(filename,MasterRanks)
        filename = ['Output_Scalability_rhoZero_compTime_P_',num2str(P),'_NumIter_',num2str(NumIter),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
        csvwrite(filename,Time_taken)
    end
end
