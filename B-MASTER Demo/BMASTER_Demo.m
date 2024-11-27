%% Tuning parameters
RandSeed = 1;      % randon generator seed
NumIter = 100;     % number of mcmc iterations, 100 works fine, 1000 recommeneded.
BurnIn = 100;      % burnin
conf = 0.90;       % sparsity cut-off
trainProp = 0.8;   % training proportion

%% Other tuning parameters

DisplayExtraUpdates = 0; % Set 1 to print detailed updates each iteration

%% Data Reading, dividing into Train and Test
rng(RandSeed)
X_raw = csvread("X.csv");
Y_raw = csvread("Y.csv");
N_all = size(X_raw,1);
indexes = 1:N_all;
shuffledIndexes = indexes(randperm(N_all));
N_Train = round(trainProp * N_all);
N_Test = N_all - N_Train;
indexesTrain = shuffledIndexes(1:N_Train);
indexesTest = shuffledIndexes(N_Train+1:end);
X_rawTrain = X_raw(indexesTrain,:);
X_rawTest = X_raw(indexesTest,:);
Y_rawTrain = Y_raw(indexesTrain,:);
Y_rawTest = Y_raw(indexesTest,:);

tic;

% Transforming X and Y to have mean = 0, sd = 1 for each coordinate
[X,Y, X_raw_mean, X_raw_std, Y_raw_mean, Y_raw_std] = XYraw2XY_withMeanSd(X_rawTrain,Y_rawTrain);
N = size(X_raw,1);
P = size(X_raw,2);
Q = size(Y_raw,2);
C = ones(P,Q);

%% MCMC hyper-parameters
sum_C = sum(C,[1,2]);
r_1 = 1;             % shape lambda1square (prior), mean = shape/rate
delta_1 = .1;        % rate lambda1square (prior)
r_2 = 1;             % shape lambda2 (prior), mean = shape/rate
delta_2 = .1;        % rate lambda2 (prior)
shape_a1 = sum_C + r_1;
shape_a2 = 0.5*sum_C + r_2;

%% Hyperparameter values and variables' initial values

B_initial = 0.01*ones(P,Q);
tauSquareMat_initial = 1*ones(P,Q);
gammapSquareArray_initial = 1*ones(P,1);
sigmaqSquareArray_initial = 1*ones(Q,1);
lambda1square_initial = 10^(1);
lambda2square_initial = 10^(1);


%% Running MCMC
[B_values,lambda1square_values,lambda2values] = MCMC_BMASTER_v4(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
    B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
    lambda1square_initial, lambda2square_initial);

%% Summary
% --> IndicatorNonzero_est = Indicator of non-zero coefficients
% --> B_post_mean = Posterior mean estimate AFTER dropping non-significant
% entries
% --> B_est = Posterior mean estimate AFTER dropping non-significant
% entries
% --> Note that, B_est is estimated based on transformed data

B_post_mean = mean(B_values(:,:,BurnIn:NumIter),3);
lambda_1_post_mean = sqrt(mean(lambda1square_values(BurnIn:NumIter)));
lambda_2_mean = mean(lambda2values(BurnIn:NumIter));
[pvalMat,IndicatorNonzero_est] = DeriveIndicatorNonzero(B_values,BurnIn,NumIter,conf);
B_est = B_post_mean.*IndicatorNonzero_est;  % Estimate on transformed data

%% Prediction on Test data
% --> Earlier, B_est is estimated based on transformed data. So to find
% prediction error on test data, we transform X test data, make prediction
% of Y in transferred domain, and then again transform that to original Y.

X_raw_mean_mat = repmat(X_raw_mean,N_Test,1);
X_raw_std_mat = repmat(X_raw_std,N_Test,1);
X_Test_scaled = (X_rawTest-X_raw_mean_mat)./X_raw_std_mat;
Y_Test_predicted = X_Test_scaled * B_est;

Y_raw_mean_mat = repmat(Y_raw_mean,N_Test,1);
Y_raw_std_mat = repmat(Y_raw_std,N_Test,1);
Y_rawTest_predicted = Y_Test_predicted.*Y_raw_std_mat + Y_raw_mean_mat;
RMSE = sqrt(mean((Y_rawTest(:) - Y_rawTest_predicted(:)).^2));
MAD = mean(abs(Y_rawTest(:) - Y_rawTest_predicted(:)));
RMSE_MAD = round([RMSE, MAD], 3);
