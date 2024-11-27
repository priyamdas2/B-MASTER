RandSeed = 1;
DisplayExtraUpdates = 0;
NumIter = 100;        % number of mcmc iterations
BurnIn = 20;          % burnin
conf = 0.90;

SampleMultFactor = 10;

P_array = [10, 20, 30, 40, 50, 60];
NumP = length(P_array);
Time_taken = nan(1,NumP);
for ii = 1:NumP
    P = P_array(ii);
    
    filename = ['Data_X_P_',num2str(P),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
    X_raw = csvread(filename);
    filename = ['Data_Y_P_',num2str(P),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
    Y_raw = csvread(filename);
    N = size(X_raw,1);
    P = size(X_raw,2);
    Q = size(Y_raw,2);
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
    tic;
    [B_values,lambda1square_values,lambda2square_values] = MCMC_BMASTER_v4(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
        B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
        lambda1square_initial, lambda2square_initial);
    Time_taken(ii) = toc;
end


Time_summary = [P_array; Time_taken];
filename = ['Time_BMASTER_Nfactor_',num2str(SampleMultFactor),'.csv'];
csvwrite(filename,Time_summary)


