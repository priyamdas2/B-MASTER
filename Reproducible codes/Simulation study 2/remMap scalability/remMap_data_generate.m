RandSeed = 1;
IsRhoNonZero = 1;
SampleMultFactor = 1;
P = 60;
N = P*SampleMultFactor;

%% Data generation parameters
MasterNonZeroProp = 0.20; % proportion of variables which are Master predictors
partialMasterNonzeroProp = 0.25; % after masters are eliminated
rlb = 0.4;
rub = 0.6;
Q = P;
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

Y_raw = X_raw*B_True + normrnd(0,0.1,[N,Q]);
filename = ['Data_X_P_',num2str(P),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
    csvwrite(filename,X_raw);
filename = ['Data_Y_P_',num2str(P),'_Nfactor_',num2str(SampleMultFactor),'.csv'];
    csvwrite(filename,Y_raw);
  



