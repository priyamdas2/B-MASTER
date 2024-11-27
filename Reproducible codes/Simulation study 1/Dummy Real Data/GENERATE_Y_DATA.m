NumExp = 10;
ErrSd = .1;
filename = ['B_BMASTER_estimated_NumIter_2000_burnin_100.csv'];
B_True = csvread(filename);
BetaScaleFactor = 1;
filename = ['X.csv'];
X_raw = csvread(filename);
filename = ['NonZeroLocation_BMASTER_estimated_NumIter_2000_burnin_100.csv'];
IndicatorNonzeroTrue = csvread(filename);


N = size(X_raw,1);
P = size(X_raw,2);
Q = size(B_True,2);

for rep = 1:NumExp
    rng(rep)
    Y_raw = X_raw*B_True*BetaScaleFactor + normrnd(0,ErrSd,[N,Q]);
    filename = ['Y_rep_',num2str(rep),'.csv'];
    csvwrite(filename,Y_raw)
end

TrueSparsity = 1 - sum(sum(IndicatorNonzeroTrue))/(P*Q);
filename = ['B_sparsity_True.csv'];
csvwrite(filename,TrueSparsity)