function [TPR_FPR_MCC_AUC_AUC20_sparsity, Final_FPR_TPR] = Summary_BMASTER(B_values, AUCgridsize, BurnIn, NumIter, IndicatorNonzeroTrue, IndicatorNonzero_est, pvalMat, X, Y) 
P = size(B_values,1);
Q = size(B_values,2);

B_post_mean = mean(B_values(:,:,BurnIn:NumIter),3);

%% TP, FP, TN, FN, TPR, FPR, MCC 

[aa, bb] = find(IndicatorNonzero_est==1);
sparsity = 1 - size(aa,1)/(P*Q);


[TrueNzRownums, TrueNzColnums] = find(IndicatorNonzeroTrue==1);

NumPos = size(TrueNzRownums,1);
NumTruePos = 0;
for i = 1:NumPos
    NumTruePos = NumTruePos + IndicatorNonzero_est(TrueNzRownums(i), TrueNzColnums(i));
end
NumFalseNeg = NumPos - NumTruePos;


[TrueZeroRownums, TrueZeroColnums] = find(IndicatorNonzeroTrue==0);
NumNeg = size(TrueZeroRownums,1);
NumFalsePos = 0;
for i = 1:NumNeg
    NumFalsePos = NumFalsePos + IndicatorNonzero_est(TrueZeroRownums(i), TrueZeroColnums(i));
end
NumTrueNeg = NumNeg - NumFalsePos;

TPR = NumTruePos/NumPos;
FPR = NumFalsePos/NumNeg;
MCC = (NumTruePos*NumTrueNeg - NumFalsePos*NumFalseNeg)/sqrt((NumTruePos + NumFalsePos)*(NumTruePos + NumFalseNeg)*...
    (NumTrueNeg + NumFalsePos)*(NumTrueNeg + NumFalseNeg));


%% AUC


pvalCuts = 0:AUCgridsize:1;

NumCuts = length(pvalCuts);
FPR_TPR = zeros(NumCuts,2);

for h = 1:NumCuts
    fprintf('\n');
    fprintf('=> AUC calculation complete: %.1f percentage.', 100*h/NumCuts);
    fprintf('\n');
    pval = h*AUCgridsize;

    IndicatorNonzero_temp = zeros(P,Q);

    for i = 1:P
        for j = 1:Q
            if(pvalMat(i,j) <= pval)
                IndicatorNonzero_temp(i,j) = 1;
            end
        end
    end

    NumTruePos_temp = 0;
    for i = 1:NumPos
        NumTruePos_temp = NumTruePos_temp + IndicatorNonzero_temp(TrueNzRownums(i), TrueNzColnums(i));
    end


    NumFalsePos_temp = 0;
    for i = 1:NumNeg
        NumFalsePos_temp = NumFalsePos_temp + IndicatorNonzero_temp(TrueZeroRownums(i), TrueZeroColnums(i));
    end

    FPR_TPR(h,:) = [NumFalsePos_temp/NumNeg, NumTruePos_temp/NumPos];
end

FPR_TPR_uniques = unique([FPR_TPR;[1,1]], 'rows');

Final_FPR_TPR = [[0,0];FPR_TPR_uniques];

AUC = trapz(Final_FPR_TPR(:,1), Final_FPR_TPR(:,2));

up20Idx = length(find(Final_FPR_TPR(:,1) <= 0.2));
Final_FPR_TPR_upto20 = Final_FPR_TPR(1:up20Idx,:);
if(Final_FPR_TPR(up20Idx,1) < 0.2)
    y1 = Final_FPR_TPR(up20Idx,2);
    y2 = Final_FPR_TPR(up20Idx+1,2);
    x1 = Final_FPR_TPR(up20Idx,1);
    x2 = Final_FPR_TPR(up20Idx+1,1);
    valAt20 = y1 + (0.2-x1)*(y2-y1)/(x2-x1);
    Final_FPR_TPR_upto20 = [Final_FPR_TPR_upto20; .2, valAt20];
end
   
AUC20 = 5*trapz(Final_FPR_TPR_upto20(:,1), Final_FPR_TPR_upto20(:,2));



%% MSE 

B_est = B_post_mean.*IndicatorNonzero_est;
Y_pred = X*B_est;
MSE = sum(sum(abs(Y - Y_pred).^2))/(P*Q);
TPR_FPR_MCC_AUC_AUC20_sparsity = round([TPR,FPR,MCC,AUC,AUC20,sparsity],4);

end
