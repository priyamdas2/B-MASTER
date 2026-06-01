function [pvalMat,IndicatorNonzero_est] = DeriveIndicatorNonzero(B_values,BurnIn,NumIter,conf)
P = size(B_values,1);
Q = size(B_values,2);
IndicatorNonzero_est = zeros(P,Q);
pvalMat = zeros(P,Q);

for i = 1:P
    for j = 1:Q
        quan_low = quantile(B_values(i,j,BurnIn:NumIter), (1-conf)/2);
        quan_up = quantile(B_values(i,j,BurnIn:NumIter), conf +(1-conf)/2);
        pvalMat(i,j) = 2*min(sum(B_values(i,j,BurnIn:NumIter) > 0), ...
            sum(B_values(i,j,BurnIn:NumIter) < 0))/(NumIter - BurnIn + 1);
        if(quan_low*quan_up > 0)
            IndicatorNonzero_est(i,j) = 1;
        end
    end
end
end