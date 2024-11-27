function [B_True, IndicatorNonzeroTrue] = Generate_B_True(RandSeed,B_True_nonsparse,MasterIndices, HowManyNonMasterEachCol)
rng(RandSeed)
P = size(B_True_nonsparse,1);
Q = size(B_True_nonsparse,2);

NonMasterIndices = setdiff(1:P, MasterIndices);
IndicatorNonzeroTrue = zeros(P,Q);
IndicatorNonzeroTrue(MasterIndices,:) = 1;
for q = 1:Q
    indices = randperm(length(NonMasterIndices), HowManyNonMasterEachCol);
    values = NonMasterIndices(indices);
    IndicatorNonzeroTrue(values,q) = 1;
end
B_True = B_True_nonsparse.*IndicatorNonzeroTrue;
end

