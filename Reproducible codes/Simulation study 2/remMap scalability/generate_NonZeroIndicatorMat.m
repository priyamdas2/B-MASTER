function IndicatorNonzeroTrue = generate_NonZeroIndicatorMat(RandSeed,...
    P, Q, MasterNonZeroProp, partialMasterNonzeroProp,...
    RangeOfNonZeroPropInPartialMaster)
rng(RandSeed)
rlb = RangeOfNonZeroPropInPartialMaster(1);
rub = RangeOfNonZeroPropInPartialMaster(2);

NumMasters = floor(P*MasterNonZeroProp); 
MasterIndices = sort(randperm(P, NumMasters));

RestIndices = setdiff(1:P,MasterIndices);

NumPartialMasters = floor(length(RestIndices)*partialMasterNonzeroProp);
tempIndices = sort(randperm(length(RestIndices), NumMasters));
partialMasterIndices = RestIndices(tempIndices);


IndicatorNonzeroTrue = zeros(P,Q);
IndicatorNonzeroTrue(MasterIndices,:) = ones(NumMasters, Q);

for ii = 1:NumPartialMasters
    idxHere = partialMasterIndices(ii);
    rand_num = rlb + (rub - rlb) * rand();
    NumNZ = floor(rand_num*Q);
    tempIndicesHere = sort(randperm(Q, NumNZ));
    IndicatorNonzeroTrue(idxHere, tempIndicesHere) = ones(1,length(tempIndicesHere));
end

end

