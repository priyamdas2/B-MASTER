function [rowIdx, colIdx, topVals] = topKAbsEntries(Mat, TopMany)
% topKAbsEntries finds the top-K largest absolute entries in a matrix.
%
% Usage:
%   [rowIdx, colIdx, topVals] = topKAbsEntries(Mat, TopMany)
%
% Inputs:
%   Mat      - input matrix (numeric)
%   TopMany  - number of top entries to return
%
% Outputs:
%   rowIdx   - row indices of top entries
%   colIdx   - column indices of top entries
%   topVals  - corresponding values (with sign preserved)

    % flatten matrix and sort by absolute value
    [sortedVals, sortedIdx] = sort(abs(Mat(:)), 'descend');
    
    % restrict to TopMany
    topLinIdx = sortedIdx(1:TopMany);
    
    % convert linear indices to row/col
    [rowIdx, colIdx] = ind2sub(size(Mat), topLinIdx);
    
    % get signed values from original matrix
    topVals = Mat(topLinIdx);
end
