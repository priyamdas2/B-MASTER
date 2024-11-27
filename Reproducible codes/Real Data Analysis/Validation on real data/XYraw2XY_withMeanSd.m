function [X,Y, X_raw_mean, X_raw_std, Y_raw_mean, Y_raw_std] = XYraw2XY_withMeanSd(X_raw,Y_raw)
N = size(X_raw,1);

X_raw_mean = mean(X_raw,1);
X_raw_std = std(X_raw,1);
X_raw_mean_mat = repmat(X_raw_mean,N,1);
X_raw_std_mat = repmat(X_raw_std,N,1);
X = (X_raw-X_raw_mean_mat)./X_raw_std_mat;

Y_raw_mean = mean(Y_raw,1);
Y_raw_std = std(Y_raw,1);
Y_raw_mean_mat = repmat(Y_raw_mean,N,1);
Y_raw_std_mat = repmat(Y_raw_std,N,1);
Y = (Y_raw-Y_raw_mean_mat)./Y_raw_std_mat;
end