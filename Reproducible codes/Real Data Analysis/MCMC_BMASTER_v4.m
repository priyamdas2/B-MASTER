function [B_values,lambda1square_values,lambda2values, ErrorCumsum] = MCMC_BMASTER_v4(NumIter,RandSeed,DisplayExtraUpdates,X,Y,C,shape_a1,delta_1,shape_a2,delta_2,...
    B_initial,tauSquareMat_initial,gammapSquareArray_initial,sigmaqSquareArray_initial,...
    lambda1square_initial, lambda2square_initial)
%% Initial parameters

N = size(X,1);
P = size(X,2);
Q = size(Y,2);
XprodX = X'*X;
non_zero_C = C == 1;
ErrorCumsum = nan(NumIter,1);
YHatNorm = nan(NumIter,1);
%% Gibbs sampling initialization

B = B_initial;
tauSquareMat = tauSquareMat_initial;
gammapSquareArray = gammapSquareArray_initial;
sigmaqSquareArray = sigmaqSquareArray_initial;
lambda1square = lambda1square_initial;
lambda2square = lambda2square_initial;

%% Other parameters
B_values = zeros(P,Q,NumIter);
lambda1square_values = zeros(NumIter,1);
lambda2values = zeros(NumIter,1);

%% Gibbs MCMC
rng(RandSeed)
tic;
for iter = 1:NumIter
    
    time_now = toc;
    
    if(iter > 1)
        remainingTime = (time_now/(iter - 1))*(NumIter - iter + 1);
        fprintf('\n');
        fprintf('=> Performing iteration: %d, remaining approx. time: %.1f seconds.',iter, remainingTime);
        fprintf('\n');
        if(DisplayExtraUpdates == 1)
            fprintf('=> norm(Y - YHat) last iter: %.3f, norm(Y - YHat)/norm(Y) norm last iter: %.3f %%.',norm(Y - YHat), 100*norm(Y - YHat)/norm(Y));
            fprintf('\n');
            fprintf('=> norm(YHat) last iter: %.3f.',norm(YHat));
            fprintf('\n');
            fprintf('=> log10(lambda1) last iter: %.4f, log10(lambda2) last iter %.4f.',log10(sqrt(lambda1square)), log10(sqrt(lambda2square)));
            fprintf('\n');
            fprintf('=> Norm B last iter: %.4f',norm(B));
            fprintf('\n');
            fprintf('=> Proportion of B whose abs > 0.00001 last iter: %.4f',length(find(abs(B)>0.00001))/(P*Q));
            fprintf('\n');
            fprintf('=> mean(abs(mvn_var)) last iter: %.4f',mean(mean(abs(mvn_var))));
            fprintf('\n');
            fprintf('=> mean(abs(mvn_mean)) last iter: %.4f',mean(abs(mvn_mean)));
            fprintf('\n');
            fprintf('=> norm(sigmaqSquareArray) last iter: %.4f',norm(sigmaqSquareArray));
            fprintf('\n');
            fprintf('=> mean(abs(inv_part)) last iter: %.4f',mean(mean(abs(inv_part))));
            fprintf('\n');
            fprintf('=> log10(tauSquareMat(P,Q)) last iter: %.3f',log10(tauSquareMat(P,Q)));
            fprintf('\n');
            fprintf('=> log10(gammapSquareArray(P)) last iter: %.3f',log10(gammapSquareArray(i)));
            fprintf('\n');
        end
    
    end
    %% Updating B
    
    DArrayMat = zeros(P,Q);  % each column is used to form diagonal matrix
    
    for i = 1:P
        for j = 1:Q
            if(C(i,j) ~= 0)
                temp_1 = (tauSquareMat(i,j))^(-2) + (gammapSquareArray(i))^(-2);
                DArrayMat(i,j) = max(10^(-6), min(temp_1,10^6));   % to avoid error
            end
        end
    end

    for j = 1:Q
        D_j = diag(DArrayMat(:,j));
        inv_part = inv(XprodX + D_j + 10^(-6)*eye(P));
        mvn_mean = inv_part*(X'*Y(:,j));
        mvn_var = sigmaqSquareArray(j)*inv_part;
        B(:,j) = mvnrnd(mvn_mean,mvn_var);
    end
    
    B_values(:,:,iter) = B;
    
    %% Updating tauSquareMat
    
    tauInvGauss_b = lambda1square;
    
    for i = 1:P
        for j = 1:Q
            if(abs(B(i,j)) > 10^(-2))
                tauInvGauss_a = sqrt(lambda1square*sigmaqSquareArray(j)/((B(i,j)+ 10^(-6))^2));
                pd = makedist('InverseGaussian','mu',tauInvGauss_a,'lambda',tauInvGauss_b + 10^(-6));  % mean = mu
                tauSquareMat(i,j) = 1/max(10^(-30),random(pd));
            else
                mg = 1;
                pd_alt = gamrnd((mg+1)/2, lambda1square/2);
                tauSquareMat(i,j) = 1/max(10^(-30),pd_alt);
            end
        end
    end
    
    %% Updating gammapSquareArray
    
    gammaInvGauss_b = lambda2square;
    
    for i = 1:P
        if(norm(B(i,:)) > Q*10^(-2))
            temp_2 = (B(i,:).^2)./sigmaqSquareArray';
            gammaInvGauss_a = sqrt( lambda2square/(sum(temp_2)+ 10^(-10)) );
            pd2 = makedist('InverseGaussian','mu',gammaInvGauss_a,'lambda',gammaInvGauss_b + 10^(-10)); % increase pd2
            gammapSquareArray(i) = 1/max(10^(-14), random(pd2)); % increase random(pd2)
        else
            mg = Q;
            inv_gamma_scale = lambda2square/2;
            pd2_alt = gamrnd((mg+1)/2, 1/inv_gamma_scale);
            gammapSquareArray(i) = 1/max(10^(-14),pd2_alt);
        end   
    end
    
    %% Updating sigmaqSquareArray
    % Inv-gamma(shape =a,scale = b) = 1/gamma(shape = a, scale = b)
    ErrorCumsum(iter) = 0;
    YHatNorm(iter) = 0;
    YHat = nan(N,Q);
    for j = 1:Q
        invGammaShape_a = (N+sum(C(:,j)))/2;                           % shape (invgamma)
        D_j = diag(DArrayMat(:,j));
        YHat(:,j) = X*B(:,j);
        temp_3 = (Y(:,j) - YHat(:,j));
        invGammaScale_b = (temp_3'*temp_3 + (B(:,j)'*D_j)*B(:,j))/2;   % scale (invgamma); invgamma(shape_invgamma,scale_invgamma) = 1/gamma(shape = shape_gamma, scale = 1/scale_invgamma)
        pd3 = gamrnd(invGammaShape_a, 1/invGammaScale_b, [1,1]);         % gamrnd(shape, scale), Mean = shape*scale
        sigmaqSquareArray(j) = 1/max(10^(-6),pd3);                                  % inv-gamma mean = scale/(shape-1)
        ErrorCumsum(iter) = ErrorCumsum(iter) + norm(temp_3);
    end
    
    %% Updating lambda1square and lambda_2
    
    rate_b1 = 0.5*sum(tauSquareMat(non_zero_C)) + delta_1;    % rate
    lambda1square = gamrnd(shape_a1, 1/rate_b1);              % gamrnd(shape, scale), Mean = shape*scale
    lambda1square_values(iter) = lambda1square;
    
    temp_4 = sum(gammapSquareArray);
    rate_b2 = 0.5*temp_4+ delta_2;                           % rate
    lambda2square = gamrnd(shape_a2, 1/rate_b2);                   % gamrnd(shape, scale), Mean = shape*scale = shape/rate
    %rate_b2
    lambda2values(iter) = lambda2square;
end
end

