function X_raw = generate_X_Raw(RandSeed,rho, N, P)
rng(RandSeed)
X_raw = nan(N,P);
mu = zeros(P,1);
Sigma = nan(P,P);

for ii = 1:P
    for jj = 1:P
        Sigma(ii,jj) = rho^abs(ii-jj);
    end
end

for n = 1:N
    X_raw(n,:) = mvnrnd(mu, Sigma, 1);
end
end