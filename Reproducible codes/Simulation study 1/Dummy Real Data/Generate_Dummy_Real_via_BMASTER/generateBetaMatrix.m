function B_True_nonsparse = generateBetaMatrix(RandSeed,P, Q)
% Generate a P x Q matrix with entries from either [-5, -1] or [1, 5]
rng(RandSeed)
B_True_nonsparse = rand(P, Q);
B = randi([0, 1], P, Q);
B_True_nonsparse(B == 0) = -5 + 4 * B_True_nonsparse(B == 0);  % Uniformly distributed in [-5, -1]
B_True_nonsparse(B == 1) = 1 + 4 * B_True_nonsparse(B == 1);   % Uniformly distributed in [1, 5]
end