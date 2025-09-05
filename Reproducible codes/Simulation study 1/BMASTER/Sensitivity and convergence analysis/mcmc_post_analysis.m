function stats = mcmc_post_analysis(B_values, rowIdx, colIdx, ...
    lambda1square_values, lambda2square_values, BurnIn)

    NumIter = size(B_values, 3);
    if any([numel(lambda1square_values), numel(lambda2square_values)] ~= NumIter)
        error('lambda vectors must have length equal to size(B_values,3).');
    end
    if numel(rowIdx) ~= 5 || numel(colIdx) ~= 5
        error('rowIdx and colIdx must have length 5.');
    end

    series = cell(1,7); labels = cell(1,7);
    series{1} = lambda1square_values(:); labels{1} = 'lambda1^2';
    series{2} = lambda2square_values(:); labels{2} = 'lambda2^2';
    for k = 1:5
        series{2+k} = squeeze(B_values(rowIdx(k), colIdx(k), :));
        labels{2+k} = sprintf('B(%d,%d)', rowIdx(k), colIdx(k));
    end

    stats = repmat(struct('name',[], 'trace',[], 'cmean',[], ...
        'mean',[], 'sd',[], 'mcse',[], 'mcse_pct',[], 'geweke_z',[]),1,7);

    for j = 1:7
        x = series{j};
        if BurnIn >= numel(x)-2, error('BurnIn too large for series %d.', j); end
        x = x(BurnIn+1:end); x = x(:);
        n = numel(x);

        cmean      = cumsum(x)./(1:n)';
        mu         = mean(x);
        sd         = std(x,0);                 % sample SD
        mcse_val   = mcse_batchmeans(x);       % <-- FIXED
        mcse_pct   = 100 * (mcse_val / sd);
        gz         = geweke_zscore(x, 0.2, 0.5);

        stats(j).name     = labels{j};
        stats(j).trace    = x;
        stats(j).cmean    = cmean;
        stats(j).mean     = mu;
        stats(j).sd       = sd;
        stats(j).mcse     = mcse_val;
        stats(j).mcse_pct = mcse_pct;
        stats(j).geweke_z = gz;
    end
end

% ===== Helpers =====

function m = mcse_batchmeans(x)
% Correct batch-means MCSE:
%   - n = a*b with batch length b and a batches
%   - batch means: x̄_b(1),...,x̄_b(a)
%   - var(mean(x)) ≈ Var_batch_means / a
%   - MCSE = sqrt(Var_batch_means / a)

    x = x(:); n = numel(x);
    b = max(20, floor(sqrt(n)));     % rule-of-thumb batch length
    a = floor(n / b);                % number of batches
    if a < 2
        m = std(x,0)/sqrt(n);        % fallback (essentially IID)
        return;
    end
    n_use = a*b;
    x = x(1:n_use);

    xb = mean(reshape(x, b, a), 1)'; % a×1 vector of batch means
    % use population variance (denominator a) for Var_batch_means
    var_batch_means = var(xb, 1);
    var_mean        = var_batch_means / a;     % <-- key correction
    m               = sqrt(max(var_mean, 0));
end

function z = geweke_zscore(x, first_frac, last_frac)
    x = x(:); n = numel(x);
    n1 = max(20, floor(first_frac * n));
    n2 = max(20, floor(last_frac  * n));
    idx2_start = n - n2 + 1;
    if idx2_start <= n1, n1 = floor(0.1*n); n2 = floor(0.5*n); idx2_start = n - n2 + 1; end

    x1 = x(1:n1); x2 = x(idx2_start:end);
    m1 = mean(x1); m2 = mean(x2);
    L1 = max(1, floor(numel(x1)^(1/3)));
    L2 = max(1, floor(numel(x2)^(1/3)));
    s0_1 = nw_spectral0(x1, L1);
    s0_2 = nw_spectral0(x2, L2);
    v = s0_1/numel(x1) + s0_2/numel(x2);
    z = (m1 - m2) / sqrt(max(v, eps));
end

function s0 = nw_spectral0(x, L)
    x = x(:); n = numel(x);
    x = x - mean(x);
    g0 = (x.'*x)/n;
    if L == 0, s0 = g0; return; end
    g = zeros(L,1);
    for ell = 1:L
        g(ell) = (x(1:n-ell)' * x(1+ell:n)) / n;
    end
    w = 1 - (1:L)'/(L+1);             % Bartlett weights
    s0 = max(g0 + 2 * sum(w .* g), 0);
end
