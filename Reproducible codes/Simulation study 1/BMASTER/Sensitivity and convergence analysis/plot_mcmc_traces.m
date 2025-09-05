function fig = plot_mcmc_traces(stats, varargin)
% plot_mcmc_traces(stats, 'File', 'traces.pdf')
% - stats: 1x7 struct from mcmc_post_analysis
% Optional name-value:
%   'File'  : filename to save (pdf/png/eps). If omitted, just shows figure.

    p = inputParser;
    addParameter(p, 'File', '', @(s)ischar(s)||isstring(s));
    parse(p, varargin{:});
    outFile = string(p.Results.File);

    % --- Figure & layout ---
    fig = figure('Color','w','Units','pixels','Position',[100 100 1200 1600]);
    tiledlayout(fig, 4, 2, 'TileSpacing','compact','Padding','compact');

    nplots = numel(stats); % should be 7

    for k = 1:nplots
        nexttile;
        x  = stats(k).trace(:);
        rm = stats(k).cmean(:);

        % ---- names (Unicode) ----
        if k == 1
            nameStr = 'λ₁²';
        elseif k == 2
            nameStr = 'λ₂²';
        else
            nameStr = sprintf('B[%d]', k-2);
        end

        % --- Special handling for λ₁² ---
        if k == 1
            % Check magnitude
            medianVal = median(abs(x));
            if medianVal < 1e-10
                scaleFactor = 1e20;
                x  = x  * scaleFactor;
                rm = rm * scaleFactor;
                ylab = 'Value ×10^{-20}';
            else
                ylab = 'Value';
            end
        else
            ylab = 'Value';
        end

        % Plot
        plot(x, '-', 'LineWidth', 0.8); hold on;
        plot(rm, 'LineWidth', 1.4);
        grid on; box on;

        % Title
        ttl = sprintf('%s  |  MCSE/SD=%.1f%%,  Geweke z=%.2f', ...
                      nameStr, stats(k).mcse_pct, stats(k).geweke_z);
        title(ttl, 'Interpreter','none');

        xlabel('Iteration (post burn-in)');
        ylabel(ylab);

        if k == 1
            legend({'Trace','Running mean'}, 'Location','best', 'Box','off');
        end
    end

    if nplots < 8
        nexttile; axis off;
    end

    % --- Save if requested ---
    if strlength(outFile) > 0
        [~,~,ext] = fileparts(outFile);
        switch lower(ext)
            case {'.pdf','.png','.tif','.tiff','.jpg','.jpeg'}
                exportgraphics(fig, outFile, 'Resolution', 300);
            case '.eps'
                print(fig, outFile, '-depsc', '-painters');
            otherwise
                exportgraphics(fig, outFile + ".png", 'Resolution', 300);
        end
    end
end
