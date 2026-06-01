function fig = plot_mcmc_traces_ITER_sensitivity(stats, varargin)
% plot_mcmc_traces_ITER_sensitivity(stats, 'File', 'traces.pdf')
%
% - No running mean plotted
% - Purple trace plots
%
% Inputs:
%   stats : 1x7 struct from mcmc_post_analysis
%
% Optional name-value:
%   'File' : filename to save (pdf/png/eps)

    p = inputParser;
    addParameter(p, 'File', '', @(s)ischar(s)||isstring(s));
    parse(p, varargin{:});
    outFile = string(p.Results.File);

    %% Figure & layout

    fig = figure('Color','w','Units','pixels', ...
                 'Position',[100 100 1200 1600]);

    tiledlayout(fig, 4, 2, ...
                'TileSpacing','compact', ...
                'Padding','compact');

    nplots = numel(stats);

    %% Purple color
    purple_col = [0.4940 0.1840 0.5560];

    for k = 1:nplots

        nexttile;

        x = stats(k).trace(:);

        %% Names

        if k == 1
            nameStr = 'λ₁²';
        elseif k == 2
            nameStr = 'λ₂²';
        else
            nameStr = sprintf('B[%d]', k-2);
        end

        %% Special handling for λ₁²

        if k == 1

            medianVal = median(abs(x));

            if medianVal < 1e-10
                scaleFactor = 1e20;
                x = x * scaleFactor;
                ylab = 'Value ×10^{-20}';
            else
                ylab = 'Value';
            end

        else

            ylab = 'Value';

        end

        %% Plot trace only

        plot(x, ...
             '-', ...
             'Color', purple_col, ...
             'LineWidth', 0.9);

        grid on;
        box on;

        %% Title

        ttl = sprintf('%s  |  MCSE/SD=%.1f%%,  Geweke z=%.2f', ...
                      nameStr, ...
                      stats(k).mcse_pct, ...
                      stats(k).geweke_z);

        title(ttl, 'Interpreter','none');

        xlabel('Iteration (post burn-in)');
        ylabel(ylab);

    end

    %% Empty tile if needed

    if nplots < 8
        nexttile;
        axis off;
    end

    %% Save if requested

    if strlength(outFile) > 0

        [~,~,ext] = fileparts(outFile);

        switch lower(ext)

            case {'.pdf','.png','.tif','.tiff','.jpg','.jpeg'}

                exportgraphics(fig, outFile, 'Resolution', 300);

            case '.eps'

                print(fig, outFile, '-depsc', '-painters');

            otherwise

                exportgraphics(fig, outFile + ".png", ...
                               'Resolution', 300);

        end

    end

end