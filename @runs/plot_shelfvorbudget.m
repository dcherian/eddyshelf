function [] = plot_shelfvorbudget(runs)

    if ~isfield(runs.vorbudget, 'time')
        error('Vorbudget terms haven''t been calculated');
    end

    time = runs.vorbudget.time/86400;

    figure;
    subplot(2,1,1)
    [ax,h1,h2] = plotyy(time, runs.vorbudget.shelf.rv, runs.csflux.time/86400, ...
                        runs.csflux.west.shelf/1e6);
    set(ax(2), 'XTick', []);
    ylabel(ax(2), 'Cross-shelfbreak Transport (Sv)')

    beautify;
    limx = xlim;

    xlabel('Time (days)');
    ylabel('Volume averaged relative vorticity (shelf water)');
    liney(0, [], [1 1 1]*0. 5);

    subplot(2,1,2)
    hold all
    plot(time, runs.vorbudget.shelf.str, 'Color', [0.68 0.85 ...
                        0.90]);
    plot(time, runs.vorbudget.shelf.hadv + runs.vorbudget.shelf.vadv, ...
         time, runs.vorbudget.shelf.tilt, ...
         time, runs.vorbudget.shelf.bfric, ...
         time, runs.vorbudget.shelf.beta);
    plot(time, runs.vorbudget.shelf.hadv, 'Color',[1 1 1]*0.6);
    plot(time, runs.vorbudget.shelf.vadv, 'Color',[1 1 1]*0.6);
    plot(time, smooth(runs.vorbudget.shelf.str, 3), 'Color', [0 0 1]);
    xlim(limx);
    liney(0, [], [1 1 1]*0.5);
    xlabel('Time (days)');
    ylabel('sec^{-2}');
    legend('str', 'adv', 'tilt', 'bfric', 'beta', 'hadv', ...
           'vadv', 'Location', 'NorthWest');
    beautify
end
