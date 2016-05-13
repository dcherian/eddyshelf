%        [handles] = plot_fluxts(runs, factor, isobath, source)
function [handles] = plot_fluxts(runs, factor, isobath, source)

    fluxin = double(runs.recalculateFlux(factor*runs.bathy.hsb, isobath, source));
    [start, stop] = runs.flux_tindices(fluxin);
    [avgflux, err] = runs.calc_avgflux(fluxin);
    [maxflux,maxloc] = runs.calc_maxflux(fluxin);

    avgflux = avgflux;
    maxflux = maxflux;
    err = err;
    tvec = runs.csflux.time / 86400;
    deltat = (tvec(stop)-tvec(start)) * 86400;
    ifluxvec = cumtrapz(tvec*86400, fluxin);
    iflux0 = ifluxvec(start);
    patchAlpha = 0.3;

    figure;
    insertAnnotation([runs.name '.plot_fluxts(' num2str(factor) ...
                      ', ' num2str(isobath) ', ' num2str(source) ')']);
    [handles.hax, handles.ts(1), handles.ts(2)] = ...
        plotyy(tvec, fluxin, tvec, ifluxvec./ifluxvec(end)*100);
    axes(handles.hax(1));
    ylabel('Flux (m^3/s)');
    handles.htitle = title(runs.name);
    xlabel('Time (day)');
    hold on;
    %handles.maxflx = plot(tvec(maxloc), fluxin(maxloc)/1000, ...
    %                      'x', 'Color', handles.ts(1).Color);
    handles.patch(1) = patch([tvec(start) tvec(stop) tvec(stop) tvec(start)], ...
                             [avgflux-err avgflux-err avgflux+err avgflux+err], ...
                             handles.ts(1).Color);
    handles.patch(1).FaceAlpha = patchAlpha;
    handles.patch(1).EdgeColor = 'none';
    handles.line(1) = plot([tvec(start) tvec(stop)], [1 1]*avgflux, ...
                           '--', 'Color', handles.ts(1).Color);
    handles.hax(1).YTick = sort(unique(round( ...
        [handles.hax(1).YTick avgflux-err avgflux avgflux+err], 2)));
    handles.hax(1).YAxis.Exponent = 3;
    handles.hax(1).YAxis.TickLabelFormat = '%.1f';

    axes(handles.hax(2));
    hold on
    handles.patch(2) = patch([tvec(start) tvec(stop) tvec(stop) tvec(start)], ...
                    [iflux0  iflux0+(avgflux - err)*deltat ...
                     iflux0+(avgflux + err)*deltat iflux0]/ifluxvec(end)*100, ...
                             handles.ts(2).Color);
    handles.patch(2).FaceAlpha = patchAlpha;
    handles.patch(2).EdgeColor = 'none';
    handles.line(2) = plot([tvec(start) tvec(stop)], ...
                  [iflux0 iflux0+avgflux*deltat]./ifluxvec(end)*100, ...
                  '--', 'Color', handles.ts(2).Color);
    ylabel('% Volume transported');
    handles.hax(2).YTick = [0 5 50 90 100];
    handles.hax(2).YAxis.TickLabelFormat = '%g\%';

    handles.leghandles = [handles.ts(1)  ...
                        handles.line(1)  ...
                        handles.patch(1)];

    [handles.hleg, handles.icons, handles.legplots, ~] = ...
        legend(handles.leghandles,  ...
               {'Instantaneous value' ; ...
                'Average'; ...
                '95% confidence interval'}, ...
               'Location', 'NorthWest');

    for ii = [1 2 3 5 7]
        handles.icons(ii).Color = 'k';
    end
    handles.icons(8).FaceColor = handles.icons(4).Color;
    handles.icons(8).FaceAlpha = patchAlpha;

    axes(handles.hax(1)); beautify([24 25 28]);
    handles.hax(1).YColor = handles.ts(1).Color;
    handles.hax(1).YLabel.Color = handles.ts(1).Color;
    linex(tvec([start stop]), [], [1 1 1]*0.85);
    LowerLines;

    axes(handles.hax(2)); beautify([24 25 28]);
    handles.hax(2).YLim(1) = 0;
    handles.hax(2).YColor = handles.ts(2).Color;
    handles.hax(2).YLabel.Color = handles.ts(2).Color;

    handles.hisolabel = runs.add_isobathlabel(isobath);
    handles.hisolabel.Position(1:2) = [0.05 0.05];
    handles.hisolabel.FontSize = 22;

    handles.hax(1).XTickLabel{find(handles.hax(1).XTick == tvec(start))} = 't_{start}';
    handles.hax(1).XTickLabel{find(handles.hax(1).XTick == tvec(stop))} = 't_{stop}';
end