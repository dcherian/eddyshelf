%        [handles] = plot_fluxts(runs, factor, isobath, source)
function [handles] = plot_fluxts(runs, factor, isobath, source)

    fluxin = double(runs.recalculateFlux(factor*runs.bathy.hsb, isobath, source));
    [start, stop] = runs.flux_tindices(fluxin);
    [avgflux, err] = runs.calc_avgflux(fluxin);
    [maxflux,maxloc] = runs.calc_maxflux(fluxin);

    avgflux = avgflux/1000;
    maxflux = maxflux/1000;
    err = err/1000;
    tvec = runs.csflux.time / 86400;
    deltat = (tvec(stop)-tvec(start)) * 86400;
    ifluxvec = cumtrapz(tvec*86400, fluxin);
    iflux0 = ifluxvec(start);
    patchAlpha = 0.3;

    figure;
    insertAnnotation([runs.name '.plot_fluxts(' num2str(factor) ...
                      ', ' num2str(isobath) ', ' num2str(source) ')']);
    [handles.hax, handles.ts(1), handles.ts(2)] = ...
        plotyy(tvec, fluxin/1e3, tvec, ifluxvec);
    axes(handles.hax(1));
    ylabel('Flux (x 10^3 m^3/s)');
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
    correct_ticks('y', '%.2f', []);

    axes(handles.hax(2));
    hold on
    handles.patch(2) = patch([tvec(start) tvec(stop) tvec(stop) tvec(start)], ...
                    [iflux0  iflux0+(avgflux - err)*1000*deltat ...
                     iflux0+(avgflux + err)*1000*deltat iflux0], handles.ts(2).Color);
    handles.patch(2).FaceAlpha = patchAlpha;
    handles.patch(2).EdgeColor = 'none';
    handles.line(2) = plot([tvec(start) tvec(stop)], ...
                  [iflux0 iflux0+avgflux*1000*deltat], ...
                  '--', 'Color', handles.ts(2).Color);
    ylabel('% Volume transported');
    handles.hax(2).YTick = sort(unique([handles.hax(2).YTick(1) ...
                        ifluxvec(start) ifluxvec(stop) ...
                        handles.hax(2).YTick(end)]));
    handles.hax(2).YTickLabel{2} = '5%';
    handles.hax(2).YTickLabel{3} = '90%';
    handles.hax(2).YTickLabel{4} = '100%';

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

    axes(handles.hax(1)); beautify;
    handles.hax(1).YColor = handles.ts(1).Color;
    handles.hax(1).YLabel.Color = handles.ts(1).Color;
    linex(tvec([start stop]), [], [1 1 1]*0.85);
    LowerLines;

    axes(handles.hax(2)); beautify;
    handles.hax(2).YColor = handles.ts(2).Color;
    handles.hax(2).YLabel.Color = handles.ts(2).Color;

    handles.hisolabel = runs.add_isobathlabel(isobath);
    handles.hisolabel.Position(1:2) = [0.05 0.05];

    handles.hinterval = plot(tvec([start stop] + 2*[1 -1]), [0 0], ...
                             'LineWidth', 10, 'Color', [1 1 1] * 0.85);
    handles.htext = text(mean(handles.hinterval.XData), handles.hax(2).YTick(2), ...
                         '[t_{start}, t_{stop}]', 'Color', [1 1 1] * 0.6, ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle');
end