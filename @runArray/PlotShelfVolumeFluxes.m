function [handles] = PlotShelfVolumeFluxes(runArray)

    if isempty(runArray.filter)
        runArray.filter = [1:runArray.len];
    end

    ls = {'-';'--'; '-.'};

    gray = [1 1 1]*0.5;
    eddyeColor = runArray.array(1).eddyeColormap;
    colors.ShelfOff = runArray.array(1).shelfSlopeColor('light');
    colors.NonShelfOn = eddyeColor(end-3,:,:);

    colors.lin = linspecer(5);

    figure;
    handles.hax = packfig(1,2);

    for ii=1:length(runArray.filter)
        ff = runArray.filter(ii);
        run = runArray.array(ff);

        names{ii} = runArray.getname(ff);
        tvec = run.eddy.t;

        % handles.hshelfoff(ii) = ...
        %     plot(run.eddy.t, run.csflux.off.slope(:,1,1), ...
        %          'Color', colors.ShelfOff, 'LineStyle', ls{ii});

        axes(handles.hax(1)); hold on;
        handles.offsh(ii) = ...
            plot(run.eddy.t, run.csflux.sb.off.shelf, ...
                 'Color', 'k', 'LineStyle', ls{ii});
        handles.onsh(ii) = ...
            plot(run.eddy.t, run.csflux.sb.on.shelf, ...
                 'Color', 'k', 'LineStyle', ls{ii});
        handles.netsh(ii) = ...
            plot(run.eddy.t, run.csflux.off.slope(:,1,1), ...
                 'Color', gray, 'LineStyle', ls{ii});

        avgflux = run.calc_avgflux(run.csflux.off.slope(:,1,1));
        hline(1,ii) = liney(avgflux);
        hline(1,ii).LineStyle = ls{ii};

        axes(handles.hax(2)); hold on;
        handles.onnonsh(ii) = ...
            plot(run.eddy.t, run.csflux.sb.on.nonshelf, ...
                 'Color', 'k', 'LineStyle', ls{ii});
        handles.offnonsh(ii) = ...
            plot(run.eddy.t, run.csflux.sb.off.nonshelf, ...
                 'Color', 'k', 'LineStyle', ls{ii});
        handles.netnonsh(ii) = ...
            plot(run.eddy.t, run.csflux.on.slope(:,1,1), ...
                 'Color', gray, 'LineStyle', ls{ii});

        avgflux = run.calc_avgflux(run.csflux.on.slope(:,1,1));
        hline(2,ii) = liney(avgflux);
        hline(2,ii).LineStyle = ls{ii};

        % handles.halongiso(ii) = ...
        %     plot(run.eddy.t, run.volume.AlongIsobathSh(2,:), ...
        %          'Color', colors.lin(5,:), 'LineStyle', ls{ii});
    end

    linkaxes(handles.hax([2 1]), 'xy');
    handles.hax(1).Title.String = 'Shelf water';
    handles.hax(2).Title.String = 'Non-shelf water';

    handles.hax(1).YLabel.String = 'Volume Transport (m^3/s)';

    [handles.hleg, objh, outh, outm]  = ...
        legend(handles.hax(1), ...
               [handles.offsh(1) handles.offsh(2)], ...
               names{1}, names{2}, 'Location', 'SouthWest');

    common(handles.hax(1));
    common(handles.hax(2));

    %objh(7).LineStyle = '-';
    %objh(9).Color = 'k';
end

function [] = common(hax)

    axes(hax)
    hax.XAxisLocation = 'origin';
    beautify([12 13 14]);
    ylim([-1 1]*max(abs(ylim)))
    liney(0);
end