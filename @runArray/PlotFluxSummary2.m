function [handles] = PlotFluxSummary2(runArray, isobath, source, factor)

    if ~exist('isobath', 'var'), isobath = 2; end
    if ~exist('source', 'var'), source = isobath; end
    if ~exist('factor', 'var'), factor = 2; end

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    corder_backup = runArray.sorted_colors;

    handles.hfig = figure; maximize;

    handles.hax(1) = subplot(4,3,[1 2 4 5]); hold all;
    handles.hax(2) = subplot(4,3,[7 8 10 11]); hold all;
    handles.hax(3) = subplot(4,3,[9 12]); hold all;

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

        run = runArray.array(ii);
        names{ii} = runArray.getname(ii);

        tvec = run.csflux.time/86400;
        fluxvec = run.recalculateFlux(-factor*run.bathy.hsb, isobath);
        ifluxvec = cumtrapz(tvec*86400, fluxvec);
        zvec = run.csflux.vertbins(:,isobath);

        axes(handles.hax(1))
        handles.hflux(ff) = plot(tvec, fluxvec/1e3);

        axes(handles.hax(2))
        handles.henv(ff) = ...
            plot(tvec, (run.csflux.off.slopewater.envelope(:,isobath,source) ...
                        - run.bathy.xsb)/1000);

        % axes(handles.hax(3))
        % handles.hiflux(ff) = plot(tvec, ifluxvec);

        axes(handles.hax(3))
        handles.hprofile(ff) = plot( ...
            run.csflux.off.slopewater.vertitrans(:,isobath,source), ...
            zvec./abs(min(zvec)));

        linkprop([handles.hflux(ff) ...
                  handles.hprofile(ff) handles.henv(ff)], 'Color');
    end

    fs = [20 22 24];

    axes(handles.hax(1))
    handles.htxt(1) = text(0.05,0.85, 'a) Flux (mSv)', ...
                           'Units', 'Normalized');
    ylim([0 max(ylim)]);
    handles.hax(1).XTickLabel = {};
    [handles.hleg, handles.legobj] = columnlegend(2, names, 'FontSize', fs(1));
    handles.hleg.Position(1) = 0.62;
    beautify(fs);

    axes(handles.hax(2))
    handles.htxt(2) = text(0.05,0.15, {'b) Distance from shelfbreak of'; ...
                        'most onshore water parcel (km)'}, ...
                           'Units', 'Normalized');
    %handles.hax(2).XTickLabel = {};
    handles.hax(2).XAxisLocation = 'top';
    handles.hax(2).XTickLabel{1} = 't = 0 day';
    beautify(fs);
    %[handles.hl, handles.hltxt] = liney(0,'shelfbreak');
    %LowerLines;

    axes(handles.hax(3));
    ylabel('Z/H_{sb}');
    xlabel('c) \int Flux dx dt (m^2)');
    handles.hax(3).XAxisLocation = 'top';
    LowerLines;
    beautify(fs);
    xlim([0 max(xlim)]);

    linkaxes(handles.hax(1:2), 'x');

    linkprop(handles.htxt, 'FontSize');
    handles.htxt(1).FontSize = fs(1);
end