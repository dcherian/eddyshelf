function [handles] = PlotCrossSBVelocity(runs, surfvar, time)

    time = runs.process_time(time);

    opt.drawtrack = 0;
    opt.addvelquiver = 1;
    opt.quiverloc = 'surf';
    opt.nocolorbar = 0;

    clim = [-50 250];

    figure; maximize;
    insertAnnotation([runs.name '.PlotCrossSBVelocity']);

    % surface dye + surface velocity
    hax(1) = subplot(2,5, [1 2]);
    handles.hsurf = runs.animate_field(surfvar, hax(1), time, 1, opt);
    caxis(clim); colorbar('off');
    title('Surface  velocity vectors');
    htxt(1) = text(0.1, 0.1, 'a) surface dye', 'Units', 'Normalized');

    % xz density + v?
    opt = [];
    opt.rhocontours = 1;
    opt.eddy0 = 0;
    hax(2) = subplot(2,5,[3 4]); cla('reset');
    handlesxz = runs.PlotSingleXZSection('v', 1, time, opt, hax(2));
    for ii=1:2
        handlesxz.hline{ii}.delete;
        handlesxz.htext{ii}.delete;
    end
    handlesxz.htime.delete;
    handlesxz.hcb.Label.String = 'Cross-shelf velocity (m/s)';
    title('Density contours at shelfbreak');
    linkprop(hax(1:2), 'XLim');
    hax(2).Position(1) = 0.48;
    htxt(2) = text(0.1, 0.1, 'b)', 'Units', 'Normalized');

    hax(3) = subplot(2,5,[6 7]);
    handles.off = pcolorcen(runs.csflux.time/86400, runs.csflux.vertbins(:,1), ...
                            runs.csflux.off.slopezt(:,:,1,1));
    xlabel('Time (days)');
    ylabel('Z (m)');
    linex(runs.csflux.time(time)/86400, [], 'k');
    title('\int (v) x (shelf water mask) dx (m^2/s)');
    center_colorbar;
    htxt(3) = text(0.1, 0.1, 'c)', 'Units', 'Normalized');
    beautify;

    hax(4) = subplot(2,5,[8 9]);
    handles.on = pcolorcen(runs.csflux.time/86400, runs.csflux.vertbins(:,1), ...
                           runs.csflux.on.slopezt(:,:,1,1));
    linex(runs.csflux.time(time)/86400, [], 'k');
    xlabel('Time (days)');
    ylabel('');
    hax(4).YTickLabel = {};
    title({'\int (v) x (non-shelf water mask) dx (m^2/s)'});
    center_colorbar;
    htxt(4) = text(0.1, 0.1, 'd)', 'Units', 'Normalized');
    hax(4).Position(1) = 0.48;
    beautify;

    hax(5) = subplot(2,5,5);
    hax(5).Position(1) = 0.85;
    hax(5).Position(2) = 0.3;
    hold on;
    % offshore
    [start, stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1));
    offflux = runs.csflux.off.slopezt(:,start:stop,1,1);
    zivec = runs.csflux.vertbins(:,1);
    profile = trapz(runs.csflux.time(start:stop)*86400, ...
                    offflux, 2);
    profile = profile./max(abs(profile));
    plot(profile, zivec);

    % onshore
    [start, stop] = runs.flux_tindices(runs.csflux.on.slope(:,1,1), 0.2, 0.9);
    onflux = runs.csflux.on.slopezt(:,start:stop,1,1);
    zivec = runs.csflux.vertbins(:,1);
    profile = trapz(runs.csflux.time(start:stop)*86400, ...
                    onflux, 2);
    profile = abs(profile./max(abs(profile)));
    plot(profile, zivec);

    onflux = runs.csflux.on.slopeztneg(:,start:stop,1,1);
    zivec = runs.csflux.vertbins(:,1);
    profile = trapz(runs.csflux.time(start:stop)*86400, ...
                    onflux, 2);
    profile = abs(profile./max(abs(profile)));
    plot(profile, zivec);
    htxt(5) = text(0.1, 0.1, 'e)', 'Units', 'Normalized');

    handles.hleg = legend('offshore shelf water', 'onshore slope/eddy water', ['onshore slope/eddy water ' ...
                        '(negative only)']);

    ylim([-runs.bathy.hsb 0]);
    ylabel('Z (m)');
    hax(5).XAxisLocation = 'top';
    xlabel({'Normalized';  'cross-shelfbreak';  'transport'});
    beautify;

    handles.hleg.FontSize = 16;
    handles.hleg.Position(1) = 0.75;
    handles.hleg.Position(2) = 0.195;

    handles.htxt = htxt;
    handles.hax = hax;
    handles.xz = handlesxz;
end