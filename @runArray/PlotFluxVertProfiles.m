function [handles] = PlotFluxVertProfiles(runArray)

    isobath = 1;

    figure; maximize;
    insertAnnotation('runArray.PlotFluxVertProfiles');
    handles.hax(1) = subplot(121); hold on;
    handles.hax(2) = subplot(122); hold on;

    phi = runArray.print_params('bathy.hsb./(V0./bathy.S_sh/sqrt(phys.N2))');
    runArray.sort(phi);
    phi = runArray.print_params('bathy.hsb./(V0./bathy.S_sh/sqrt(phys.N2))');

    corder_backup = runArray.sorted_colors;

    axes(handles.hax(1));
    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);
        run = runArray.array(ii);

        hsb = run.bathy.hsb;

        [start, stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
        offflux = run.csflux.off.slopezt(:,start:stop,1,1);
        zivec = run.csflux.vertbins(:,isobath);
        profile = trapz(run.csflux.time(start:stop)*86400, ...
                        offflux, 2);
        profile = profile./max(abs(profile));

        handles.hplt1(ff) = plot(profile, zivec./hsb, ...
                                 'DisplayName', ['\phi = ' ...
                            num2str(phi(ff), '%.2f')]);
    end

    chi = runArray.print_params(['(2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
                        'V0/Lz0/(bathy.S_sl*sqrt(phys.N2))']);
    runArray.sort(chi);
    chi = runArray.print_params(['(2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
                        'V0/Lz0/(bathy.S_sl*sqrt(phys.N2))']);

    axes(handles.hax(2));
    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);
        run = runArray.array(ii);

        hsb = run.bathy.hsb;

        [start, stop] = run.flux_tindices(run.csflux.on.slope(:,1,1), 0.05, 0.5);
        onflux = run.csflux.on.slopezt(:,start:stop,1,1);
        zivec = run.csflux.vertbins(:,isobath);
        profile = trapz(run.csflux.time(start:stop)*86400, ...
                        onflux, 2);
        profile = abs(profile./max(abs(profile)));

        handles.hplt2(ff) = plot(profile, zivec./hsb, ...
                                 'DisplayName', ['\chi = ' ...
                            num2str(chi(ff), '%.2f')]);
    end

    axes(handles.hax(1));
    set(gca, 'XAxisLocation', 'Top');
    liney(-1);%ylim([-1 0]);
    limx = xlim;
    xlim([0 limx(2)]);
    xlabel('a) Shelf water outflow (m^2/s)');
    ylabel('Z / H_{sb}');
    handles.hleg(1) = legend('Location', 'NorthWest');
    axis square;
    beautify;

    axes(handles.hax(2));
    set(gca, 'XAxisLocation', 'Top');
    liney(-1);
    limx = xlim;
    xlim([0 limx(2)]);
    xlabel('b) Eddy & slope water inflow (m^2/s)');
    handles.hax(2).YTickLabels = {};
    axis square;
    handles.hleg(2) = legend('Location', 'NorthWest');
    beautify;
    linkaxes(handles.hax, 'y');
    handles.hax(2).Position(1) = 0.5;
end
