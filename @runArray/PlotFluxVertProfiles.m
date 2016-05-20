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

    kk = 1;
    axes(handles.hax(1));
    for ii=1:length(runArray.array)
        run = runArray.array(ii);

        hsb = run.bathy.hsb;

        if run.params.misc.rdrg ~= 0
            continue;
        end

        [start, stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
        offflux = run.csflux.off.slopezt(:,start:stop,1,1);
        zivec = run.csflux.vertbins(:,isobath);
        profile = trapz(run.csflux.time(start:stop)*86400, ...
                        offflux, 2);
        profile = profile./max(abs(profile));

        handles.hplt1(kk) = plot(profile, zivec./hsb, ...
                                 'DisplayName', ['\phi = ' ...
                            num2str(phi(ii), '%.2f')]);
        if run.bathy.sl_shelf == 0
            handles.hplt1(kk).Color = [1 1 1]*0;
            hflat(1) = handles.hplt1(kk);
        end
        kk = kk+1;
    end

    %uistack(hflat, 'top');

    chi = runArray.print_params(['1./((2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
                        'V0/Lz0/(bathy.S_sl*sqrt(phys.N2)))']);
    runArray.sort(chi);
    chi = runArray.print_params(['1./((2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
                        'V0/Lz0/(bathy.S_sl*sqrt(phys.N2)))']);

    kk = 1;
    axes(handles.hax(2));
    for ii=1:length(runArray.array)
        run = runArray.array(ii);

        if run.params.misc.rdrg ~= 0
            continue;
        end

        hsb = run.bathy.hsb;

        [start, stop] = run.flux_tindices(run.csflux.on.slope(:,1,1), 0.05, 0.5);
        onflux = run.csflux.on.slopeztneg(:,start:stop,1,1);
        zivec = run.csflux.vertbins(:,isobath);
        profile = trapz(run.csflux.time(start:stop)*86400, ...
                        onflux, 2);
        profile = abs(profile./max(abs(profile)));

        handles.hplt2(kk) = plot(profile, zivec./hsb, ...
                                 'DisplayName', ['\chi = ' ...
                            num2str(chi(ii), '%.2f')]);
        if run.bathy.sl_shelf == 0
            handles.hplt2(kk).Color = [1 1 1]*0;
            hflat(2) = handles.hplt2(kk);
        end
        kk = kk+1;
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

    uistack(hflat(1), 'top');
    uistack(hflat(2), 'top');
end
