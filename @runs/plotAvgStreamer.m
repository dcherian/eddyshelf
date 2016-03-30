function [handles] = plotAvgStreamer(runs, isobath)

    debug = 0;
    %pint = runs.csflux.off.slopewater.vertitrans(:,isobath,source);
    vmean = runs.streamer.vmean(:,:,isobath);
    xivec = runs.streamer.xivec;
    zvec = runs.streamer.zvec(:,isobath);
    xprof = squeeze(trapz(zvec, vmean, 2)); % full vel : x-profile
    [zrho,~] = runs.predict_zpeak(isobath, 'detect');
    if isobath == 1, zrho = []; end

    hf = figure; maximize();
    insertAnnotation([runs.name '.plotAvgStreamer']);
    ax(1) = subplot(3,3,[4 5 7 8]);
    handles.hfield = pcolorcen(xivec, zvec, vmean');
    shading interp;
    xlabel('Along-isobath, X - X_{eddy} (km)'); ylabel('Depth, Z (m)');
    if isobath ~= 1
        [hl,ht] = liney([-runs.bathy.hsb -runs.eddy.Lgauss(1) zrho], ...
                        {'h_{sb}'; 'L_z'; 'z_\rho'});
        for ii=1:length(hl)
            hl{ii}.XData(2) = ax(1).XTick(end-2);
            ht{ii}.VerticalAlignment = 'middle';
            %ht3{ii}.Units = 'normalized';
            ht{ii}.Position(1) = mean(ax(1).XTick(end-2));
            %ht3{ii}.Units = 'data';
        end
        hl{end+1} = linex(0);
    else
        hl = linex(0);
    end

    handles.isolabel = runs.add_isobathlabel(isobath);
    beautify;
    title([runs.name ' | mean streamer velocity (m/s)']);
    hcb = center_colorbar;
    hcb.Position(1) = 0.58;

    ax(2) = subplot(3,3,[1 2]);
    hx = plot(xivec, xprof);
    title('(a) Horizontal profile (\int (b) dz, m^2/s)');
    hold on;
    hl2(1) = linex(0); hl2(2) = liney(0);
    uistack(hl2, 'bottom');
    beautify;
    htxt(1) = text(0.1, 0.85, 'Offshore', 'Units', 'Normalized');
    htxt(2) = text(0.8, 0.15, 'Onshore', 'Units', 'Normalized');
    linkprop([handles.isolabel htxt], {'FontSize', 'Color'});
    htxt(1).Color = [1 1 1]*0.45;

    ax(3) = subplot(3,3,[6 9]);
    hz = plot(runs.streamer.off.zprof(:,isobath), zvec);
    [hl3,ht3] = liney([-runs.bathy.hsb -runs.eddy.Lgauss(1) zrho], ...
                      {'h_{sb}'; 'L_z'; 'z_\rho'});
    for ii=1:3
        ht3{ii}.VerticalAlignment = 'middle';
        %ht3{ii}.Units = 'normalized';
        ht3{ii}.Position(1) = 0.98 * max(xlim);
        %ht3{ii}.Units = 'data';
    end
    title({'(c) Vertical profile of'; 'offshore transport'; '(\int (b) dx, m^2/s)'});
    beautify;
    ax(3).XAxisLocation = 'top';
    if isobath == 1, xlim([0 max(xlim)*1.15]); end

    if debug
        error('debug plots not fixed yet.');
        L = runs.eddy.rhovor.lmaj(1)/1000;
        R = runs.csflux.R/1000;

        a = 3; Ln = L/3;
        yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
        y0oL = R/L * (1 - yoR); % y0/L - used in derivation
        ideal = runs.streamer_ideal_profile(isobath);
        idealx = trapz(zvec, ideal) *  ...
                 diff(exp(-abs(xivec'/Ln).^a))./diff(xivec'/Ln) ...
                 * exp(-y0oL.^2);

        axes(ax(2));
        hx.YData = hx.YData / max(actualx);
        plot(avg1(xivec), idealx./max(idealx), 'k-');
        ylabel('Flux / max flux');

        axes(ax(3));
        hold on;
        hz.YData = hz.YData / max(pmean);
        plot(pint./max(pint), zvec);
        plot(ideal, zvec);
        legend('Mean', 'Integrated', 'Idealized', 'Location', 'SouthEast');
        xlabel('Flux / max flux');
    end

    linkaxes(ax(1:2), 'x');
    linkaxes(ax([1 3]), 'y');

    linkprop(ax(1:2), 'XTick');
    linkprop(ax([1 3]), 'YTick');

    handles.ax = ax;
    handles.hx = hx;
    handles.hz = hz;
    handles.hcb = hcb;
    handles.htxt = htxt;

    handles.ax(2).XTickLabel = {};
    handles.ax(3).YTickLabel = {};
end