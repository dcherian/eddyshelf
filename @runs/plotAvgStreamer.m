function [handles] = plotAvgStreamer(runs, isobath)

    debug = 0;
    %pint = runs.csflux.off.slopewater.vertitrans(:,isobath,source);
    vmean = runs.streamer.vmean(:,:,isobath);
    xivec = runs.streamer.xivec;
    zvec = runs.streamer.zvec(:,isobath);

    hf = figure; maximize();
    insertAnnotation([runs.name '.plotAvgStreamer']);
    ax(1) = subplot(3,3,[4 5 7 8]);
    handles.hfield = pcolorcen(xivec, zvec, vmean');
    xlabel('X - X_{eddy} (km)'); ylabel('Z (m)');
    hl = liney([-runs.bathy.hsb -runs.eddy.Lgauss(1)]);
    hl{1}.XData(2) = ax(1).XTick(end-1);
    hl{2}.XData(2) = ax(1).XTick(end-1);
    hl{3} = linex(0);
    handles.isolabel = runs.add_isobathlabel(isobath);
    beautify;
    title([runs.name ' | mean streamer velocity (m/s)']);
    hcb = center_colorbar;
    hcb.Position(1) = 0.58;

    ax(2) = subplot(3,3,[1 2]);
    hx = plot(xivec, runs.streamer.xprof(:,isobath));
    title('\int dz');
    hold on;
    hl2(1) = linex(0); hl2(2) = liney(0);
    uistack(hl2, 'bottom');
    beautify;
    ylabel('Flux (m^2/s)');

    ax(3) = subplot(3,3,[6 9]);
    hz = plot(runs.streamer.off.zprof(:,isobath), zvec);
    title('Offshore transport (\int dx)');
    xlabel('Flux (m^2/s)');
    beautify;
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

    handles.ax(2).XTickLabel = {};
    handles.ax(3).YTickLabel = {};
end