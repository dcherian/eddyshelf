function [] = plotAvgStreamer(runs, isobath)

    debug = 0;
    %pint = runs.csflux.off.slopewater.vertitrans(:,isobath,source);
    vmean = runs.streamer.off.vmean(:,:,isobath);
    xivec = runs.streamer.xivec;
    zvec = runs.streamer.zvec(:,isobath);

    hf = figure; maximize();
    insertAnnotation([runs.name '.plotAvgStreamer']);
    ax1 = subplot(3,3,[4 5 7 8]);
    pcolorcen(xivec, zvec, vmean');
    xlabel('X - X_{eddy} (km)'); ylabel('Z (m)');
    liney(-runs.bathy.hsb);
    linex(0);
    title([runs.name ' | mean streamer velocity | y/R = ' num2str(runs.csflux.ndloc(isobath))]);
    hcb = center_colorbar;
    hcb.Position(1) = 0.5;

    ax2 = subplot(3,3,[1 2]);
    hx = plot(xivec, runs.streamer.off.xprof(:,isobath));
    title('\int dz');
    hold on;
    linex(0); liney(0);
    ylabel('Flux (m^2/s)');

    ax3 = subplot(3,3,[6 9]);
    hz = plot(runs.streamer.off.zprof(:,isobath), zvec);
    title('Offshore transport (\int dx)');
    xlabel('Flux (m^2/s)');

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

        axes(ax2);
        hx.YData = hx.YData / max(actualx);
        plot(avg1(xivec), idealx./max(idealx), 'k-');
        ylabel('Flux / max flux');

        axes(ax3);
        hold on;
        hz.YData = hz.YData / max(pmean);
        plot(pint./max(pint), zvec);
        plot(ideal, zvec);
        legend('Mean', 'Integrated', 'Idealized', 'Location', 'SouthEast');
        xlabel('Flux / max flux');
    end

    linkaxes([ax1 ax2], 'x');
    linkaxes([ax1 ax3], 'y');
end