function [] = plot_tracpy(runs)

    lonp = runs.tracpy.x/1000;
    latp = runs.tracpy.y/1000;

    xsb = runs.bathy.xsb/1000;
    %    inds = find(latp(end,:) > xsb & latp(1,:) < xsb);
    inds = runs.tracpy.reenter_inds;

    figure;
    subplot(2,1,1);
    runs.plot_bathy('contour', 'k');
    plot_map(lonp, latp, runs.tracpy.winds, runs.bathy.xsb);
    ylim([0 runs.rgrid.y_rho(end,1)/2000])
    xlim([0 runs.rgrid.x_rho(1,end)/1000]);
    title([runs.name ' | Day ' ...
           num2str(runs.tracpy.time(1)/86400) ' to ' ...
           num2str(runs.tracpy.time(end)/86400) ' | S_{sh} = ', ...
           num2str(runs.bathy.S_sh) ' | Streamer n = ' ...
           num2str(length(runs.tracpy.winds))]);
    beautify

    subplot(2,1,2);
    runs.plot_bathy('contour', 'k');
    plot_map(lonp, latp, runs.tracpy.reenter_inds, runs.bathy.xsb);
    ylim([0 runs.rgrid.y_rho(end,1)/2000])
    xlim([0 runs.rgrid.x_rho(1,end)/1000]);
    title([runs.name ' | Day ' ...
           num2str(runs.tracpy.time(1)/86400) ' to ' ...
           num2str(runs.tracpy.time(end)/86400) ' | S_{sh} = ', ...
           num2str(runs.bathy.S_sh) ' | Re-entry n = ' ...
           num2str(length(runs.tracpy.reenter_inds))]);
    beautify
end

function [] = plot_map(lonp, latp, inds, xsb);
    hold on;
    % all initial locations
    plot(lonp(1,:), latp(1,:), '.', 'Color', [158 202 225]/255);

    % all tracks
    for ii=1:length(inds)
        index = inds(ii);
        color = [1 1 1]*0.90;
        plot(lonp(:, index), latp(:, index), '-', 'Color', color);
    end

    % some randomly selected tracks
    for ii=1:100:length(inds)
        index = inds(ii);
        color = [1 1 1]*0.75;
        plot(lonp(:, index), latp(:, index), '-', 'Color', color);
    end

    % initial and final locations for drifters that crossed the
    % shelfbreak
    for ii=1:length(inds)
        index = inds(ii);
        plot(lonp(end, index), latp(end, index), '.', 'Color', [227 74 51]/255);
        plot(lonp(1, index), latp(1, index), '.', 'Color', [49,130,189]/ ...
             255);
    end
    axis image
    xlabel('X (km)');
    ylabel('Y (km)');

    miny = min(latp(1,inds));
    liney(miny, ['Width = ' num2str(xsb - miny) 'km'], 'k');
end