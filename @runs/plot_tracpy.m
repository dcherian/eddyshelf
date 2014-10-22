function [] = plot_tracpy(runs)

    lonp = runs.tracpy.x/1000;
    latp = runs.tracpy.y/1000;

    xsb = runs.bathy.xsb/1000;
    inds = find(latp(end,:) > xsb & latp(1,:) < xsb);

    figure;
    hold on;
    runs.plot_bathy('contour', 'k');

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
    ylim([0 runs.rgrid.y_rho(end,1)/1000])
    xlim([0 runs.rgrid.x_rho(1,end)/1000]);
    xlabel('X (km)');
    ylabel('Y (km)');

    miny = min(latp(1,inds));
    liney(miny, ['Width = ' num2str(xsb - miny) 'km'], 'k');

    title([runs.name ' | Day ' ...
           num2str(runs.tracpy.time(1)/86400) ' to ' ...
           num2str(runs.tracpy.time(end)/86400) ' | S_{sh} = ', ...
           num2str(runs.bathy.S_sh)]);


    beautify
end