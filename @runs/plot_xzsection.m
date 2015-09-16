function [] = plot_xzsection(runs, days, hfig)

    if ~exist('hfig', 'var'), hfig = figure; end

    tindex = vecfind(runs.time, days*86400);
    yind = find_approx(runs.rgrid.y_rho(:,1), ...
                       runs.eddy.my(tindex)');

    L = runs.eddy.rhovor.dia(tindex)/2;

    zvec = runs.rgrid.z_r(:, yind+1, 1);
    xvec = (runs.rgrid.x_rho(1,2:end-1) - runs.eddy.mx(tindex)) ...
           ./ L;
    zmin = min(runs.rgrid.z_r(:,max(yind),1));

    % copied from csfluxes.m
    csvel = squeeze(avg1( ...
        dc_roms_read_data(runs.dir, 'v', tindex, ...
                          {runs.bathy.axis yind-1 yind}, [], runs.rgrid, ...
                          'his', 'single'), 2));
    csvel = csvel(2:end-1,:,:,:);

    % process cross-shelf dye
    csdye = dc_roms_read_data(runs.dir, runs.csdname, ...
                              tindex, {runs.bathy.axis yind yind}, ...
                              [], runs.rgrid, 'his', 'single');
    csdye = csdye(2:end-1,:,:);

    eddye = dc_roms_read_data(runs.dir, runs.eddname, ...
                              tindex, {runs.bathy.axis yind yind}, ...
                              [], runs.rgrid, 'his', 'single');
    eddye = eddye(2:end-1,:,:);

    rho = dc_roms_read_data(runs.dir, 'rho', ...
                            tindex, {runs.bathy.axis yind yind}, ...
                            [], runs.rgrid, 'his', 'single');
    rho = rho(2:end-1,:,:);

    rback = dc_roms_read_data(runs.dir, 'rho', ...
                              1, {runs.bathy.axis yind yind}, ...
                              [], runs.rgrid, 'his', 'single');
    rho = rho - rback(2:end-1,:,:);

    %mask = fillnan(bsxfun(@times, csdye < runs.csflux.x(isobath), ...
    %                      runs.csflux.westmask(:,tindex,isobath)),0)';

    figure(hfig)
    insertAnnotation([runs.name '.plot_xzsection']);

    ax(1) = subplot(221);
    [~,h1] = contourf(xvec, zvec, csvel', 20);
    hold on
    caxis([-1 1]);
    % contour(xvec, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
    runs.add_timelabel(tindex);
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title(['v (m/s) | ' runs.name]);
    % linex(xfrac);
    center_colorbar;
    xlim([-1 1]*max(abs(xlim))/2);

    ax(2) = subplot(222);
    pcolorcen(xvec, zvec, csdye'/1000); % .* mask
    hold on; shading interp
    % contour(xvec, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
    colorbar;
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title(['Cross-shelf dye (km) | ' runs.name]);
    runs.add_timelabel(tindex);
    linkaxes(ax, 'xy');
    % linex(xfrac);
    liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
    liney(-1 * runs.bathy.hsb, 'h_{sb}');

    ax(3) = subplot(223);
    pcolorcen(xvec, zvec, rho'); % .* mask
    hold on; shading interp;
    % contour(xvec, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
    colorbar;
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title(['Cross-shelf dye (km) | ' runs.name]);
    runs.add_timelabel(tindex);
    linkaxes(ax, 'xy');
    % linex(xfrac);
    liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
    liney(-1 * runs.bathy.hsb, 'h_{sb}');
    caxis([-0.05 0]);

    ax(4) = subplot(224);
    pcolorcen(xvec, zvec, eddye'); % .* mask
    hold on; shading interp;
    % contour(xvec, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
    colorbar;
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title(['Eddy dye (km) | ' runs.name]);
    runs.add_timelabel(tindex);
    linkaxes(ax, 'xy');
    % linex(xfrac);
    liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
    liney(-1 * runs.bathy.hsb, 'h_{sb}');
    caxis([0 1]);
end
