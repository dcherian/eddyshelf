function [] = plot_eddytrack(runs)

    if runs.bathy.axis == 'y'
        x0 = runs.eddy.mx(1);
        y0 = runs.bathy.xsb;
    else
        x0 = runs.bathy.xsb;
        y0 = runs.eddy.my(1);
    end

    xnorm = runs.eddy.vor.dia(1)/2;
    ynorm = runs.eddy.vor.dia(1)/2;

    plotx = (runs.eddy.mx - x0)/xnorm;
    ploty = (runs.eddy.my - y0)/ynorm;

    figure;
    insertAnnotation('runs.plot_eddytrack');
    plot(plotx, ploty);
    xlabel('Along isobath distance / Initial radius');
    ylabel('Cross isobath dist / Initial radius');
    hold all

    %limy = ylim; ylim([-2 max(limy)]);

    % telescoping line
    %liney((runs.rgrid.y_rho([runs.params.grid.iyp],1) - runs.bathy.xsb)/ynorm, ...
    %      'telescope');
    %linex((runs.rgrid.x_rho(1, [runs.params.grid.ixn ...
    %                    runs.params.grid.ixp]) - runs.eddy.mx(1))/xnorm, ...
    %      'telescope');

    % sponge line
    sz = ceil(size(runs.sponge)./2);
    xx = sz(1); yy = sz(2);
    iy1 = find(runs.sponge(xx, 1:yy) == 1, 1, 'last');
    iy2 = find(runs.sponge(xx, yy:end) == 1, 1, 'first') + yy;
    ix1 = find(runs.sponge(1:xx, yy) == 1, 1, 'last');
    ix2 = find(runs.sponge(xx:end, yy) == 1, 1, 'first') + xx;

    liney((runs.rgrid.y_rho([iy1 iy2], 1) - y0)/ynorm, ...
          'sponge');
    linex((runs.rgrid.x_rho(1, [ix1 ix2]) - x0)/xnorm, ...
          'sponge');

    if runs.bathy.axis == 'y'
        liney((runs.bathy.xsl - y0)/ynorm, ...
              'slopebreak');
    else
        linex((runs.bathy.xsl - x0)/xnorm, ...
              'slopebreak');
    end

    %[~,~,tind] = runs.locate_resistance;
    runs.fit_traj(1);
    tind = runs.traj.tind;
    plot(plotx(tind), ploty(tind), 'x', 'MarkerSize', 12);

    liney(0);
    %xlim([0 max(xlim)]);
end