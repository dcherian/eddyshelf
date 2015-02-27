function [] = plot_eddytrack(runs)

    plotx = (runs.eddy.mx - runs.eddy.mx(1))/runs.rrdeep;
    ploty = (runs.eddy.my - runs.bathy.xsb)/runs.rrdeep;

    figure;
    insertAnnotation('runArray.plot_eddytrack');
    plot(plotx, ploty);
    hold all

    limy = ylim; ylim([-2 max(limy)]);

    % telescoping line
    liney((runs.rgrid.y_rho([runs.params.grid.iyp],1) - runs.bathy.xsb)/runs.rrdeep, ...
          'telescope');
    linex((runs.rgrid.x_rho(1, [runs.params.grid.ixn ...
                        runs.params.grid.ixp]) - runs.eddy.mx(1))/runs.rrdeep, ...
          'telescope');

    % sponge line
    sz = ceil(size(runs.sponge)./2);
    xx = sz(1); yy = sz(2);
    iy1 = find(runs.sponge(xx, 1:yy) == 1, 1, 'last');
    iy2 = find(runs.sponge(xx, yy:end) == 1, 1, 'first') + yy;
    ix1 = find(runs.sponge(1:xx, yy) == 1, 1, 'last');
    ix2 = find(runs.sponge(xx:end, yy) == 1, 1, 'first') + xx;

    liney((runs.rgrid.y_rho([iy1 iy2], 1) - runs.bathy.xsb)/runs.rrdeep, ...
          'sponge');
    linex((runs.rgrid.x_rho(1, [ix1 ix2]) - runs.eddy.mx(1))/runs.rrdeep, ...
          'sponge');

    liney((runs.bathy.xsl - runs.bathy.xsb)/runs.rrdeep, ...
          'slopebreak');

    [~,~,tind] = runs.locate_resistance;
    plot(plotx(tind), ploty(tind), 'x', 'MarkerSize', 12);

    liney(0);
end