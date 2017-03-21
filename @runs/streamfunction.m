% returns streamfunction at requested depth
%     [psi, iynan] = streamfunction(runs, depth)
% works only if topo is southern/northern coast.

function [psi, iynan] = streamfunction(runs, tt, depth)
    if ~exist('depth', 'var'), depth = 0; end

    if runs.bathy.axis == 'x'
        error('Streamfunction calc not setup for western coast');
    end

    if depth == 0
        runs.read_velsurf(tt);

        u = avg1(runs.usurf(:,2:end-1,tt), 1);
        v = avg1(runs.vsurf(2:end-1,:,tt), 2);

        iynan = 1;
    else
        [ugrid.xax, ugrid.yax, ugrid.zax,~,~,~] = ...
            dc_roms_var_grid(runs.rgrid,'u');
        [vgrid.xax, vgrid.yax, vgrid.zax,~,~,~] = ...
            dc_roms_var_grid(runs.rgrid,'v');

        % get z-slices of (u,v)
        u = dc_roms_zslice_var(...
            dc_roms_read_data(runs.dir, 'u', tt, ...
                              {}, [], runs.rgrid), ...
            depth, ugrid);
        v = dc_roms_zslice_var(...
            dc_roms_read_data(runs.dir, 'v', tt, ...
                              {}, [], runs.rgrid), ...
            depth, vgrid);

        u = avg1(u(:,2:end-1), 1);
        v = avg1(v(2:end-1,:), 2);

        % exclude nans for psi computation
        iynan = find(isnan(v(1,:)) == 1, 1, 'last') + 1;
    end

    xr = runs.rgrid.x_rho(1, 2:end-1)';
    yr = runs.rgrid.y_rho(2:end-1, 1)';

    psi = nan([size(runs.eddy.mask, 1) ...
               size(runs.eddy.mask, 2)]);

    [~, psi(:,iynan:end)] = ...
        flowfun(xr(:,1), yr(1, iynan:end), ...
                u(:,iynan:end), v(:,iynan:end), 'psi');
end