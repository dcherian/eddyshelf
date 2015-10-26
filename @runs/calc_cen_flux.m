function [] = calc_cen_flux(runs)

% calculate transport as fn of vertical depth - west of
% eddy only - f(z,t)
    if runs.bathy.axis == 'y'
        dx = 1./runs.rgrid.pm(1,2:end-1)';
    else
        dx = 1./runs.rgrid.pn(2:end-1,1);
    end

    t0 = 1;
    tinf = size(runs.csflux.slopext, 2);
    if runs.bathy.axis == 'y'
        spongemask = ~runs.sponge(2:end-1,1);
    else
        spongemask = ~runs.sponge(1,2:end-1)';
    end

    if runs.bathy.axis == 'y'
        cxi = runs.eddy.mx(t0:end);
        xvec = runs.rgrid.x_rho(1,2:end-1)';
    else
        cxi = runs.eddy.my(t0:end);
        xvec = runs.rgrid.y_rho(2:end-1,1);
    end

    if runs.bathy.axis == 'y'
        westmask = bsxfun(@times, ...
                          bsxfun(@lt, runs.eddy.xr(:,1), cxi(t0:tinf)), ...
                          spongemask);
    else
        westmask = bsxfun(@times, ...
                          bsxfun(@gt, runs.eddy.yr(1,:)', cxi(t0:tinf)), ...
                          spongemask);
    end
    eastmask = bsxfun(@times, 1 - westmask, spongemask);

    runs.csflux.cenflux.slopext = bsxfun(@times, runs.csflux.slopext, ...
                                         westmask);
    runs.csflux.cenflux.slope = squeeze( ...
        trapz(xvec, runs.csflux.cenflux.slopext, 1));
end
