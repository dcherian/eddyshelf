% calculate fluxes focusing on one radius on either side of eddy center
% avoids having to modify the clusterfuck that is csfluxes.m

function FluxOneRadius(runs, factor)

    ticstart = tic;
    if ~exist('factor', 'var'), factor=1.2; end

    kk = 1; % choose shelfbreak
    bathyax = 2;
    precision = 'single';
    ftype = 'his';

    [t0, tinf] = runs.flux_tindices(runs.csflux.off.slope(:,1,1));

    volr = {runs.bathy.axis runs.csflux.ix(kk) runs.csflux.ix(kk)};
    volv = {runs.bathy.axis runs.csflux.ix(kk)-1 runs.csflux.ix(kk)};

    % read along-shore section of cross-shore vel.
    % dimensions = (x/y , z , t )
    % average carefully to get values at RHO point
    csvel = squeeze(avg1(dc_roms_read_data(runs.dir, runs.csvelname, [t0 tinf], ...
                                           volv, [], runs.rgrid, ftype, ...
                                           precision), bathyax));
    csvel = csvel(2:end-1,:,:,:);

    % process cross-shelf dye
    csdye = dc_roms_read_data(runs.dir, runs.csdname, [t0 tinf], ...
                              volr, [], runs.rgrid, ftype, 'single');
    csdye = csdye(2:end-1,:,:);
    slopemask = csdye < runs.bathy.xsb;
    clear csdye

    xvec = runs.rgrid.x_rho(1, :);

    radius = [];
    npts = ceil(runs.eddy.fitx.L(1)/1000 * factor);
    radius.npts = npts;

    radius.off.shelf = ...
        nan([2*npts+1 runs.rgrid.N tinf-t0+1]);
    radius.on.nonshelf = ...
        nan([2*npts+1 runs.rgrid.N tinf-t0+1]);

    for tt=1:tinf-t0+1
        ix = find_approx(xvec, runs.eddy.mx(t0+tt), 1);
        range = (ix-npts):(ix+npts);
        if tt == 1
            radius.xvec = xvec(range)' - runs.eddy.mx;
        end

        % (x,z,t)
        radius.off.shelf(:,:,tt) = ...
            csvel(range, :, tt) .* slopemask(range, :, tt);
        radius.off.shelfpos(:,:,tt) = ...
            csvel(range, :, tt) .* slopemask(range, :, tt) ...
            .* (csvel(range, :, tt) > 0);

        radius.on.nonshelf(:,:,tt) = ...
            csvel(range, :, tt) .* ~slopemask(range, :, tt);

        radius.on.nonshelfneg(:,:,tt) = ...
            csvel(range, :, tt) .* ~slopemask(range, :, tt) ...
            .* (csvel(range, :, tt) < 0);
    end


    radius.comment = ['vertical profiles of cross-shelfbreak ' ...
                      'flux calculated within factor x radius of eddy center'];
    radius.factor = factor;
    hash = githash([mfilename('fullpath') '.m']);
    save([runs.dir '/radius.mat'], 'radius');

    runs.radius = radius;

    toc(ticstart);

end