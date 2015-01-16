function [] = bottom_torque(runs)

    tind = length(runs.eddy.t);
    tind = [tind-80 tind];

    rho0 = runs.params.phys.rho0;
    g = runs.params.phys.g;

    mask = runs.eddy.vormask(:,:,tind(1):tind(2));

    % indices of eddy extremes - based on mask
    indx = repmat([1:size(mask, 1)]', [1 size(mask,2)]);
    indy = repmat([1:size(mask, 2)], [size(mask,1) 1]);

    mask = fillnan(mask, 0);

    ixmax = squeeze(nanmax(nanmax(bsxfun(@times, mask, indx), [], ...
                                  1), [], 2));
    ixmin = squeeze(nanmin(nanmin(bsxfun(@times, mask, indx), [], ...
                                  1), [], 2));

    iymax = squeeze(nanmax(nanmax(bsxfun(@times, mask, indy), [], ...
                                  1), [], 2));
    iymin = squeeze(nanmin(nanmin(bsxfun(@times, mask, indy), [], ...
                                  1), [], 2));

    % check edge detection
    %for ind = 1:size(mask, 3)
    %    clf;
    %    pcolorcen(mask(:,:,ind)');
    %    linex([ixmin(ind) ixmax(ind)]);
    %    liney([iymin(ind) iymax(ind)]);
    %    title(num2str(ind));
    %    pause(1);
    %end

    %%%%%%% first, bottom pressure
    imnx = min(ixmin(:)); imny = min(iymin(:));
    imxx = max(ixmax(:)); imxy = max(iymax(:));

    volume = {'x' imnx imxx; ...
              'y' imny imxy};

    if isempty(runs.zeta)
        runs.read_zeta;
    end
    % subsample to size(mask);
    % then subsample to region I'm interested in.
    zeta = runs.zeta(2:end-1,2:end-1, tind(1):tind(2));
    zeta = zeta(imnx:imxx, imny:imxy, :);

    % subsample mask
    mask = mask(imnx:imxx, imny:imxy, :);

    % now read density and eddye fields
    rho = dc_roms_read_data(runs.dir, 'rho', tind, volume, [], ...
                            runs.rgrid, 'his', 'single');
    eddye = dc_roms_read_data(runs.dir, runs.eddname, tind, volume, [], ...
                            runs.rgrid, 'his', 'single') > runs.eddy_thresh;
    if runs.bathy.axis == 'y'
        % (y,z)
        rback = dc_roms_read_data(runs.dir, 'rho', [1 1], {'x' 1 1}, [], ...
                                  runs.rgrid, 'his', 'single');
        rback = rback(2:end-1, :);

        % subsample and make (x,y,z)
        rback = permute(rback(imny:imxy,:), [3 1 2]);
    else
        rback = dc_roms_read_data(runs.dir, 'rho', [1 1], {'y' Inf Inf}, [], ...
                                  runs.rgrid, 'his', 'single');
        error('not implemented for NS isobaths yet');
    end

    % bathymetry
    H = runs.bathy.h(2:end-1, 2:end-1);
    H = H(imnx:imxx, imny:imxy);

    % grid vectors - referenced at each time level to location of
    % eddy center
    xrmat = bsxfun(@minus, runs.rgrid.xr(imnx:imxx, imny:imxy), ...
                   permute(runs.eddy.mx(:,tind(1):tind(2)), [3 1 2]));
    yrmat = bsxfun(@minus, runs.rgrid.yr(imnx:imxx, imny:imxy), ...
                   permute(runs.eddy.my(:,tind(1):tind(2)), [3 1 2]));
    zrmat = runs.rgrid.z_r(:,2:end-1,2:end-1);
    zrmat = zrmat(:,imny:imxy, imnx:imxx);

    % subtract out density anomaly, and integrate from bottom to
    % surface
    rho = bsxfun(@minus, rho, rback) .* eddye;

    % depth-integrate density field
    tic;
    disp('integrating vertically');
    irho = nan([size(rho,1) size(rho,2) size(rho,4)]);
    for ii=1:size(rho, 1)
        for jj=1:size(rho,2)
            irho(ii,jj,:) = trapz(zrmat(:, jj, ii), rho(ii, jj, :, ...
                                                        :), 3);
        end
    end
    toc;

    % calculate bottom pressure (x,y,t)
    % note that in Flierl (1987) the 1/ρ0 is absorbed into the
    % pressure variable
    pbot = mask .* (rho0 .* g .* zeta +  irho .* g) ./ rho0;

    %%%%%%%%% now, angular momentum
    ubar = dc_roms_read_data(runs.dir, 'ubar', tind, volume, [], runs.rgrid, ...
                          'his', 'single');
    vbar = dc_roms_read_data(runs.dir, 'vbar', tind, volume, [], runs.rgrid, ...
                             'his', 'single');

    % vertically integrated angular momentum
    iam = mask .* bsxfun(@times, H, (vbar.*xrmat - ubar.*yrmat));


    %%%%%%%%% Summarize
    bottom.pressure = P;
    bottom.angmom = AM;
    bottom.pbtorque = P .* runs.bathy.sl_slope;
    bottom.betatorque = AM .* runs.params.phys.beta;

    bottom.comment = ['(pressure, angmom) = volume integrated ' ...
                      'pressure, angular momentum | pbtorque = slope ' ...
                      '* pressure | betatorque = beta .* angmom'];

    bottom.hash = githash;

    runs.bottom = bottom;
    save([runs.dir '/bottom.mat'], 'bottom', '-v7.3');
end