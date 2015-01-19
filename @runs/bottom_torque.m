function [] = bottom_torque(runs)

    tind = length(runs.eddy.t);
    tind = [tind-80 tind];
    %tind = [1 length(runs.eddy.t)];

    rho0 = runs.params.phys.rho0;
    g = runs.params.phys.g;

    % eddy-based mask
    mask = runs.eddy.mask(:,:,tind(1):tind(2));
    maskstr = 'sshmask';

    % vorticity mask
    %mask = runs.eddy.vormask(:,:,tind(1):tind(2));
    %maskstr = 'vormask';

    % topography based mask
    %mask = (runs.rgrid.y_rho(2:end-1,2:end-1)' > runs.bathy.xsb) & ...
    %       (runs.rgrid.y_rho(2:end-1,2:end-1)' < runs.bathy.xsl);
    %mask = mask .* ~runs.sponge(2:end-1,2:end-1);
    % mask = 'slopemask';

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
    %eddye = dc_roms_read_data(runs.dir, runs.eddname, tind, volume, [], ...
    %                        runs.rgrid, 'his', 'single') > runs.eddy_thresh;
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

    % looks like (eddy.mx, eddy.my) isn't totally accurate, so
    % re-detect that.
    xrmat = runs.rgrid.xr(imnx:imxx, imny:imxy);
    yrmat = runs.rgrid.yr(imnx:imxx, imny:imxy);
    clear mx my
    for tt=1:size(zeta, 3)
        mzeta = mask(:,:,tt) .* zeta(:,:,tt);
        maxz = nanmax(nanmax(mzeta, [], 1), [], 2);
        ind = find(mzeta == maxz);
        [a,b] = ind2sub([size(mzeta,1) size(mzeta,2)], ind);
        mx(tt) = xrmat(a,b);
        my(tt) = yrmat(a,b);
    end

    % grid vectors - referenced at each time level to location of
    % eddy center
    xrmat = bsxfun(@minus, runs.rgrid.xr(imnx:imxx, imny:imxy), ...
                   permute(mx, [3 1 2]));
    yrmat = bsxfun(@minus, runs.rgrid.yr(imnx:imxx, imny:imxy), ...
                   permute(my, [3 1 2]));
    zrmat = runs.rgrid.z_r(:,2:end-1,2:end-1);
    zrmat = zrmat(:,imny:imxy, imnx:imxx);

    % subtract out density anomaly, and integrate from bottom to
    % surface
    rho = bsxfun(@minus, rho, rback);

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
    pbot = (rho0 .* g .* zeta +  irho .* g) ./ rho0;

    %%%%%%%%% now, angular momentum
    % depth averaged velocities (m/s)
    ubar = dc_roms_read_data(runs.dir, 'ubar', tind, volume, [], runs.rgrid, ...
                          'his', 'single');
    vbar = dc_roms_read_data(runs.dir, 'vbar', tind, volume, [], runs.rgrid, ...
                             'his', 'single');

    % convert to depth integrated velocities (m^2/s)
    U = bsxfun(@times, H, ubar);
    V = bsxfun(@times, H, vbar);

    % vertically integrated angular momentum
    %iam = 1/2 .* (V .* xrmat - U .* yrmat); % if ψ ~ O(1/r²)
    iam = bsxfun(@times, U, yrmat); % if ψ ~ O(1/r)

    % debug plots
    %tt = 20;
    %mzeta = mask .* zeta;
    %for tt =1:size(mzeta,3)
    %    clf;
    %    contourf(xrmat(:,:,tt), yrmat(:,:,tt), mzeta(:,:,tt), 60);
    %    hold on;
    %    plot(runs.eddy.cx(tind(1)+tt) - mx(tt), ...
    %         runs.eddy.cy(tind(1)+tt) - my(tt), 'k*', 'MarkerSize', 16);
    %    shading flat;
    %    linex(0); liney(0);
    %    pause(0.5);
    %end

    %%%%%%%%% Translation term
    %c = runs.eddy.cvx(tind(1):tind(2)) .* 1000/86400; % convert to m/s
    c = smooth(runs.eddy.mvx(tind(1):tind(2)), 10) .* 1000/86400; % convert to m/s
    f = runs.rgrid.f(2:end-1,2:end-1)';
    f = f(imnx:imxx, imny:imxy);

    % height anomaly for eddy is zeta
    h = bsxfun(@minus, zeta, mean(zeta, 2));

    iv = bsxfun(@times, bsxfun(@times, h, f), permute(c, [3 2 1]));
    %iv2 = bsxfun(@times, bsxfun(@times, irho, f), permute(c, [3 1 2]));
    %iv = runs.params.phys.f0 .* U;

    %%%%%%%%% mask?
    pbot = mask .* pbot;
    iv = mask .* iv;
    iam = mask .* iam;

    clear V P AM
    %%%%%%%%% area-integrate
    for tt=1:size(iam,3)
        P(tt) = squeeze(trapz(yrmat(1,:,tt), trapz(xrmat(:,1,tt), repnan(pbot(:,:,tt),0), ...
                                            1), 2));
        AM(tt) = squeeze(trapz(yrmat(1,:,tt), trapz(xrmat(:,1,tt), repnan(iam(:,:,tt),0), ...
                                                    1), 2));
        V(tt) = squeeze(trapz(yrmat(1,:,tt), trapz(xrmat(:,1,tt), repnan(iv(:,:,tt),0), ...
                                             1), 2));
    end

    %%%%%%%%% Summarize
    bottom.pressure = P;
    bottom.angmom = AM;
    bottom.pbtorque = P .* runs.bathy.sl_slope;
    bottom.betatorque = AM .* runs.params.phys.beta;
    bottom.transtorque = V;
    bottom.time = runs.eddy.t(tind(1):tind(2))*86400;
    bottom.maskstr = maskstr;

    % plots
    figure; hold all
    plot(bottom.pbtorque);
    plot(bottom.betatorque);
    plot(bottom.transtorque);
    legend('\alpha \int\int P_{bot}', '\beta \int\int \Psi', ['c\' ...
                        'int\int fh']);

    bottom.comment = ['(pressure, angmom) = volume integrated ' ...
                      'pressure, angular momentum | pbtorque = slope ' ...
                      '* pressure | betatorque = beta .* angmom'];

    bottom.hash = githash;

    runs.bottom = bottom;
    save([runs.dir '/bottom.mat'], 'bottom', '-v7.3');
end