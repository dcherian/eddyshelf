function [] = bottom_torque(runs)

    tindices = [1 5 length(runs.eddy.t)];

    slab = Inf;
    [iend,tind,dt,nt,~] = roms_tindices(tindices, Inf, ...
                                        length(runs.eddy.t));

    rho0 = runs.params.phys.rho0;
    g = runs.params.phys.g;
    beta = runs.params.phys.beta;
    f0 = runs.params.phys.f0;

    % use some mask to determine edges of domain
    % that I want to analyze?
    use_mask = 0;
    if use_mask
        % eddy-based mask
        mask = runs.eddy.mask(:,:,tind(1):dt:tind(2));
        maskstr = 'sshmask';

        % vorticity mask
        %mask = runs.eddy.vormask(:,:,tind(1):dt:tind(2));
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

        di = 40;
        imnx = min(ixmin(:)) - di; imny = min(iymin(:)) - di;
        imxx = max(ixmax(:)) + di; imxy = max(iymax(:)) + di;

        mask = mask(imnx:imxx, imny:imxy, :);
    else
        sz = size(runs.rgrid.x_rho') - 4;
        imnx = 2; imxx = sz(1);
        imny = 2; imxy = sz(2);

        maskstr = '';
    end

    volumer = {'x' imnx imxx; ...
               'y' imny imxy};
    volumeu = {'x' imnx-1 imxx; ...
               'y' imny imxy};
    volumev = {'x' imnx imxx; ...
               'y' imny imxy+1};

    % eddy center
    %mx = runs.eddy.vor.cx(tind(1):dt:tind(2));
    %my = runs.eddy.vor.cy(tind(1):dt:tind(2));
    %imy = vecfind(runs.rgrid.yr(1,imny:imxy), my);

    % read free-surface
    if isempty(runs.zeta)
        runs.read_zeta;
    end
    % subsample to size(mask);
    % then subsample to region I'm interested in.
    zeta = runs.zeta(2:end-1,2:end-1, tind(1):dt:tind(2));
    zeta = zeta(imnx:imxx, imny:imxy, :);

    % subsample f
    f = runs.rgrid.f(2:end-1,2:end-1)';
    f = f(imnx:imxx, imny:imxy);
    % f - f @ center of eddy
    % f = bsxfun(@minus, f, permute(f(1,imy),[3 1 2]));
    f = repmat(f, [1 1 nt]);

    % subsample bathymetry
    H = runs.bathy.h(2:end-1, 2:end-1);
    H = H(imnx:imxx, imny:imxy);

    % subsample bottom slope
    slbot = diff(runs.rgrid.h',1,2)./diff(runs.rgrid.y_rho',1,2);
    slbot = repmat(slbot(imnx:imxx, imny:imxy), [1 1 nt]);

    % get background density field for initial time instant
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

    % now read density and eddye fields
    rho = dc_roms_read_data(runs.dir, 'rho', tindices, volumer, [], ...
                            runs.rgrid, 'his', 'single');
    %eddye = dc_roms_read_data(runs.dir, runs.eddname, tind, volumer, [], ...
    %                        runs.rgrid, 'his', 'single') > runs.eddy_thresh;

    % grid vectors and matrices
    xrmat = repmat(runs.rgrid.x_rho(imny:imxy, imnx:imxx)', [1 1 nt]);
    yrmat = repmat(runs.rgrid.y_rho(imny:imxy, imnx:imxx)', [1 1 nt]);

    zrmat = runs.rgrid.z_r(:,2:end-1,2:end-1);
    zrmat = zrmat(:,imny:imxy, imnx:imxx);
    zwmat = runs.rgrid.z_w(:,2:end-1,2:end-1);
    zwmat = zwmat(:,imny:imxy, imnx:imxx);
    dzmat = diff(permute(zwmat, [3 2 1]), 1, 3);
    xvec = runs.rgrid.x_rho(1,imnx:imxx);
    yvec = runs.rgrid.y_rho(imny:imxy);

    %%%%%%% first, bottom pressure

    % subtract out background density to get anomaly
    % see Flierl (1987)
    rho = bsxfun(@minus, rho, rback);

    % depth-integrate density anomaly field from surface to bottom
    % tic;
    % disp('integrating vertically');
    % irho = nan(size(rho));
    % frho = flipdim(rho, 3); % flipped to integrate from _surface_
    %                         % to bottom
    % fzrmat = flipdim(zrmat, 1);
    % for ii=1:size(rho, 1)
    %     for jj=1:size(rho,2)
    %         irho(ii,jj,:,:) = cumtrapz(fzrmat(:, jj, ii), ...
    %                                    frho(ii, jj, :, :), 3);
    %     end
    % end
    % toc;
    % irho = flipdim(irho, 3);
    % clear frho fzrmat

    % calculate bottom pressure (x,y,t)
    % note that in Flierl (1987) the 1/ρ0 is absorbed into the
    % pressure variable
    %pres = bsxfun(@plus, g./rho0 .* irho, g.*permute(zeta,[1 2 4
    %3]));
    % pres = -g/ρ0 ∫_{z}^{ζ} ρ dz
    % -1*dzmat to integrate from _surface_ to bottom cumulatively
    irho = flipdim(cumsum(flipdim(bsxfun(@times, rho, -1*dzmat),3),3),3);
    pres = -g./rho0 .* irho;
    % remove some more background signal
    pres = bsxfun(@minus, pres, pres(1,:,:,1));
    pbot = squeeze(pres(:,:,1,:));

    %%%%%%%%% now, angular momentum

    % depth averaged velocities (m/s)
    use_davg = 1;
    mom_budget = 0;
    if use_davg
        ubar = dc_roms_read_data(runs.dir, 'ubar', tindices, volumer, [], runs.rgrid, ...
                                 'his', 'single');
        vbar = dc_roms_read_data(runs.dir, 'vbar', tindices, volumer, [], runs.rgrid, ...
                                 'his', 'single');

        % convert to depth integrated velocities (m^2/s)
        U = bsxfun(@times, H, ubar);
        V = bsxfun(@times, H, vbar);

     else
        % read depth dependent velocity fields and integrate
        u = avg1(dc_roms_read_data(runs.dir, 'u', tindices, volumeu, [], ...
                              runs.rgrid, 'his', 'single'), 1);
        % read depth dependent velocity fields and integrate
        v = avg1(dc_roms_read_data(runs.dir, 'v', tindices, volumev, [], ...
                                   runs.rgrid, 'his', 'single'), 2);

        vormask = runs.eddy.vormask(imnx:imxx, imny:imxy, :);
        sshmask = runs.eddy.mask(imnx:imxx, imny:imxy, :);

        % find what density corresponds to 0 vorticity contour
        rhothreshvor = squeeze(nanmax(nanmax(rho(:,:,1) .* ...
                                             fillnan(vormask(:,:,1),0), ...
                                             [], 1), [], 2));
        rhothreshssh = squeeze(nanmax(nanmax(rho(:,:,1) .* ...
                                             fillnan(sshmask(:,:,1),0), ...
                                             [], 1), [], 2));

        % mask out velocities
        masked = bsxfun(@times, rho < rhothreshssh, ...
                        permute(sshmask(:,:,tind(1):dt:tind(2)), [1 2 4 3]));
        u = u .* masked;
        v = v .* masked;

        % depth-integrate quantities
        tic;
        U = squeeze(sum(bsxfun(@times, u, dzmat), 3));
        if mom_budget
             V = squeeze(sum(bsxfun(@times, v, dzmat), 3));
            UV = squeeze(sum(bsxfun(@times, u.*v, dzmat), 3));
            U2 = squeeze(sum(bsxfun(@times, u.^2, dzmat), 3));
            V2 = squeeze(sum(bsxfun(@times, v.^2, dzmat), 3));
             P = squeeze(sum(bsxfun(@times, pres, dzmat), 3));
        end
        toc;

        Y = max(runs.rgrid.y_rho(:));

        % try depth integrated momentum budget
        if mom_budget

            % pressure gradients
            dpdx = integrate(avg1(xvec), yvec, ...
                             bsxfun(@rdivide, diff(P,1,1), diff(xvec')));
            dpdy = integrate(xvec, avg1(yvec), ...
                             bsxfun(@rdivide, diff(P,1,2), diff(yvec)));

            % coriolis terms
            fv = integrate(xvec, yvec, f .* V);
            fu = integrate(xvec, yvec, f .* U);
            f0u = integrate(xvec, yvec, f0 .* U);
            byu = integrate(xvec, yvec, beta .* (yrmat-Y/2) .* U);
            f0v = integrate(xvec, yvec, f0 .* V);
            byv = integrate(xvec, yvec, beta .* (yrmat-Y/2) .* V);

            % non-linear terms
            dv2dy = integrate(xvec, avg1(yvec), ...
                              bsxfun(@rdivide, diff(V2,1,2), diff(yvec)));
            duvdx = integrate(avg1(xvec), yvec, ...
                              bsxfun(@rdivide, diff(UV,1,1), diff(xvec')));
            % tendency term - THIS IS A BAD ESTIMATE
            %dvdt = squeeze(trapz(trapz(diff(V,1,3)./86400,1),2));
            % bottom torque
            btq = integrate(xvec, yvec, pbot .* slbot);

            total = duvdx + dv2dy + fu + dpdy - btq;
            figure; hold all;
            plot(-1*f0u./total);
            plot(-1*byu./total);
            plot(dpdy./total);
            plot(duvdx./total);
            plot(dv2dy./total);
            plot(btq./total);
            legend('-f_0u','\beta yu', 'dpdy','duvdx','dv2dy', ...
                   'btq');

            time = runs.eddy.t(tind);
            save([runs.dir '/mombudget.mat'], 'dpdx', 'dpdy', 'fu', ...
                 'fv', 'f0u', 'byu', 'dv2dy', 'duvdx', 'btq', 'total', ...
                 'time');
        end
    end

    % vertically integrated angular momentum
    %iam = 1/2 * beta .* (V .* xrmat - U .* yrmat); % if ψ ~ O(1/r²)
    iam = f0 .* U + beta .* U .* yrmat; % if ψ ~ O(1/r)

    %%%%%%%%% Translation term
    %c = runs.eddy.cvx(tind(1):dt:tind(2)) .* 1000/86400; % convert to m/s
    c = smooth(runs.eddy.mvx(tind(1):dt:tind(2)), 10) .* 1000/86400; % convert to m/s

    % height anomaly for eddy is zeta
    h = bsxfun(@minus, zeta, mean(zeta, 2));

    iv = bsxfun(@times, bsxfun(@times, h, f), permute(c, [3 2 1]));
    %iv2 = bsxfun(@times, bsxfun(@times, irho, f), permute(c, [3 1 2]));
    %iv = runs.params.phys.f0 .* U;

    %%%%%%%%% mask?
    botmask = pbot < 0.1*min(pbot(:));
    mask_rho = botmask; %irho < -1;
    mpbot = mask_rho .* pbot .* slbot;
    miv = mask_rho .* iv;
    miam = mask_rho .* iam;

    clear V P AM
    %%%%%%%%% area-integrate
    for tt=1:size(iam,3)
        P(tt) = squeeze(trapz(yrmat(1,:,tt), ...
                              trapz(xrmat(:,1,tt), repnan(mpbot(:,:,tt),0), ...
                                    1), 2));
        AM(tt) = squeeze(trapz(yrmat(1,:,tt), ...
                               trapz(xrmat(:,1,tt), repnan(miam(:,:,tt),0), ...
                                                    1), 2));
        V(tt) = squeeze(trapz(yrmat(1,:,tt), ...
                              trapz(xrmat(:,1,tt), repnan(miv(:,:,tt),0), ...
                                             1), 2));
    end

    %%%%%%%%% Summarize
    bottom.pressure = P;
    bottom.angmom = AM;
    bottom.pbtorque = P;
    bottom.betatorque = AM;
    bottom.transtorque = V;
    bottom.time = runs.eddy.t(tind(1):dt:tind(2))*86400;
    bottom.maskstr = maskstr;

    % plots
    figure; hold all
    plot(bottom.time/86400, bottom.pbtorque);
    plot(bottom.time/86400, bottom.betatorque);
    plot(bottom.time/86400, bottom.transtorque);
    legend('\alpha \int\int P_{bot}', '\beta \int\int \Psi', ['c\' ...
                        'int\int fh'], 'Location', 'NorthWest');
    beautify;

    bottom.comment = ['(pressure, angmom) = volume integrated ' ...
                      'pressure, angular momentum | pbtorque = slope ' ...
                      '* pressure | betatorque = beta .* angmom'];

    bottom.hash = githash([mfilename('fullpath') '.m']);

    runs.bottom = bottom;
    save([runs.dir '/bottom.mat'], 'bottom', '-v7.3');
end

function [out] = integrate(xvec, yvec, in)
    out = squeeze(trapz(yvec, ...
                        trapz(xvec, in, 1), 2));
end


    % looks like (eddy.mx, eddy.my) isn't totally accurate, so
    % re-detect that.
    %xrmat = runs.rgrid.xr(imnx:imxx, imny:imxy);
    %yrmat = runs.rgrid.yr(imnx:imxx, imny:imxy);
    %clear mx my
    %for tt=1:size(zeta, 3)
    %    mzeta = mask(:,:,tt) .* zeta(:,:,tt);
    %    maxz = nanmax(nanmax(mzeta, [], 1), [], 2);
    %    ind = find(mzeta == maxz);
    %    [a,b] = ind2sub([size(mzeta,1) size(mzeta,2)], ind);
    %    mx(tt) = xrmat(a,b);
    %    my(tt) = yrmat(a,b);
    %end
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


        % check edge detection
        %for ind = 1:size(mask, 3)
        %    clf;
        %    pcolorcen(mask(:,:,ind)');
        %    linex([ixmin(ind) ixmax(ind)]);
        %    liney([iymin(ind) iymax(ind)]);
        %    title(num2str(ind));
        %    pause(1);
        %end

        % grid vectors - referenced at each time level to location of
    % eddy center
    %xrmat = bsxfun(@minus, runs.rgrid.xr(imnx:imxx, imny:imxy), ...
    %               permute(mx, [3 1 2]));
    %yrmat = bsxfun(@minus, runs.rgrid.yr(imnx:imxx, imny:imxy), ...
    %               permute(my, [3 1 2]));
