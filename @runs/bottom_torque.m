function [] = bottom_torque(runs)

    tind = length(runs.eddy.t);
    tind = [tind-80 tind];
    %tind = [1 length(runs.eddy.t)];

    rho0 = runs.params.phys.rho0;
    g = runs.params.phys.g;
    beta = runs.params.phys.beta;

    dx = 1000; dy = 1000;

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
    di = 10;
    imnx = min(ixmin(:)) - di; imny = min(iymin(:)) - di;
    imxx = max(ixmax(:)) + di; imxy = max(iymax(:)) + di;

    volumer = {'x' imnx imxx; ...
              'y' imny imxy};
    volumeu = {'x' imnx-1 imxx; ...
               'y' imny imxy};
    volumev = {'x' imnx imxx; ...
               'y' imny imxy+1};

    % eddy center
    mx = runs.eddy.vor.cx(tind(1):tind(2));
    my = runs.eddy.vor.cy(tind(1):tind(2));
    imy = vecfind(runs.rgrid.yr(1,imny:imxy), my);

    if isempty(runs.zeta)
        runs.read_zeta;
    end
    % subsample to size(mask);
    % then subsample to region I'm interested in.
    zeta = runs.zeta(2:end-1,2:end-1, tind(1):tind(2));
    zeta = zeta(imnx:imxx, imny:imxy, :);

    % subsample mask and f
    mask = mask(imnx:imxx, imny:imxy, :);
    f = runs.rgrid.f(2:end-1,2:end-1)';
    f = f(imnx:imxx, imny:imxy);
    % f - f @ center of eddy
    % f = bsxfun(@minus, f, permute(f(1,imy),[3 1 2]));
    f = repmat(f, [1 1 size(mask,3)]);

    % subsample bottom slope
    slbot = diff(runs.rgrid.h',1,2)./diff(runs.rgrid.y_rho',1,2);
    slbot = repmat(slbot(imnx:imxx, imny:imxy), [1 1 size(mask,3)]);

    % now read density and eddye fields
    rho = dc_roms_read_data(runs.dir, 'rho', tind, volumer, [], ...
                            runs.rgrid, 'his', 'single');
    %eddye = dc_roms_read_data(runs.dir, runs.eddname, tind, volumer, [], ...
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

    % grid vectors - referenced at each time level to location of
    % eddy center
    xrmat = bsxfun(@minus, runs.rgrid.xr(imnx:imxx, imny:imxy), ...
                   permute(mx, [3 1 2]));
    yrmat = bsxfun(@minus, runs.rgrid.yr(imnx:imxx, imny:imxy), ...
                   permute(my, [3 1 2]));
    zrmat = runs.rgrid.z_r(:,2:end-1,2:end-1);
    zrmat = zrmat(:,imny:imxy, imnx:imxx);
    zwmat = runs.rgrid.z_w(:,2:end-1,2:end-1);
    zwmat = zwmat(:,imny:imxy, imnx:imxx);
    zumat = runs.rgrid.z_u(:,2:end-1,2:end-1);
    zumat = zumat(:,imny:imxy, imnx:imxx);
    zvmat = runs.rgrid.z_v(:,2:end-1,2:end-1);
    zvmat = zvmat(:,imny:imxy, imnx:imxx);
    dzmat = diff(permute(zwmat, [3 2 1]), 1, 3);

    % subtract out background density to get anomaly
    rho = bsxfun(@minus, rho, rback) .* eddye;
    maskstr = [maskstr ' + rho.*eddye'];

    % depth-integrate density anomaly field from surface to bottom
    tic;
    disp('integrating vertically');
    irho = nan(size(rho));
    frho = flipdim(rho, 3); % flipped to integrate from _surface_
                            % to bottom
    fzrmat = flipdim(zrmat, 1);
    for ii=1:size(rho, 1)
        for jj=1:size(rho,2)
            irho(ii,jj,:,:) = cumtrapz(fzrmat(:, jj, ii), ...
                                       frho(ii, jj, :, :), 3);
        end
    end
    toc;
    irho = flipdim(irho, 3);
    clear frho fzrmat

    % calculate bottom pressure (x,y,t)
    % note that in Flierl (1987) the 1/ρ0 is absorbed into the
    % pressure variable
    pres = bsxfun(@plus, g./rho0 .* irho, g.*permute(zeta,[1 2 4 3]));
    pbot = squeeze(pres(:,:,1,:));

    %%%%%%%%% now, angular momentum
    % depth averaged velocities (m/s)
    use_davg = 1;
    mom_budget = 0;
    if use_davg
        ubar = dc_roms_read_data(runs.dir, 'ubar', tind, volumer, [], runs.rgrid, ...
                                 'his', 'single');
        vbar = dc_roms_read_data(runs.dir, 'vbar', tind, volumer, [], runs.rgrid, ...
                                 'his', 'single');

        % convert to depth integrated velocities (m^2/s)
        U = bsxfun(@times, H, ubar);
        V = bsxfun(@times, H, vbar);

        % vertically integrated angular momentum
        %iam = 1/2 * beta .* (V .* xrmat - U .* yrmat); % if ψ ~ O(1/r²)
        iam = beta .* bsxfun(@times, U, yrmat); % if ψ ~ O(1/r)
    else
        % read depth dependent velocity fields and integrate
        u = avg1(dc_roms_read_data(runs.dir, 'u', tind, volumeu, [], ...
                              runs.rgrid, 'his', 'single'), 1);
        % read depth dependent velocity fields and integrate
        v = avg1(dc_roms_read_data(runs.dir, 'v', tind, volumev, [], ...
                              runs.rgrid, 'his', 'single'), 2);

        um = u; % .* eddye; % .* bsxfun(@times, eddye, permute(mask, [1 2 4 3]));
        vm = v; % .* eddye; % .* bsxfun(@times, eddye, permute(mask, [1 2 4 3]));
        pm = pres; % .* bsxfun(@times, eddye, permute(mask, [1 2 4 3]));
        maskstr = [maskstr ' + u.*eddye'];

        % depth-integrate quantities
        U = nan([size(um, 1) size(um, 2) size(um, 4)]);
        if mom_budget
            V = U; V2 = V; UV = V; U2 = U; P = U;
        end
        tic;
        for ii=1:size(um, 1)
            for jj=1:size(um, 2)
                U(ii,jj,:) = trapz(zumat(:,jj,ii), um(ii,jj,:,:), 3);
                if mom_budget
                    V(ii,jj,:) = trapz(zvmat(:,jj,ii), vm(ii,jj,:,:), 3);
                    UV(ii,jj,:) = trapz(zvmat(:,jj,ii), um(ii,jj,:,:) ...
                                        .* vm(ii,jj,:,:), 3);
                    V2(ii,jj,:) = trapz(zvmat(:,jj,ii), ...
                                        vm(ii,jj,:,:).^2, 3);
                    U2(ii,jj,:) = trapz(zumat(:,jj,ii), um(ii,jj,:,:).^2, 3);
                    P(ii,jj,:) = trapz(zrmat(:,jj,ii), pm(ii,jj,:,:),3);
                end
            end
        end
        toc;

        % vertically integrated angular momentum
        iam = beta .* yrmat .* U;

        % try depth integrated momentum budget
        if mom_budget
            tic;
            dPdx = nan([size(um,1)-1 size(um,2) size(um,4)]);
            dpresdx = diff(pres,1,1)./dx;
            for ii=1:size(dpresdx,1)
                for jj=1:size(dpresdx, 2)
                    dPdx(ii,jj,:) = trapz(zrmat(:,jj,ii), dpresdx(ii,jj,:,:), 3);
                end
            end
            toc;


            maskstr = [maskstr ' + P.*botmask'];
            dPdx = squeeze(trapz(trapz(dPdx .* (avg1(botmask) > 0),1),2)*dx*dy);

            dv2dy = squeeze(trapz(trapz(diff(V2,1,2)./dy,1)*dx,2)*dy);
            du2dx = squeeze(trapz(trapz(diff(U2,1,1)./dx,1)*dx,2)*dy);
            duvdx = squeeze(trapz(trapz(diff(UV,1,1)./dx,1)*dx,2)*dy);
            duvdy = squeeze(trapz(trapz(diff(UV,1,2)./dy,1)*dx,2)*dy);
            dpdx = squeeze(trapz(trapz(diff(P,1,1)./dx .* ...
                                       (avg1(botmask,1)>0),1)*dx)*dy);
            dpdy = squeeze(trapz(trapz(diff(P,1,2)./dy .* ...
                                       (avg1(botmask,2)>0),1)*dx)*dy);
            fu = squeeze(trapz(trapz(f.*U, 1)*dx,2)*dx);
            fv = squeeze(trapz(trapz(f.*V, 1)*dx,2)*dy);
            dUdt = diff(squeeze(trapz(trapz(U,1)*dx, 2)*dy))./86400;
            dVdt = diff(squeeze(trapz(trapz(V,1)*dx, 2)*dy))./86400;

            % bottom torques
            btq = squeeze(trapz(trapz(pbot .* botmask .* runs.bathy.sl_slope, 1)*dx, ...
                                2)*dy);
            figure;
            subplot(211); hold all;
            plot(dUdt);
            plot(du2dx); plot(duvdy);
            plot(fv);
            plot(dpdx);
            plot(dUdt + avg1(du2dx + duvdy - fv + dpdx));
            title('x-mom');
            legend('dudt', 'd/dx(u^2)', 'd/dy(uv)', 'fv', 'dPdx', 'total');
            subplot(212); hold all;
            plot(dVdt);
            plot(duvdx); plot(dv2dy);
            plot(fu); plot(dpdy); plot(btq);
            plot(dVdt + avg1(dv2dy + duvdx + fu + dpdy + btq));
            title('y-mom');
            legend('dvdt', 'd/dx(uv)', 'd/dy(v^2)', 'fu', 'dPdy', ['p_0 ' ...
                                's_{0y}'],'total');
        end
    end

    %%%%%%%%% Translation term
    %c = runs.eddy.cvx(tind(1):tind(2)) .* 1000/86400; % convert to m/s
    c = smooth(runs.eddy.mvx(tind(1):tind(2)), 10) .* 1000/86400; % convert to m/s

    % height anomaly for eddy is zeta
    h = bsxfun(@minus, zeta, mean(zeta, 2));

    iv = bsxfun(@times, bsxfun(@times, h, f), permute(c, [3 2 1]));
    %iv2 = bsxfun(@times, bsxfun(@times, irho, f), permute(c, [3 1 2]));
    %iv = runs.params.phys.f0 .* U;

    %%%%%%%%% mask?
    botmask = pbot > 0.1*max(pbot(:));
    mask_rho = 1; botmask; %irho < -1;
    mpbot = mask_rho .* pbot;
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
    bottom.pbtorque = P .* runs.bathy.sl_slope;
    bottom.betatorque = AM;
    bottom.transtorque = V;
    bottom.time = runs.eddy.t(tind(1):tind(2))*86400;
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
