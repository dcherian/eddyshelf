function [] = bottom_torque(runs)

    ticstart = tic;
    tindices = [1 length(runs.eddy.t)];

    pcrit = 0.2;
    amcrit = 0.2;

    flags.subtract_edge = 1;
    flags.subtract_mean = 0;
    flags.use_time_varying_dz = 1;

    assert(flags.subtract_edge ~= flags.subtract_mean);

    % deprecated
    flags.use_prsgrd = 0;
    flags.calc_angmom = 0;
    flags.use_davg = 0;
    flags.mom_budget = 0;
    flags.use_thermal_wind = 0;
    % use some mask to determine edges of domain
    % that I want to analyze?
    flags.use_mask = 0;
    flags.use_masked = 0;

    slab = 10; % 5 at a time
    [iend,tind,dt,nt,~] = roms_tindices(tindices, slab, ...
                                        length(runs.eddy.t));

    rho0 = runs.params.phys.rho0;
    g = runs.params.phys.g;
    beta = runs.params.phys.beta;
    f0 = runs.params.phys.f0;

    if flags.use_mask
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
        imnx = runs.spng.sx1+2; imxx = runs.spng.sx2-2;
        imny = 2; imxy = runs.spng.sy2-2;

        maskstr = '';
    end

    volumer = {'x' imnx imxx; ...
               'y' imny imxy};
    volumeu = {'x' imnx-1 imxx; ...
               'y' imny imxy};
    volumev = {'x' imnx imxx; ...
               'y' imny-1 imxy};

    % grid vectors and matrices
    xrmat = repmat(runs.rgrid.x_rho(imny:imxy, imnx:imxx)', [1 1]);
    yrmat = repmat(runs.rgrid.y_rho(imny:imxy, imnx:imxx)', [1 1]);

    xvec = xrmat(:,1,1);
    yvec = yrmat(1,:,1);

    %xrmat = bsxfun(@minus, xrmat, permute(runs.eddy.mx, [3 1 2]));
    %yrmat = bsxfun(@minus, yrmat, permute(runs.eddy.my, [3 1 2]));

    % eddy center
    mx = runs.eddy.vor.cx(tind(1):dt:tind(2));
    my = runs.eddy.vor.cy(tind(1):dt:tind(2));
    imx = vecfind(runs.rgrid.xr(imnx:imxx,1), mx);
    imy = vecfind(runs.rgrid.yr(1,imny:imxy), my);

    % read free-surface to be sure I'm not screwing up.
    zeta = dc_roms_read_data(runs.dir, 'zeta', [tind(1) tind(2)], ...
                             volumer, [], runs.rgrid);

    % subsample f
    f = single(repmat(runs.rgrid.f(imny:imxy, imnx:imxx)', [1 1]));
    % f - f @ center of eddy
    % f = bsxfun(@minus, f, permute(f(1,imy),[3 1 2]));
    % This is so that I don't have trouble finding out the
    % reference latitude
    %bymat = single(f - f0);

    % subsample bathymetry
    H = runs.bathy.h(imnx:imxx, imny:imxy);

    % subsample bottom slope
    slbot = diff(runs.rgrid.h',1,2)./diff(runs.rgrid.y_rho',1,2);
    slbot = single(slbot .* (slbot > 0.95 * runs.bathy.sl_slope));
    slbot = slbot(imnx:imxx, imny:imxy);

    %vormask = runs.eddy.vormask(imnx-1:imxx-1, imny-1:imxy-1, :);
    %sshmask = runs.eddy.mask(imnx-1:imxx-1, imny-1:imxy-1, :);

    if flags.use_masked
        if ~isfield(runs.eddy, 'drhothreshssh')
            % find what density corresponds to 0 vorticity contour
            rhothreshvor = squeeze(nanmax(nanmax(rho(:,:,1) .* ...
                                                 fillnan(vormask(:,:,1),0), ...
                                                 [], 1), [], 2));
            rhothreshssh = squeeze(nanmax(nanmax(rho(:,:,1) .* ...
                                                 fillnan(sshmask(:,:,1),0), ...
                                                 [], 1), [], ...
                                          2));
        else
            rhothreshssh = runs.eddy.drhothreshssh;
            rhothreshvor = runs.eddy.drhothresh(1);
        end
    end

    % % get background density field for initial time instant
    % if runs.bathy.axis == 'y'
    %     % (y,z)
    %     rback = dc_roms_read_data(runs.dir, 'rho', [1 1], ...
    %                               {'x' 1 1; 'y' imny imxy}, [], ...
    %                               runs.rgrid, 'his') + 1000;

    %     % make (x,y,z)
    %     rback = permute(rback, [3 1 2]);
    % else
    %     rback = dc_roms_read_data(runs.dir, 'rho', [1 1], {'y' Inf Inf}, [], ...
    %                               runs.rgrid, 'his', 'single');
    %     error('not implemented for NS isobaths yet');
    % end

    % dzmat0 = diff(set_depth(2,4,runs.rgrid.theta_s,runs.rgrid.theta_b, ...
    %                      runs.rgrid.Tcline,runs.rgrid.N,5,H,...
    %                      zeta(:,:,1), 0), 1, 3);
    % % pressure due to background stratification = pstrat(y)
    % irback = bsxfun(@plus, rho0 .* zeta(1,:,1), ...
    %                 sum( (rback-rho0) .* dzmat0(1,:,:), 3));
    % pstrat = g./rho0 .* irback;
    % clear dzmat0;

    % read data from start
    pbot = single(nan(size(zeta)));
    AM = pbot; ipres = pbot;
    masku = logical(zeros(size(zeta)));
    maskp = masku;

    for i=0:iend-1
        disp(['==== Iteration : ' num2str(i+1) '/' num2str(iend)  ...
                   ' ====']);

        tstart = tindices(1) + i*slab*dt;
        tend = tindices(1) -1 + (i+1)*slab*dt;

        % now read density and eddye fields
        rho = dc_roms_read_data(runs.dir, 'rho', [tstart tend], volumer, [], ...
                                runs.rgrid, 'his') + 1000;

        tsave = (1+i*slab) + (0:size(rho,4)-1);
        assert(size(rho,4) == length(tsave));

        % pretty certain that this is correct. zwmat equals zeta at
        % surface and H at bottom.
        if flags.use_time_varying_dz
            tic;
            disp('Calculating time varying dz');
            zwmat = (nan(size(rho) + [0 0 1 0]));
            if flags.use_prsgrd
                zrmat = (nan(size(rho)));
            end
            for tt=1:size(rho,4)
                zwmat(:,:,:,tt) = (( ...
                    set_depth(2,4,runs.rgrid.theta_s,runs.rgrid.theta_b, ...
                                runs.rgrid.Tcline,runs.rgrid.N,5,H,...
                                zeta(:,:,tsave(tt)), 0)));
                if flags.use_prsgrd
                    zrmat(:,:,:,tt) = (( ...
                        set_depth(2,4,runs.rgrid.theta_s,runs.rgrid.theta_b, ...
                                    runs.rgrid.Tcline,runs.rgrid.N,1,H,...
                                    zeta(:,:,tsave(tt)), 0)));
                end
            end
            toc;
        end

         % U is provided by here
        %%%%%%% first, bottom pressure
        %irho = squeeze(sum(bsxfun(@times, rho, diff(zwmat,1,3)), 3));
        % avoid some roundoff errors?
        irhofull = bsxfun(@plus, rho0 .* permute(zeta(:,:,tsave), [1 2 4 3]), ...
                          flipdim(cumsum(flipdim(bsxfun(@times, rho-rho0, ...
                                                        diff(zwmat,1,3)), 3),3),3));
        % full pressure field
        %pres = g./rho0 .* bsxfun(@minus, irhofull,
        %irhofull(1,:,:,1));
        pres = g./rho0 .* irhofull;

        % removing mean zeta changes pbot by 1e-9 only.
        % using p_η = gη -> p_η η_y = gη η_y ~ O(1e-8), so not
        % much difference.
        % trapezoidal integration makes no difference
        %irhotrap = squeeze(sum(dzmat .* avg1(rho,3),3));
        %%%%%%%%% now, angular momentum
        ubar = dc_roms_read_data(runs.dir, 'ubar', [tstart tend], ...
                                 volumer, [], runs.rgrid, 'his', 'single');

        vbar = dc_roms_read_data(runs.dir, 'vbar', [tstart tend], ...
                                 volumer, [], runs.rgrid, 'his', 'single');

        U = ubar .* bsxfun(@plus, zeta(:,:,tsave), H);
        V = vbar .* bsxfun(@plus, zeta(:,:,tsave), H);

        % u = avg1(dc_roms_read_data(runs.dir, 'u', [tstart tend], ...
        %                       volumeu, [], runs.rgrid, 'his', ...
        %                       'single'),1);
        % v = avg1(dc_roms_read_data(runs.dir, 'v', [tstart tend], ...
        %                       volumev, [], runs.rgrid, 'his', ...
        %                       'single'),2);

        % U2 = squeeze(sum(u.*u.*diff(zwmat,1,3),3));
        % V2 = squeeze(sum(v.*v.*diff(zwmat,1,3),3));
        % UV = squeeze(sum(u.*v.*diff(zwmat,1,3),3));

        clear ubar vbar u v

        [~, AM(:,:,tsave)] = flowfun(xvec, yvec, U, V);

        AM = single(AM);
        if i == 0
            amcrit = amcrit .* squeeze(max(max(AM(:,:,1),[],1),[], ...
                                           2));
        end

        % get proper pressure & velocity regions
        masku(:,:,tsave) = find_mask(AM(:,:,tsave), amcrit, imx(tsave), ...
                                     imy(tsave));

        % duvdx(tsave) = integrate(avg1(xvec),yvec, ...
        %                   bsxfun(@rdivide, ...
        %                          diff(UV.*masku(:,:,tsave),1,1), ...
        %                          diff(xvec)));

        % duvdy(tsave) = integrate(xvec,avg1(yvec), ...
        %                   bsxfun(@rdivide, ...
        %                          diff(UV.*masku(:,:,tsave),1,2), ...
        %                          diff(yvec)));

        % du2dx(tsave) = integrate(avg1(xvec),yvec, ...
        %                   bsxfun(@rdivide, ...
        %                          diff(U2.*masku(:,:,tsave),1,1), ...
        %                          diff(xvec)));

        % dv2dy(tsave) = integrate(xvec,avg1(yvec), ...
        %                   bsxfun(@rdivide, ...
        %                          diff(V2.*masku(:,:,tsave),1,2), ...
        %                          diff(yvec)));

        % somehow remove background gradient post-boundary layer adjustment
        if flags.subtract_edge % remove western / eastern edge signal
            pres = bsxfun(@minus, pres, ...
                          (pres(1,:,:,:) + pres(end,:,:,:))/2);
        end
        if flags.subtract_mean
            % determine mean outside the AM contour
            pmean = nanmean(fillnan(bsxfun(@times, pres, ...
                                           permute(~masku(:,:,tsave), [1 2 4 3])), 0), 1);
            pres = bsxfun(@minus, pres, pmean);
        end

        % bottom pressure
        pbot(:,:,tsave) = single(squeeze(pres(:,:,1,:)));
        % integrated pressure
        ipres(:,:,tsave) = squeeze(sum(pres .* diff(zwmat,1,3), ...
                                       3));

        if i == 0
            pcrit = pcrit .* squeeze(max(max(ipres(:,:,1),[],1),[], ...
                                           2));
        end

        % get proper pressure & velocity regions
        maskp(:,:,tsave) = find_mask(ipres(:,:,tsave), pcrit, imx(tsave), ...
                                     imy(tsave));

        dipdy = avg1(maskp(:,:,tsave),2) .* ...
                bsxfun(@rdivide, -diff(ipres(:,:,tsave),1,2), diff(yvec));
        % d/dx ∫P ~ d/dy ∫P ~ 1e3
        %dipdx = avg1(maskp,1) .* bsxfun(@rdivide, -diff(ipres,1,1), diff(xvec));
        %dipresdx(tsave) = integrate(avg1(xvec), yvec, dipdx);
        dipresdy(tsave) = integrate(xvec, avg1(yvec), dipdy);
    end

    clear rho zwmat dipdy pres
    %iU = cumtrapz(yvec, U, 2); % crude streamfunction estimate
    %iV = cumtrapz(xvec, V, 1); % crude streamfunction estimate

    uarea = integrate(xvec, yvec, masku);
    parea = integrate(xvec, yvec, bsxfun(@times, maskp, slbot > 0));

    btrq = integrate(xvec, yvec, ...
                     bsxfun(@times, pbot .* maskp, slbot))./parea;

    f0u = integrate(xvec, avg1(yvec), ...
                    f0 .* bsxfun(@rdivide, diff(AM,1,2), diff(yvec)) ...
                    .* avg1(masku,2))./uarea;

    f0v = integrate(avg1(xvec), yvec, ...
                    f0 .* bsxfun(@rdivide, diff(AM,1,1), diff(xvec)) ...
                    .* avg1(masku,1))./uarea;

    byu = integrate(xvec, yvec, beta .* AM .* masku)./uarea;

    %%%%%%%%%% plots
    nsmooth = runs.eddy.turnover./mean(diff(runs.time));
    figure; maximize(); pause(0.2);
    insertAnnotation([runs.name '.bottom_torque']);
    hold all
    plot(smooth(byu, nsmooth));
    plot(smooth(btrq, nsmooth));
    %plot(abs(runs.angmom.sym_betatrq)./(pi*runs.eddy.Lfit.^2));
    plot(smooth(f0u, nsmooth)); %plot(f0v); %plot(dipresdy'./parea(1:end-1,:));
    legend('\beta yu', 'p_{bot}', 'f_0 u', 'f_0 v');
    linex(runs.traj.tind); liney(0);
    ylim([-0.5 1]*max(byu(:)));
    title([runs.name ' |  subtract\_mean = ' num2str(flags.subtract_mean) ...
          ' | subtract\_edge = ' num2str(flags.subtract_edge)]);
    beautify;

    export_fig('-painters', ['images/angmom-' runs.name '-2.png']);

    %%%%%%%%% Summarize
    save([runs.dir '/pbot.mat'], 'pbot', 'slbot', 'masku', 'maskp', ...
         'AM', 'ipres', 'xvec', 'yvec');

    bottom.uarea = uarea;
    bottom.parea = parea;
    bottom.f0u = f0u;
    bottom.byu = byu;
    bottom.dipresdy = dipresdy./parea(1:end-1)';
    bottom.btrq = btrq;
    bottom.pcrit = pcrit;
    bottom.amcrit = amcrit;
    bottom.time = runs.eddy.t(tind(1):dt:tind(2))*86400;
    bottom.maskstr = maskstr;
    bottom.flags = flags;

    bottom.comment = ['(pressure, angmom) = volume integrated ' ...
                      'pressure, angular momentum | btrq = slope ' ...
                      '* pressure | byu = beta .* angmom'];

    bottom.hash = githash([mfilename('fullpath') '.m']);

    runs.bottom = bottom;
    save([runs.dir '/bottom.mat'], 'bottom', '-v7.3');

    %keyboard;

    %keyboard;

    animation = 1;
    if animation
        %umask = AM .* masku;
        %pmask = pbot .* maskp;
        %var = avg1(runs.ubot(volumeu{1,2}:volumeu{1,3}, ...
        %                volumeu{2,2}:volumeu{2,3}, :),1);
        var = bsxfun(@times, pbot .* maskp, slbot);
        t0 = 20;
        tt = t0;
        hp = pcolor(xvec, yvec, double(var(:,:,tt)'));
        cbfreeze;center_colorbar; shading flat; hold all;
        plot(runs.eddy.mx, runs.eddy.my, 'k');
        hc = plot(runs.eddy.mx(tt), runs.eddy.my(tt), 'k*');
        plot(runs.eddy.mx, runs.eddy.my - runs.eddy.vor.dia/2, 'k');
        liney(runs.bathy.xsb);
        for tt=t0+1:4:size(var,3)
            set(hp,'CData', double(var(:,:,tt)'));
            set(hc,'XData', runs.eddy.mx(tt), ...
                   'YData', runs.eddy.my(tt));
            pause(0.5);
        end
    end
    toc(ticstart);

end

function [out] = integrate(xvec, yvec, in)
    out = squeeze(trapz(yvec, ...
                        trapz(xvec, double(in), 1), 2));
end

function [out] = find_mask(in,crit,imx,imy)

    out = logical(zeros(size(in)));
    for kk=1:size(in,3)
        masktemp = in(:,:,kk) > crit;

        % first find simply connected regions
        regions = bwconncomp(masktemp, 8);

        clear masktemp;

        for rr = 1:regions.NumObjects
            maskreg = logical(zeros(regions.ImageSize));
            maskreg(regions.PixelIdxList{rr}) = 1;

            % center location
            if maskreg(imx(kk), imy(kk)) == 1
                out(:,:,kk) = maskreg;
                break;
            end
        end
    end
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


    %%%%%%%%% Translation term
    %c = runs.eddy.cvx(tind(1):dt:tind(2)) .* 1000/86400; % convert to m/s
    % c = smooth(runs.eddy.mvx(tind(1):dt:tind(2)), 10) .* 1000/86400; % convert to m/s

    % % height anomaly for eddy is zeta
    % h = bsxfun(@minus, zeta, mean(zeta, 2));

    % iv = bsxfun(@times, bsxfun(@times, h, f), permute(c, [3 2 1]));
    % %iv2 = bsxfun(@times, bsxfun(@times, irho, f), permute(c, [3 1 2]));
    % %iv = runs.params.phys.f0 .* U;

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


%%%%%%%%% mask?
    % changing this mask threshold gives me larger pressures
    % pcrit = 0.1;
    % botmask = pbot < pcrit*min(pbot(:));
    % mask_rho = sshmask; %botmask; %irho < -1;
    % mpbot = mask_rho .* pbot .* slbot;
    % mpbotneg = mpbot .* (mpbot < 0);
    % miv = mask_rho .* iv;
    % if flags.calc_angmom
    %     miam = mask_rho .* iam2;
    % end

    % clear V P AM
    % %%%%%%%%% area-integrate - axes referenced to center
    % for tt=1:size(pbot,3)
    %     P(tt) = squeeze(trapz(yrmat(1,:,tt), ...
    %                           trapz(xrmat(:,1,tt), repnan(mpbot(:,:,tt),0), ...
    %                                 1), 2));
    %     Pneg(tt) = squeeze(trapz(yrmat(1,:,tt), ...
    %                           trapz(xrmat(:,1,tt), repnan(mpbotneg(:,:,tt),0), ...
    %                                 1), 2));
    %     if flags.calc_angmom
    %         AM(tt) = squeeze(trapz(yrmat(1,:,tt), ...
    %                                trapz(xrmat(:,1,tt), repnan(miam(:,:,tt),0), ...
    %                                      1), 2));
    %     end
    %     %V(tt) = squeeze(trapz(yrmat(1,:,tt), ...
    %     %                     trapz(xrmat(:,1,tt), repnan(miv(:,:,tt),0), ...
    %     %                           1), 2));
    % end

    % figure;
    % hold all
    % plot(P);
    % plot(Pneg);
    % plot(runs.angmom.sym_betatrq);
    % title(runs.name);

    % calculate bottom pressure (x,y,t)
    % note that in Flierl (1987) the 1/ρ0 is absorbed into the
    % pressure variable
    %pres = bsxfun(@plus, g./rho0 .* irho, g.*permute(zeta,[1 2 4
    %3]));
    % pres = -g/ρ0 ∫_{z}^{ζ} ρ dz
    % -1*dzmat to integrate from _surface_ to bottom cumulatively

    % figure out initial error.
    % u1 = avg1(dc_roms_read_data(runs.dir, 'u', 1, volumeu, [], ...
    %                             runs.rgrid, 'his', 'single'), 1);
    % rho1 = dc_roms_read_data(runs.dir, 'rho', 1, volumer, [], ...
    %                          runs.rgrid, 'his', 'single') + 1000;
    % rho1 = bsxfun(@minus, rho1, rback);
    % masked = bsxfun(@and, rho1 < rhothreshvor, ...
    %                 permute(sshmask(:,:,1), [1 2 4 3]));

    % % angular momentum correction
    % U1full = sum(u1 .* dzmat0(:,:,:,1), 3);
    % U1ed = sum(u1.*masked .* dzmat0(:,:,:,1), 3);
    % amfactor =  integrate(xvec, yvec, bymat(:,:,1) .* U1full) ...
    %     ./ integrate(xvec, yvec, bymat(:,:,1) .* U1ed);
    % % bottom pressure correction
    % pbfactor = integrate(xvec,yvec, sum(rho1 .* dzmat0, 3)) ./ ...
    %     integrate(xvec,yvec, sum(rho1 .* masked .* dzmat0, 3));
    % clear u1 rho1 masked U1full U1ed


        % if ~flags.use_thermal_wind
        %     if flags.use_davg
        %         % depth averaged velocities (m/s)
        %ubar = dc_roms_read_data(runs.dir, 'ubar', [], ...
        %                                  volumer, [], runs.rgrid, 'his', 'single');
        %         vbar = dc_roms_read_data(runs.dir, 'vbar', [tstart tend], ...
        %                                  volumer, [], runs.rgrid, 'his', 'single');

        %         % convert to depth integrated velocities (m^2/s)
        %        U = bsxfun(@times, H, ubar);
        %         V = bsxfun(@times, H, vbar);
        %     else
        %         % read depth dependent velocity fields and integrate
        %         u = avg1(dc_roms_read_data(runs.dir, 'u', [tstart tend], volumeu, [], ...
        %                                    runs.rgrid, 'his', 'single'), 1);
        %         if mom_budget
        %             v = avg1(dc_roms_read_data(runs.dir, 'v', [tstart tend], volumev, [], ...
        %                                        runs.rgrid, 'his', 'single'), 2);
        %         end

        %         if flags.use_masked
        %             disp('Using rho based eddy mask.');

        %             % mask out velocities
        %             tic;
        %             masked = bsxfun(@and, bsxfun(@minus,rho,rback) < rhothreshvor, ...
        %                             permute(sshmask(:,:,tstart:dt:tend), ...
        %                                     [1 2 4 3]));
        %             toc;
        %             %masked = rho > rhothreshssh;

        %             u = u .* masked;
        %             rho = rho .* masked;
        %             if mom_budget
        %                 v = v .* masked;
        %             end
        %         end

        %         % depth-integrate quantities
        %         tic;
        %         U(:,:,tsave) = squeeze(sum(bsxfun(@times, u, dzmat), 3));
        %         if mom_budget
        %             V = squeeze(sum(bsxfun(@times,    v, dzmat), 3));
        %             UV = squeeze(sum(bsxfun(@times, u.*v, dzmat), 3));
        %             U2 = squeeze(sum(bsxfun(@times, u.^2, dzmat), 3));
        %             V2 = squeeze(sum(bsxfun(@times, v.^2, dzmat), 3));
        %             P = squeeze(sum(bsxfun(@times, pres, dzmat), 3));
        %         end
        %         toc;

        %         % try depth integrated momentum budget
        %         if mom_budget
        %             % pressure gradients
        %             dpdx = integrate(avg1(xvec), yvec, ...
        %                              bsxfun(@rdivide, diff(P,1,1), diff(xvec')));
        %             dpdy = integrate(xvec, avg1(yvec), ...
        %                              bsxfun(@rdivide, diff(P,1,2), diff(yvec)));

        %             % coriolis terms
        %             fv = integrate(xvec, yvec, f .* V);
        %             fu = integrate(xvec, yvec, f .* U);
        %             f0u = integrate(xvec, yvec, f0 .* U);
        %             byu = integrate(xvec, yvec, bymat .* U);
        %             f0v = integrate(xvec, yvec, f0 .* V);
        %             byv = integrate(xvec, yvec, bymat .* V);

        %             % non-linear terms
        %             dv2dy = integrate(xvec, avg1(yvec), ...
        %                               bsxfun(@rdivide, diff(V2,1,2), diff(yvec)));
        %             duvdx = integrate(avg1(xvec), yvec, ...
        %                               bsxfun(@rdivide, diff(UV,1,1), diff(xvec')));
        %             % tendency term - THIS IS A BAD ESTIMATE
        %             %dvdt = squeeze(trapz(trapz(diff(V,1,3)./86400,1),2));
        %             % bottom torque
        %             btq = integrate(xvec, yvec, pbot .* slbot);

        %             total = duvdx + dv2dy + fu + dpdy + btq;
        %             figure; hold all;
        %             plot(-1*f0u./total);
        %             plot(-1*byu./total);
        %             plot(dpdy./total);
        %             plot(duvdx./total);
        %             plot(dv2dy./total);
        %             plot(btq./total);
        %             legend('-f_0u','\beta yu', 'dpdy','duvdx','dv2dy', ...
        %                    'btq');

        %             time = runs.eddy.t(tind);
        %             save([runs.dir '/mombudget.mat'], 'dpdx', 'dpdy', 'fu', ...
        %                  'fv', 'f0u', 'byu', 'dv2dy', 'duvdx', 'btq', 'total', ...
        %                  'time');
        %         end
        %     end
        % end

        % if flags.use_thermal_wind
        %     % estimate velocity field associated with rho
        %     sz = flip(size(zrmat));
        %     grd.xmat = repmat(xvec', [1 sz(2) sz(3)]);
        %     grd.ymat = repmat(yvec , [sz(1) 1 sz(3)]);
        %     grd.zmat = permute(zrmat, [3 2 1]);
        %     dRdx = diff_cgrid(grd, rho, 1);
        %     dRdy = diff_cgrid(grd, rho, 2);
        %     uzest = -g./rho0 .* dRdy / f0;
        %     vzest = g./rho0 .* dRdx / f0;
        %     % geostrophic velocity
        %     ugest = cumsum(bsxfun(@times,uzest, avg1(avg1(dzmat,2),3)), ...
        %                    3);
        %     vgest = cumsum(bsxfun(@times,vzest, avg1(avg1(dzmat,1),3)), ...
        %                    3);
        %     % gradient wind
        % end
        % % prsgrd32.h
        % if flags.use_prsgrd

        %     dR = nan(size(zwmat)); dZ = nan(size(zwmat));
        %     dR(:,:,2:end-1,:) = diff(rho,1,3);
        %     dZ(:,:,2:end-1,:) = diff(zrmat,1,3);

        %     dR(:,:,end,:) = dR(:,:,end-1,:);
        %     dZ(:,:,end,:) = dZ(:,:,end-1,:);

        %     dR(:,:,1,:) = dR(:,:,2,:);
        %     dZ(:,:,1,:) = dZ(:,:,2,:);

        %     N = runs.rgrid.N; tic;
        %     for kk=N+1:-1:2
        %         dZ(:,:,kk,:) = 2 * dZ(:,:,kk,:) .* dZ(:,:,kk-1,:) ...
        %             ./ (dZ(:,:,kk,:) + dZ(:,:,kk-1,:));
        %         cff = 2*dR(:,:,kk,:) .* dR(:,:,kk-1,:);
        %         cff(cff < 1e-10) = 0;
        %         dR(:,:,kk,:) = cff ./ (dR(:,:,kk,:) + dR(:,:,kk-1,:));
        %     end
        %     toc;

        %     tic;
        %     P = nan(size(rho));
        %     P(:,:,end,:) = 1000 * g./rho0 .* zwmat(:,:,end,:) + ...
        %         g/rho0 * (zwmat(:,:,end,:)-zrmat(:,:,end,:)) .* ...
        %         ( -1000 + rho(:,:,end,:) + ...
        %           1./(zrmat(:,:,end,:)-zrmat(:,:,end-1,:)) .* ...
        %           (rho(:,:,end,:)-rho(:,:,end-1,:)) .* ...
        %           (zwmat(:,:,end,:) - zrmat(:,:,end,:)));

        %     for kk=N-1:-1:1
        %         P(:,:,kk,:) = P(:,:,kk+1,:) + ...
        %             1/2*g/rho0 .* ( ...
        %                 (-2000 + rho(:,:,kk+1,:) + rho(:,:,kk,:)) ...
        %                 .* (zrmat(:,:,kk+1,:) - zrmat(:,:,kk,:)) ...
        %                 - 1/5 * ( (dR(:,:,kk+1,:) - dR(:,:,kk,:)) .* ...
        %                           (zrmat(:,:,kk+1,:) - zrmat(:,:,kk,:) - ...
        %                            1/12 * (dZ(:,:,kk+1,:) + dZ(:,:,kk,:))) ...
        %                           - (dZ(:,:,kk+1,:)-dZ(:,:,kk,:)) ...
        %                           .* ( rho(:,:,kk+1,:) - rho(:,:,kk,:) ...
        %                                - 1/12 * (dR(:,:,kk+1,:) + dR(:,:,kk,:)))));
        %     end
        %     P = bsxfun(@minus, P, P(1,:,:,:));
        %     pbot1(:,:,tsave) = squeeze(P(:,:,1,:));
        %     toc;
        % end
        % % check balance
        % if flags.use_thermal_wind
        %     ranom = bsxfun(@minus, rho, rback);
        %     ranom = bsxfun(@minus, ranom, ranom(1,:,:,:));

        %     % iranom = squeeze(sum(bsxfun(@times, ranom, dzmat0),3));
        %     % %diRdy = bsxfun(@rdivide, diff(iranom,1,2), diff(yvec));
        %     % %dizdy = -1/rho0 * diRdy;
        %     % %dRdy = bsxfun(@rdivide, diff(ranom, 1, 2), diff(yvec));
        %     % %dzdy = -1 * squeeze(sum(bsxfun(@times, dRdy, avg1(dzmat,2)), ...
        %     % %                   3))/rho0;
        %     % %dzetady = (bsxfun(@rdivide, diff(zeta(:,:,tsave),1,2), ...
        %     % %                       diff(yvec)));
        %     % %error = dzdy - dzetady;

        %     % pbc1 = g/rho0 * iranom;
        %     % pbt1 = g * zeta(:,:,tsave);

        %     % pbot1(:,:,tsave) = pbc1+pbt1;

        %     % THIS IS NOT HOW YOU DIFFERENTIATE ON A C-GRID
        %     drady = bsxfun(@rdivide, diff(ranom,1,2), diff(yvec));
        %     U(:,:,tsave) = -g./f0/rho0 .* squeeze(sum( bsxfun(@times, ...
        %                                                       cumsum( bsxfun(@times, drady, avg1(dzmat,2)), 3), ...
        %                                                       avg1(dzmat,2)), 3));
        %     %ubot = bsxfun(@rdivide, diff(pbot,1,2), diff(yvec))./f0;
        % end
