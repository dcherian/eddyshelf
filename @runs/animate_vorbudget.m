function [] = animate_vorbudget(runs,tind, plotflag)

    vorbudgetstart = tic;

    runs.vorbudget = [];

    if ~exist('tind','var')
        tind = 1;
    end
    %             if ~exist([runs.dir '/ocean_vor.nc'],'file')
    %                 dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
    %             end

    runs.video_init('vor');

    if ~exist('plotflag', 'var')
        plotflag = 1;
    end

    debug = 0;

    %%
    zmin = -1 * runs.bathy.hsb;min(runs.rgrid.z_r(:));
    zmax = max(runs.rgrid.z_r(:,end,end));
    zwnew = unique([linspace(zmin, -1*runs.bathy.hsb, 70) ...
                    linspace(-1*runs.bathy.hsb, zmax-0.01, 36)]');
    zrnew = avg1(zwnew);
    % for integrating quantities later
    zint = avg1(zrnew);

    % prepare grids for differentiation
    xvor = avg1(avg1(runs.rgrid.xr,1),2);
    yvor = avg1(avg1(runs.rgrid.yr,1),2);
    N = runs.rgrid.N;
    Nnew = length(zrnew);

    % setup grid
    [sx sy] = size(runs.rgrid.x_rho');
    gridu.xmat = repmat(runs.rgrid.x_u',[1 1 Nnew]);
    gridu.ymat = repmat(runs.rgrid.y_u',[1 1 Nnew]);
    gridu.zmat = permute(runs.rgrid.z_u,[3 2 1]);
    gridu.znew = repmat(permute(zrnew,[3 2 1]),[sx-1 sy 1]);
    %gridu.s = runs.rgrid.s_rho;
    %gridu.zw = runs.rgrid.z_w;
    %gridu.s_w = runs.rgrid.s_w;

    gridv.xmat = repmat(runs.rgrid.x_v',[1 1 Nnew]);
    gridv.ymat = repmat(runs.rgrid.y_v',[1 1 Nnew]);
    gridv.zmat = permute(runs.rgrid.z_v,[3 2 1]);
    gridv.znew = repmat(permute(zrnew,[3 2 1]),[sx sy-1 1]);
    %gridv.s = runs.rgrid.s_rho;
    %gridv.zw = runs.rgrid.z_w;
    %gridv.s_w = runs.rgrid.s_w;

    gridr.xmat = repmat(runs.rgrid.x_rho',[1 1 Nnew]);
    gridr.ymat = repmat(runs.rgrid.y_rho',[1 1 Nnew]);
    gridr.zmat = permute(runs.rgrid.z_r,[3 2 1]);
    gridr.znew = repmat(permute(zrnew,[3 2 1]),[sx sy 1]);
    %gridr.s = runs.rgrid.s_rho;
    %gridr.zw = runs.rgrid.z_r;
    %gridr.s_w = runs.rgrid.s_w;

    gridw.xmat = repmat(runs.rgrid.x_rho',[1 1 Nnew+1]);
    gridw.ymat = repmat(runs.rgrid.y_rho',[1 1 Nnew+1]);
    gridw.zmat = permute(runs.rgrid.z_w,[3 2 1]);
    gridw.znew = repmat(permute(zwnew,[3 2 1]),[sx sy 1]);
    %gridw.s = runs.rgrid.s_w;
    %gridw.zw = runs.rgrid.z_w;
    %gridw.s_w = runs.rgrid.s_w;

    gridrv.xmat = repmat(xvor,[1 1 Nnew]);
    gridrv.ymat = repmat(yvor,[1 1 Nnew]);
    gridrv.zmat = avg1(avg1(avg1(permute(runs.rgrid.z_r,[3 2 1]),1),2),3);
    gridrv.znew = repmat(permute(zrnew,[3 2 1]),[sx-1 sy-1 1]);
    gridrv.s = avg1(runs.rgrid.s_rho);
    gridrv.zw = avg1(avg1(avg1(permute(runs.rgrid.z_w,[3 2 1]),1),2),3);
    gridrv.s_w = avg1(runs.rgrid.s_w);

    % for depth integration - banas code
    h = runs.bathy.h(2:end-1,2:end-1);
    %        csr = runs.rgrid.Cs_r(2:end-1);
    %csw = runs.rgrid.Cs_w(2:end-1);

    % for depth-averaging
    hmax = max(abs(zint)); % max. depth of integration
    hmat = h .* (h <= hmax) + hmax .* (h > hmax);
    % add in sponge mask
    hmat = hmat .* fillnan(~runs.sponge(2:end-1, 2:end-1), 0);

    % for bottom friction I need to mask out the area that
    % doesn't touch the bottom
    hbfric = fillnan(hmat .* (hmat == runs.bathy.h(2:end-1,2:end-1)), ...
                     0);

    xavg = avg1(avg1(xvor,1),2)/1000; yavg = avg1(avg1(yvor,1),2)/1000;

    % AREA AVERAGING - for bottom friction terms
    dA = 1./runs.rgrid.pm(2:end-1,2:end-1)' .* 1./runs.rgrid.pn(2:end-1, ...
                                                      2:end-1)';
    dA = dA .* ~runs.sponge(2:end-1, 2:end-1);
    area = sum(dA(:));

    % VOLUME AVERAGING
    % 2D array - water column volume for each (x,y) - masked
    %dVxy = 1./runs.rgrid.pm(2:end-1,2:end-1)' .* 1./runs.rgrid.pn(2:end-1, 2:end-1)' ...
    %.* hmat;

    % 3D array - cell volumes for each (x,y,z)
    % nansum(dV(:))  ~= nansum(dVxy(:)) since,
    % i'm integrating to a level just
    % above the bottom.
    zmat = repmat(permute(zrnew, [3 2 1]), [size(hmat,1) ...
                        size(hmat,2)]);
    zmat(bsxfun(@lt, zmat, -1 * hmat)) = NaN;
    dV = bsxfun(@times, ...
                bsxfun(@times, dA, diff(zmat, 1, 3)), ...
                ~isnan(hmat));
    vol = nansum(dV(:));
    %disp(['error in volumes = ' num2str((vol - nansum(dVxy(:)))./vol ...
    %                                    * 100) ' percent']);

    % time range and file reading parameters
    slab = 12;
    stride = 1;

    timehis = dc_roms_read_data(runs.dir, 'ocean_time', [], {}, ...
                                [], runs.rgrid, 'his');
    trange = tind:stride:length(timehis);
    disp(['starting from t instant = ' num2str(trange(1))]);

    % save vorticity budget for whole domain
    runs.vorbudget.hadv = nan([length(trange) 1]);
    runs.vorbudget.vadv = runs.vorbudget.hadv;
    runs.vorbudget.tilt = runs.vorbudget.hadv;
    runs.vorbudget.str  = runs.vorbudget.hadv;
    runs.vorbudget.beta = runs.vorbudget.hadv;
    runs.vorbudget.bfric = runs.vorbudget.hadv;
    %runs.vorbudget.sol = runs.vorbudget.hadv;
    %runs.vorbudget.budget = runs.vorbudget.hadv;

    % vorticity budget for shelf water
    runs.vorbudget.shelf.hadv = runs.vorbudget.hadv;
    runs.vorbudget.shelf.vadv = runs.vorbudget.hadv;
    runs.vorbudget.shelf.str = runs.vorbudget.hadv;
    runs.vorbudget.shelf.tilt = runs.vorbudget.hadv;
    runs.vorbudget.shelf.beta = runs.vorbudget.hadv;
    runs.vorbudget.shelf.bfric = runs.vorbudget.hadv;

    runs.vorbudget.comment = ['hadv + vadv + beta = str + tilt  ' ...
                        '+ bfric'];
    %runs.vorbudget.conthis = runs.vorbudget.hadv;
    %%

    for kk=1:slab:length(trange)
        tt = trange(kk);
        disp(['kk = ' num2str(kk/slab) '/' num2str(length(trange)/slab) ...
              ' | tt = ' num2str(tt/2) ' days | plotflag = ' ...
              num2str(plotflag) ' | run = ' runs.name]);
        %zeta = runs.zeta(2:end-1,2:end-1,tt);
        % read data
        %fname = [runs.dir '/ocean_his.nc.new2'];
        %fname = runs.out_file;
        %w  = double(ncread(fname,'w',[1 1 1 tt],[Inf Inf Inf 1]));
        %zeta = double(ncread(fname,'zeta',[1 1 tt],[Inf Inf 1]));

        if stride ~= 1
            error('stride does not work');
        end

        % read in history file data
        tindices = [tt tt+stride*slab-1]
        if tt+stride*slab-1 > trange(end)
            tindices(end) = trange(end)
        end
        uh = dc_roms_read_data(runs.dir,'u',tindices,{},[],runs.rgrid, ...
                               'his', 'single');
        vh = dc_roms_read_data(runs.dir,'v',tindices,{},[],runs.rgrid, ...
                               'his', 'single');
        wh = dc_roms_read_data(runs.dir,'w',tindices,{},[],runs.rgrid, ...
                               'his', 'single');
        csdye = dc_roms_read_data(runs.dir, runs.csdname, tindices, ...
                                  {}, [], runs.rgrid, 'his', 'single');
        %zeta = dc_roms_read_data(runs.dir, 'zeta', tt, {}, [], ...
        %                         runs.rgrid, 'his', 'single');

        %rhoh = dc_roms_read_data(runs.dir,'rho',tt,{},[],runs.rgrid, ...
        %                         'his');

        %ubar = dc_roms_read_data(runs.dir, 'ubar', tt, {}, [], ...
        %                         runs.rgrid, 'his');
        %vbar = dc_roms_read_data(runs.dir, 'vbar', tt, {}, [], ...
        %                         runs.rgrid, 'his');

        % interpolate to znew depths
        disp('interpolating variables');
        u = single(interpolate(uh, gridu.zmat, zrnew));
        v = single(interpolate(vh, gridv.zmat, zrnew));
        w = single(interpolate(wh, gridw.zmat, zwnew));
        csd = single(interpolate(csdye, gridr.zmat, zrnew));
        % rho = interpolate(rhoh, gridr.zmat, zrnew);

        ux = bsxfun(@rdivide, diff(u,1,1), diff(gridu.xmat,1,1));
        uy = bsxfun(@rdivide, diff(u,1,2), diff(gridu.ymat,1,2));
        uz = bsxfun(@rdivide, diff(u,1,3), diff(gridu.znew,1,3));

        vx = bsxfun(@rdivide, diff(v,1,1), diff(gridv.xmat,1,1));
        vy = bsxfun(@rdivide, diff(v,1,2), diff(gridv.ymat,1,2));
        vz = bsxfun(@rdivide, diff(v,1,3), diff(gridv.znew,1,3));

        wx = bsxfun(@rdivide, diff(w,1,1), diff(gridw.xmat,1,1));
        wy = bsxfun(@rdivide, diff(w,1,2), diff(gridw.ymat,1,2));
        wz = bsxfun(@rdivide, diff(w,1,3), diff(gridw.znew,1,3));

        %rx = diff(rho,1,1)./diff(gridr.xmat,1,1);
        %ry = diff(rho,1,2)./diff(gridr.ymat,1,2);
        %rz = diff(rho,1,3)./diff(gridr.znew,1,3);

        %cont = ux(:, 2:end-1, :) + vy(2:end-1, :, :) + ...
        %       wz(2:end-1, 2:end-1, :);

        % check cont
        %ix = 150; iy = 164;
        %ix = 240; iy = 164
        %figure; hold all;
        %hold all;
        %plot(squeeze(ux(ix, iy+1,:)) + squeeze(vy(ix+1, iy,:)), zrnew)
        %plot(-1*squeeze(wz(ix+1, iy+1,:)), zrnew)
        %plot(squeeze(cont(ix, iy, :)), zrnew);
        %legend('ux + vy', 'wz', 'ux + vy + wz');

        % tendency term code - not really needed since it is probably a
        % bad estimate when using daily snapshots .
        %             if debug
        %                 u1 = interpolate(u1, gridu.zmat, znew);
        %                 v1 = interpolate(v1, gridv.zmat, znew);
        %                 v1x = diff(v1,1,1)./diff(gridv.xmat,1,1);
        %                 u1y = diff(u1,1,2)./diff(gridu.ymat,1,2);
        %                 rv1 = v1x-u1y;
        %             end
        rv = vx-uy;
        rvx = bsxfun(@rdivide, diff(rv,1,1), diff(gridrv.xmat,1,1));
        rvy = bsxfun(@rdivide, diff(rv,1,2), diff(gridrv.ymat,1,2));
        rvz = bsxfun(@rdivide, diff(rv,1,3), diff(gridrv.znew,1,3));

        rvavg = avg1(avg1(avg1(rv, 1), 2), 3);

        if debug
            u1h = double(ncread(runs.dir,'u',[1 1 1 tt+1],[Inf Inf Inf 1]));
            v1h = double(ncread(runs.dir,'v',[1 1 1 tt+1],[Inf Inf ...
                                Inf 1]));
            t1 = double(ncread(runs.dir, 'ocean_time'));

            u1 = interpolate(u1h, gridu.zmat, zrnew);
            v1 = interpolate(v1h, gridv.zmat, zrnew);

            u1y = diff(u1,1,2)./diff(gridu.ymat,1,2);
            v1x = diff(v1,1,1)./diff(gridv.xmat,1,1);

            rv1 = v1x-u1y;

            % calculate term and average to agree with 'budget' size
            drvdt = avg1(avg1(avg1( ...
                (rv1-rv)./(t1(tt+1)-t1(tt)), 1), 2), 3);

            DRVDT = trapz(zint, repnan(drvdt, 0), 3)./hmat;
        end

        str = avg1(avg1(avg1(bsxfun(@plus, rv, ...
                                    avg1(avg1(runs.rgrid.f',1),2)),1) ...
                        ,2) .* -1 .* (ux(:,2:end-1,:,:) + vy(2:end-1,:,:,:)),3);
        %wz(2:end-1,2:end-1,:), 3);%

        tilt = -1 * avg1(avg1( avg1(wx(:,:,2:end-1,:),2) .* avg1(vz,1) + ...
                               avg1(wy(:,:,2:end-1,:),1) .* avg1(uz,2) ,1),2);
        beta = avg1(avg1(runs.params.phys.beta * v(2:end-1,:,:,:),2),3);
        hadv = avg1( avg1(u(:,2:end-1,:,:),1) .* avg1(rvx,2) + ...
                     avg1(v(2:end-1,:,:,:),2) .* avg1(rvy,1),3);
        vadv = avg1(avg1( avg1(avg1(w(:,:,2:end-1,:),1),2) .* rvz ...
                          ,1),2);

        budget = str + tilt - hadv - vadv - beta;

        % shelf water budget
        % shelf water mask defined with csdye + I remove sponge
        % region based on filtering already done in hmat
        shelfmask = bsxfun(@times, (avg1(csd(2:end-1, 2:end-1, :, :),3) < ...
                                    runs.bathy.xsb), ~isnan(hmat));
        %shelfmaskrv = bsxfun(@times, avg1(avg1(csd,1),2) < ...
        %                            runs.bathy.xsb, ~isnan(hmat));
        %  sol = -runs.params.phys.g/runs.params.phys.rho0 .* ...
        %          ( avg1(rx,2) .* avg1(zy,1) - avg1(ry,1) .* avg1(zx,2));

        % need bottom vorticity for bfric calculation
        rvbot = nan(size(squeeze(rvavg(:,:,1,:))));
        shelfmaskbot = rvbot;

        if runs.params.misc.rdrg ~= 0
            tic;
            disp('calculating bottom vorticity');
            if kk == 1
                % valid cells never change, so save mask
                % (botmask) that when multiplied with field
                % gives me the bottom values.
                botmask  = nan(size(rvavg(:,:,:,1)));
                for kkk = 1:size(rvbot, 3)
                    for iii = 1:size(rvbot,1)
                        for jjj = 1:size(rvbot,2)
                            % locate first 0  since z=1 is bottom
                            zind = find(isnan(squeeze(rvavg(iii,jjj,:,kkk))) ...
                                        == 0, 1, 'first');
                            if ~isempty(zind)
                                botmask(iii,jjj,zind) = 1;
                            end
                        end
                    end
                end
            end
            rvbot = squeeze(nansum(bsxfun(@times, rvavg, botmask),3));
            shelfmaskbot = squeeze(nansum(bsxfun(@times, shelfmask, ...
                                                 botmask),3));
            ubot = squeeze(nansum(bsxfun(@times, avg1(avg1(u(:, ...
                                                             2: ...
                                                             end-1,:,:),3),1), botmask), 3)) .* shelfmaskbot;
            ;
            vbot = squeeze(nansum(bsxfun(@times, avg1(avg1(v(2: ...
                                                             end-1,:,:,:),3),2), botmask), 3)) .* shelfmaskbot;;

            toc;
        end

        % depth INTEGRATED QUANTITIES
        RV   = avg1(avg1(trapz(zrnew, repnan(rv,0), 3),1), 2);
        %RVSHELF = trapz(zrnew, repnan(rvavg,0) .* ...
        %                shelfmask, 3);

        % depth - AVERAGED quantities for plotting
        STR  = squeeze(bsxfun(@rdivide, trapz(zint, repnan(str,0),  3), hmat));
        TILT = squeeze(bsxfun(@rdivide, trapz(zint, repnan(tilt,0), 3), hmat));
        BETA = squeeze(bsxfun(@rdivide, trapz(zint, repnan(beta,0), 3), hmat));
        HADV = squeeze(bsxfun(@rdivide, trapz(zint, repnan(hadv,0), 3), hmat));
        VADV = squeeze(bsxfun(@rdivide, trapz(zint, repnan(vadv,0), 3), hmat));
        ADV = HADV + VADV;

        % FRICTION only when integrating to bottom surface
        BFRIC = bsxfun(@times, bsxfun(@rdivide, -runs.params.misc.rdrg .* rvbot, ...
                                      hmat), hmat == h);
        BFRICSHELF = BFRIC .* shelfmaskbot;
        bfric = bsxfun(@times, -runs.params.misc.rdrg .* rvbot, ...
                       hmat == h);
        bfricshelf = bfric .* shelfmaskbot;

        % BUDGET = TEND = d(RV)/dt
        %BUD = STR + BFRIC + TILT - BETA - ADV;

        % ubar, vbar calculated for depth averaged interval
        % only
        % ubar = bsxfun(@rdivide, trapz(zrnew, repnan(avg1(u(:,2:end-1,:,:),1),0), 3) ...
        %            ,hmat);
        % vbar = bsxfun(@rdivide, trapz(zrnew, repnan(avg1(v(2:end-1,:,:,:),2),0), ...
        %              3), hmat);
        if debug
            BUD = BUD - DRVDT;
            imagesc(BUD');
        end
        %BUD = trapz(zint, repnan( str+tilt - beta - hadv
        %-vadv,0), 3);

        % volume of shelfwater - shelfmask has sponge taken out
        shelfvol = bsxfun(@times, shelfmask, dV);
        shelfvol = squeeze(nansum(nansum(nansum(shelfvol, 1), ...
                                         2), 3));

        % area of shelfwater in contact with bottom
        shelfarea = bsxfun(@times, shelfmaskbot, dA);
        shelfarea = squeeze(nansum(nansum(shelfarea, 1), 2));

        % reshape for volume averaging
        sz4d = size(rvavg);
        if length(sz4d) == 3
            sz4d(4) = 1;
        end
        sz2d = [sz4d(1)*sz4d(2)*sz4d(3) sz4d(4)];

        % calculate vorticity eqn terms - with shelfmask -
        % volume averaged
        indices = [tindices(1):tindices(end)] - trange(1) + 1;
        runs.vorbudget.shelf.vol(indices) = shelfvol;
        %            runs.vorbudget.shelf.area(indices) = shelfarea;
        runs.vorbudget.shelf.rv(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, rvavg .* shelfmask, dV),1), 2), 3)) ...
            ./ shelfvol;
        runs.vorbudget.shelf.str(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, str .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
        runs.vorbudget.shelf.tilt(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, tilt .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
        runs.vorbudget.shelf.hadv(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, hadv .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
        runs.vorbudget.shelf.vadv(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, vadv .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
        runs.vorbudget.shelf.beta(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, beta .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
        runs.vorbudget.shelf.bfric(indices) = squeeze(nansum(nansum( ...
            bsxfun(@times, bfricshelf, dA), 1), 2)) ./ ...
            shelfvol;

        bfricold = squeeze(nansum(nansum( ...
            bsxfun(@times, BFRICSHELF, dA), 1), 2)) ./ shelfarea;


        % save volume averaged quantities for whole domain
        runs.vorbudget.rv(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, rvavg, dV),1), 2), 3)) ...
            ./ vol;
        runs.vorbudget.str(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, str, dV),1), 2), 3)) ./ vol;
        runs.vorbudget.tilt(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, tilt, dV),1), 2), 3)) ./ vol;
        runs.vorbudget.hadv(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, hadv, dV),1), 2), 3)) ./ vol;
        runs.vorbudget.vadv(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, vadv, dV),1), 2), 3)) ./ vol;
        runs.vorbudget.beta(indices) = squeeze(nansum(nansum(nansum( ...
            bsxfun(@times, beta, dV),1), 2), 3)) ./ vol;
        runs.vorbudget.bfric(indices) = squeeze(nansum(nansum( ...
            bsxfun(@times, bfric, dA), 1), 2)) ./ vol;

        if plotflag
            limc = [-1 1] * nanmax(abs(ADV(:)));
            limy = [0 150];
            limx = [xvor(find(~runs.sponge(:,1) == 1, 1, 'first'),1) ...
                    xvor(find(~runs.sponge(:,1) == 1, 1, 'last'),1)]/1000;
            titlestr = 'Depth integrated rvor';
            % plot
            if kk == 1
                figure; maximize();
                ax(1) = subplot(2,4,[1:2]);
                hvor = pcolor(xavg, yavg, RV); hold on; shading flat;
                axis image;
                ht = runs.set_title('Depth int rvor', ceil(tt/2));
                he(1) = runs.plot_eddy_contour('contour',ceil(tt/2));
                hbathy = runs.plot_bathy('contour','k');
                shading flat
                caxis([-1 1] * nanmax(abs(RV(:))));
                colorbar;
                ylim(limy); xlim(limx);

                ax(2) = subplot(2,4,3);

                if runs.params.misc.rdrg == 0
                    hbet = pcolor(xavg,yavg,-BETA);
                    title('- \beta V');
                else
                    hbet = pcolor(xavg, yavg, BFRIC);
                    title('Bottom Friction');
                end
                colorbar; shading flat;
                he(2) = runs.plot_eddy_contour('contour', ceil(tt/2));
                hbathy = runs.plot_bathy('contour','k');
                caxis(limc); %caxis([-1 1] * nanmax(abs(BETA(:))));

                ylim(limy); xlim(limx);

                ax(3) = subplot(2,4,4); cla
                xran = 1:6:size(xavg,1); yran = 1:4:size(yavg,2);
                hquiv = quiver(xavg(xran,yran),yavg(xran,yran), ...
                               ubar(xran,yran), vbar(xran,yran),1.5);
                title('(ubar,vbar)');
                he(3) = runs.plot_eddy_contour('contour', ceil(tt/2));
                hbathy = runs.plot_bathy('contour','k');
                ylim(limy); xlim(limx);

                %                 ax(4) = subplot(2,4,5);
                %                 htend = pcolor(xavg,yavg,TEND); colorbar; shading flat;
                %                 he(4) = runs.plot_eddy_contour('contour', ceil(tt/2));
                %                 hbathy = runs.plot_bathy('contour','k');
                %                 caxis([-1 1] * max(abs(TEND(:))));
                %                 title('d\xi/dt');

                ax(5) = subplot(2,4,7);
                hgadv = pcolor(xavg,yavg,-ADV); colorbar; shading flat;
                he(5) = runs.plot_eddy_contour('contour', ceil(tt/2));
                hbathy = runs.plot_bathy('contour','k');
                caxis(limc); %caxis([-1 1] * max(abs(ADV(:))));
                title('-Advection');

                ax(6) = subplot(2,4,8);
                htilt = pcolor(xavg,yavg,TILT); colorbar; shading flat;
                he(6) = runs.plot_eddy_contour('contour', ceil(tt/2));
                hbathy = runs.plot_bathy('contour','k');
                caxis(limc/10); %caxis([-1 1] * max(abs(TILT(:))));
                title('Tilting');

                ax(7) = subplot(2,4,[5 6]);
                hstr = pcolor(xavg,yavg,STR); colorbar; hold on; shading flat;
                %hquiv = quiverclr(xavg(xran,yran),yavg(xran,yran), ...
                %    ubar(xran,yran),vbar(xran,yran),0.3,STR(xran,yran), ...
                %    [-1 1]*1e-11);
                %set(gca,'color',[0 0 0]);
                he(7) = runs.plot_eddy_contour('contour',ceil(tt/2));
                hbathy = runs.plot_bathy('contour','k');
                caxis(limc); %caxis([-1 1] * max(abs(STR(:))));
                title('Stretching = (f+\xi)w_z')
                spaceplots(0.06*ones([1 4]),0.05*ones([1 2]))
                linkaxes(ax,'xy');
                runs.video_update();
                pause();
            else
                set(hvor ,'cdata',RV);
                set(hgadv ,'cdata',-ADV);
                if runs.params.misc.rdrg == 0
                    set(hbet ,'cdata',-BETA);
                else
                    set(hbet, 'cdata', BFRIC);
                end
                set(hstr ,'cdata',STR);
                set(htilt,'cdata',TILT);
                %set(htend,'cdata',TEND);
                try
                    set(hquiv,'udata',ubar(xran,yran),'vdata',vbar(xran,yran));
                catch ME
                end

                runs.update_eddy_contour(he,ceil(tt/2));
                runs.update_title(ht,titlestr,ceil(tt/2));
                runs.video_update();
                pause(0.01);
            end
        end
    end

    if plotflag
        runs.video_write();
    end
    runs.vorbudget.time = timehis(trange);

    figure;
    plot(runs.vorbudget.time,-runs.vorbudget.hadv,'r'); hold on
    plot(runs.vorbudget.time,-runs.vorbudget.vadv,'g');
    plot(runs.vorbudget.time,runs.vorbudget.tilt,'b');
    plot(runs.vorbudget.time,runs.vorbudget.str,'c');
    %plot(runs.vorbudget.time,runs.vorbudget.sol,'m');
    plot(runs.vorbudget.time,-runs.vorbudget.beta,'y');
    %       plot(runs.vorbudget.time,runs.vorbudget.budget,'k');
    title('signs so that all terms are on RHS and tendency is LHS');
    legend('hadv','vadv','tilt','str','beta');

    vorbudget = runs.vorbudget;
    vorbudget.hash = githash;
    save([runs.dir '/vorbudget.mat'],'vorbudget');

    toc(vorbudgetstart);
end
