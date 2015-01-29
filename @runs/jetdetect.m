% study along-shore jet
function [] = jetdetect(runs)
    tictotal = tic;
    % rossby radius
    rr = runs.rrshelf;
    % number of rossby radii east of eddy to plot section
    nrr = 8;

    debug = 1;

    t0 = runs.tscaleind - 30;
    ix0 = vecfind(runs.rgrid.x_u(1,:),runs.eddy.vor.cx(t0:end));

    ilo = runs.bathy.isb - 10; % low index
    ihi = runs.bathy.isl; % high index

    % bottom dye concentration
    eddye = dc_roms_read_data(runs.dir, runs.eddname, [t0 Inf], ...
                              {runs.bathy.axis ilo ihi; ...
                        'z' runs.rgrid.N runs.rgrid.N}, [], runs.rgrid, 'his', 'single');
    % along-shore velocity on s=0
    asbot = dc_roms_read_data(runs.dir, 'u', [t0 Inf], ...
                              {runs.bathy.axis ilo ihi; ...
                        'z' 1 1}, [], runs.rgrid, 'his', 'single');

    % allocate variables
    runs.jet.xnose = nan(size(runs.time));
    runs.jet.ixnose = nan(size(runs.time));
    runs.jet.vscale = nan(size(runs.time));
    runs.jet.yscale = nan(size(runs.time));
    runs.jet.zscale = nan(size(runs.time));
    runs.jet.h = nan(size(runs.time));
    runs.jet.bc = nan(size(runs.time));
    run.jet.uprof = cell(size(runs.time));
    runs.jet.width = nan(size(runs.time));

    %% diagnostics
    % let's find location of nose
    thresh = 0.75
    sz = size(eddye);
    if runs.bathy.axis == 'y'
        xd = runs.rgrid.xr(:,ilo:ihi);
        yd = runs.rgrid.yr(:,ilo:ihi);
        edge = runs.eddy.we;
        bound = runs.eddy.ne;
        xdvec = xd(:,1);
    else
        xd = runs.rgrid.yr(ilo:ihi,:);
        edge = runs.eddy.se;
        xdvec = xd(1,:)';
    end
    xmask = reshape(xd, [sz(1)*sz(2) 1]);

    if runs.bathy.axis == 'y'
        % jet is east of eddy
        % find locations to the east of the eddy where dye > thresh
        masked = reshape((eddye .*  ...
                          bsxfun(@gt, xd, permute(edge(t0:end), [3 1 2])) .* ...
                          bsxfun(@lt, yd, permute(bound(t0:end), [3 1 2]))) ...
                          > thresh, [sz(1)*sz(2) sz(3)]);

    else
        % jet is south of eddy
        masked = reshape((eddye .*  bsxfun( ...
            @lt, xd, permute(edge(t0:end), [3 1 2])) ...
                          > thresh), [sz(1)*sz(2) sz(3)]);
    end
    [dmax, idmax] = max(bsxfun(@times, masked, xmask), [], 1);
    [ix, iy] = ind2sub([sz(1) sz(2)], idmax);
    runs.jet.xnose(t0:end) = fillnan(dmax,min(xmask(:)));
    runs.jet.ixnose(t0:end) = ix;
    runs.jet.thresh = thresh;

    % width at nose
    index = vecfind(xdvec, runs.jet.xnose(t0:end));
    index(runs.jet.xnose(t0:end) == 0) = NaN;
    tstart = find(~isnan(index) == 1, 1, 'first'); % W.R.T
                                                   % t0!!!!!

    if debug

        if isempty(runs.usurf), runs.read_velsurf; end
        cmap = flipud(cbrewer('div','RdBu',32));
        maskfull = reshape(masked, sz);
        svel = runs.usurf(:,ilo:ihi,t0:end);
        us = avg1(svel,1);
        umask = maskfull(2:end-1,:,:) .* us;
        if isempty(runs.vorsurf), runs.calc_vorsurf; end
        vor = avg1(runs.vorsurf(:,ilo:ihi,t0:end));
        vormask = maskfull(2:end-1,:,:) .* vor;
        dyemask = maskfull(2:end-1,:,:) .* eddye(2:end-1,:,:);

        tt = 52;
        figure;
        subplot(3,1,1); pcolorcen(umask(:,:,tt)'); center_colorbar; ...
            colormap(cmap); title('u'); axis image; freezeColors; cbfreeze;
        subplot(3,1,2); pcolorcen(vormask(:,:,tt)'); center_colorbar; ...
            colormap(cmap); title('vor'); axis image; freezeColors; ...
            cbfreeze;
        subplot(3,1,3); pcolorcen(dyemask(:,:,tt)'); caxis([0 1]); ...
            colorbar; colormap((cbrewer('seq','Reds', 32))); ...
            title('dye'); axis image; freezeColors
        %figure; pcolorcen((vormask(:,:,tt)<0)' .* dyemask(:,:,tt)'); caxis([0 1]); ...
        %    colorbar; colormap((cbrewer('seq','Reds', 32))); title('dye');

    end

    % read in data
    if runs.bathy.axis == 'y'
        [uprof,~,yu,zu,~] = dc_roms_read_data(runs.dir, 'u', [t0+tstart Inf], ...
                                              {'x' min(index)-1 max(index)-1; ...
                            'y' ilo ihi}, ...
                                              [], runs.rgrid, ...
                                              'his', 'single');
        dprof = dc_roms_read_data(runs.dir, runs.eddname, [t0+tstart Inf], ...
                                  {'x' min(index) max(index); ...
                            'y' ilo ihi}, ...
                                  [], runs.rgrid, 'his', 'single');
    else
        [uprof,yu,~,zu,~] = dc_roms_read_data(runs.dir, 'v', [t0+tstart Inf], ...
                                              {'y' min(index)-1 max(index)-1; ...
                            'x' ilo ihi}, ...
                                              [], runs.rgrid, ...
                                              [], 'single');
        dprof = dc_roms_read_data(runs.dir, runs.eddname, [t0+tstart Inf], ...
                                  {'y' min(index) max(index); ...
                            'x' ilo ihi}, ...
                                  [], runs.rgrid, [], 'single');
    end

    % 1 : take vertical profile of along-shore vel at index
    % 2 : find level of maximum velocity = velocity scale
    % 3 : then take cross shore section of velocity at that level
    %     (interpolated) and figure out scale.
    if runs.bathy.axis == 'y'
        yu = squeeze(yu(1,:,:));
        zu = squeeze(zu(1,:,:));
        h = runs.bathy.h(1, ilo:ihi);
    else
        yu = squeeze(yu(:,1,:));
        zu = squeeze(zu(:,1,:));
        h = runs.bathy.h(ilo:ihi,1);
    end
    ixmin = min(index); % needed for indexing

    % loop in time
    for ii=1:size(uprof,4)
        if isnan(index(tstart+ii-1)), continue; end
        % get y-z cross-section

        tind = t0 + tstart + ii - 1;
        if runs.bathy.axis == 'y'
            uvel = squeeze(uprof(index(tstart+ii-1)-ixmin+1,:,:, ...
                                 ii));
            dye  = squeeze(dprof(index(tstart+ii-1)-ixmin+1,:,:, ...
                                 ii));
        else
            uvel = squeeze(uprof(:,index(tstart+ii-1)-ixmin+1,:, ...
                                 ii));
            dye  = squeeze(dprof(:,index(tstart+ii-1)-ixmin+1,:, ...
                                 ii));
        end
        % find max. velocity
        [runs.jet.vscale(tind), ivmax] = max(uvel(:) .* ...
                                             (dye(:) > thresh));
        [iy,iz] = ind2sub(size(uvel), ivmax);
        % location of max. NOSE velocity in vertical
        runs.jet.zscale(tind) = zu(iy,iz);
        % location of max. NOSE velocity in cross-shore co-ordinate
        runs.jet.yscale(tind) = yu(iy,iz);

        % depth of water at location of max NOSE velocity
        runs.jet.h(tind) = h(iy);

        % baroclinicty of vertical profile at location of max
        % NOSE velocity
        runs.jet.bc(tind) = baroclinicity(zu(iy,:), uvel(iy,:));

        % width of jet at NOSE
        % first interpolate to get velocity at constant
        % z-level. this level is the location of
        % max. along-shore velocity i.e., jet.zscale(tind)
        ynew = yu(:,1);
        znew = ones(size(ynew)) .* runs.jet.zscale(tind);
        F = scatteredInterpolant(yu(:), zu(:), double(uvel(:)));
        unew = F(ynew, znew);
        % calculate auto-covariance, find first zero crossing
        % and multiply by 4 to get width
        ucov = xcov(unew);
        % symmetric, so discard first half
        ucov = ucov(length(ynew):end);
        iu = find(ucov < 0, 1, 'first');
        iu = iu-1;
        dy = min(1./runs.rgrid.pn(:));
        runs.jet.uprof{tind} = unew;
        runs.jet.width(tind) = 4 * dy * iu;
    end

    jet = runs.jet;
    jet.hash  = githash([mfilename('fullpath') '.m']);
    save([runs.dir '/jet.mat'], 'jet');

    if debug
        %% animation
        if isempty(runs.usurf), runs.read_velsurf; end
        svel = runs.usurf(:,ilo:ihi,t0:end);
        figure;
        ii=60;

        cmap = flipud(cbrewer('div','RdBu',32));

        ax(1) = subplot(311);
        hsv = pcolorcen(runs.rgrid.x_u(1,:)/1000, ...
                        runs.rgrid.y_u(ilo:ihi,1)/1000, ...
                        svel(:,:,ii)');
        hold on
        he1 = runs.plot_eddy_contour('contour',t0+ii-1);
        colorbar; caxis([-0.1 0.1]); colormap(cmap);%axis image
        hl1 = linex(runs.jet.xnose(t0+ii-1)/1000, '');
        title('along-shore surface velocity');

        ax(2) = subplot(3,1,2);
        hbv = pcolorcen(runs.rgrid.x_u(1,:)/1000, ...
                        runs.rgrid.y_u(ilo:ihi,1)/1000, ...
                        asbot(:,:,ii)');
        he2 = runs.plot_eddy_contour('contour',t0+ii-1);
        colorbar; caxis([-0.1 0.1]); colormap(cmap); %axis image
        hl2 = linex(runs.jet.xnose(t0+ii-1)/1000, '');
        title('Along-shore vel on s = 0');

        ax(3) = subplot(313);
        [hd] = pcolorcen(runs.rgrid.xr(:,1)/1000, ...
                         runs.rgrid.yr(1,ilo:ihi)/1000, ...
                         eddye(:,:,ii)');
        he3 = runs.plot_eddy_contour('contour',t0+ii-1);
        colorbar; caxis([0 1]); colormap(cmap); %axis image
        hl3 = linex(runs.jet.xnose(t0+ii-1)/1000, '');

        linkaxes(ax, 'xy');

        for ii=ii+1:size(eddye,3)
            set(hsv, 'cdata', double(svel(:,:,ii)'));
            set(hbv, 'cdata', double(asbot(:,:,ii)'));
            set(hd, 'cdata', double(eddye(:,:,ii)'));
            runs.update_eddy_contour(he1, t0+ii-1);
            runs.update_eddy_contour(he2, t0+ii-1);
            runs.update_eddy_contour(he3, t0+ii-1);
            set(hl1, 'xdata', [1 1]*runs.jet.xnose(t0+ii-1)/1000);
            set(hl2, 'xdata', [1 1]*runs.jet.xnose(t0+ii-1)/1000);
            set(hl3, 'xdata', [1 1]*runs.jet.xnose(t0+ii-1)/1000);
            title(['Dye on s=0 | day no = ' num2str(t0+ii) ', ii=', num2str(ii)])
            pause(0.05);
        end

        %% older animation showing cross-shore sections of
        %% along-shore velocity
        % first section moves with eddy
        %xind = ix0(1) + [nan 10 60] *ceil(rr/runs.rgrid.dx);

        %tt = 1;
        %xind(1) = ix0(tt) + nrr * ceil(rr/runs.rgrid.dx);
        %subplot(2,3,[1 2 3])
        %hzeta = runs.plot_zeta('contourf',t0+tt-1);
        %hlines = linex(xind*runs.rgrid.dx/1000);
        %colorbar

        %subplot(234)
        %[~,huas1] = contourf(yz/1000,squeeze(zas(xind(1),:,:)), ...
        %                     squeeze(uas(xind(1),:,:,tt)));
        %colorbar; caxis([-1 1]*0.1); ylim([-300 0]);

        %subplot(235)
        %[~,huas2] = contourf(yz/1000,squeeze(zas(xind(2),:,:)), ...
        %                     squeeze(uas(xind(2),:,:,tt)));
        %colorbar; caxis([-1 1]*0.1); ylim([-1000 0]);

        %subplot(236)
        %[~,huas3] = contourf(yz/1000,squeeze(zas(xind(3),:,:)), ...
        %                     squeeze(uas(xind(3),:,:,tt)));
        %colorbar; caxis([-1 1]*0.1); ylim([-1000 0]);

        %for tt=2:size(uas,4)
        %    runs.update_zeta(hzeta,t0+tt-1);
        %
        %    xind(1) = ix0(tt) + nrr * ceil(rr/runs.rgrid.dx);
        %    set(hlines(1),'XData',[1 1]*xind(1)*runs.rgrid.dx/1000);
        %
        %    set(huas1,'ZData',squeeze(uas(xind(1),:,:,tt)));
        %    set(huas2,'ZData',squeeze(uas(xind(2),:,:,tt)));
        %    set(huas3,'ZData',squeeze(uas(xind(3),:,:,tt)));
        %    pause(0.2);
        %end
    end
    toc(tictotal);
end