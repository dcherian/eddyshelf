function [handles] = animate_3d(runs, tind, opt, hax)
    stride = [1 1 1 1];

    if ~isfield(opt, 'csdreducepatch')
        opt.csdreducepatch = 0.1;
    end
    if ~isfield(opt, 'eddreducepatch')
        opt.eddreducepatch = 0.1;
    end
    if ~isfield(opt, 'finalize')
        opt.finalize = 0;
    end
    if ~isfield(opt, 'linefilter')
        opt.linefilter = 0;
    end
    if ~isfield(opt, 'csdcontours')
        opt.csdcontours = [runs.bathy.xsb + runs.sgntamp*5]*1000;
    end
    if ~isfield(opt, 'eddthresh')
        opt.eddthresh = 0.8;
    end
    if ~isfield(opt, 'sect')
        opt.sect = 'y';
    end
    if ~isfield(opt, 'MoveToZLevel')
        opt.MoveToZLevel = 0;
    end
    if ~isfield(opt, 'nolabels')
        opt.nolabels = 0;
    end

    tind = runs.process_time(tind);

    imx = runs.eddy.imx;
    imy = runs.eddy.imy;

    if length(tind) == 1
        xrange = imx(tind)-100:stride(1):imx(tind)+100;
        yrange = 1:stride(2):min(imy(tind)+100, size(runs.rgrid.xr,2));
    else
        xrange = 1:stride(1):size(runs.rgrid.xr,1);
        yrange = 1:stride(2):size(runs.rgrid.xr,2);
    end

    clear imx imy
    volume = {'x' xrange(1) xrange(end);
              'y' yrange(1) yrange(end)};

    xrmat = repmat(runs.rgrid.xr(xrange, yrange)', [1 1 runs.rgrid.N]);
    yrmat = repmat(runs.rgrid.yr(xrange, yrange)', [1 1 runs.rgrid.N]);
    zrmat = permute(runs.rgrid.z_r(:, yrange, xrange),[2 3 1]);

    csdye = dc_roms_read_data(runs, runs.csdname, tind, volume);
    eddye = dc_roms_read_data(runs, runs.eddname, tind, volume);
    csdye = permute(csdye,[2 1 3 4]);
    eddye = permute(eddye,[2 1 3 4]);

    mask = ones([size(runs.eddy.vor.mask,1)+2 ...
                  size(runs.eddy.vor.mask, 2)+2 length(tind)]);
    %mask(2:end-1,2:end-1,:) = runs.eddy.vor.mask(:,:,tind);
    mask = permute(mask(xrange,yrange,:), [2 1 3]); %mask = ones(size(xrmat));

    % chop off extra eddy dye by specifying a mask
    if opt.linefilter
        x = opt.x;
        y = opt.y;
        m = diff(y)./diff(x);
        c = y(2) - m * x(2);
        %assert(c == (y(1) - m*x(1)));
        mask(yrmat(:,:,1) - m*xrmat(:,:,1) - c > 0) = 0;
    end

    %% make isosurface plot
    eddlevel = zrmat; %opt.eddthresh;
    xsb = runs.bathy.xsb/1000;
    cslevel = opt.csdcontours;
    if ~isempty(cslevel)
        sbcolors = cbrewer('seq', 'Blues',8);
        sbcolors = sbcolors(end-length(cslevel):end,:);
    end

    if ~exist('hax', 'var') | isempty(hax)
        figure; maximize; hax = gca;
    else
        axes(hax);
    end

    clear hcsd hedd;
    hold on
    % first bathymetry
    hbathy = surf(runs.rgrid.xr/1000,runs.rgrid.yr/1000,-runs.bathy.h);
    set(hbathy,'FaceColor','Flat','EdgeColor','None');
    hbathy.FaceColor = [1 1 1];
    hbathy.FaceAlpha = 0.3;

    % add eddy. mask the dye field, mooth it and create a patch object.
    eddColorMap = runs.eddyeColormap;
    eddVolume = smooth3(bsxfun(@times,eddye(:,:,:,1),mask(:,:,1)));
    hedd = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, eddVolume, opt.eddthresh), ...
                 'EdgeColor', 'none', 'AmbientStrength', 0.5, 'DiffuseStrength', 0.6);
    % decimate patch
    reducepatch(hedd, opt.eddreducepatch, 'verbose');
    if opt.finalize
        % makes lighting better and smoother
        isonormals(eddVolume, hedd);
    end
    hedd.FaceColor = eddColorMap(end,:);
    hedd.FaceAlpha = 1;

    % add a cap to the top of the eddy.
    heddcap = patch(isocaps(xrmat/1000,yrmat/1000,zrmat, eddVolume, opt.eddthresh), ...
                    'EdgeColor', 'none');
    linkprop([hedd heddcap], {'FaceColor', 'FaceAlpha'});
    hedd.FaceColor = brighten([215 48 31]/255, 0.1);

    % Position light in dataunits
    hlight = camlight;
    hlight.Position(3) = -200;
    hlight.Position(2) = 1400;
    hlight.Position(1) = -738;

    sz = size(xrmat);

    % cross-sectional planes with what look like pcolor/contourf plots
    if ~isempty(cslevel)
        % combine the dyes to make one field that I'll plot as cross-section
        mergedye = eddVolume;
        mergedye(eddVolume < 0.6) = 0;
        mergedye = fillnan(smooth3(mergedye + -1 * (csdye <= cslevel(1)), 'gaussian', 5),0);

        if opt.sect == 'y'
            % x-z cross-secion
            imy = runs.eddy.imy(tind);

            handles.hsect = ...
                surface(squeeze(xrmat(imy,:,:)/1000), ...
                        squeeze(max(yrmat(:)/1000) * ones(size(yrmat(imy,:,:)))), ...
                        squeeze(zrmat(imy,:,:)), ...
                        squeeze(mergedye(imy,:,:)), ...
                        'EdgeColor', 'none', 'FaceColor', 'interp', ...
                        'AmbientStrength', 0.8);

            xvec = xrmat(imy,:,1);
            yvec = yrmat(imy,:,1);
            zvec = squeeze(zrmat(imy,1,:));

            handles.hsectoutline = plot3(xvec([1 sz(2) sz(2) 1 1])/1000, ...
                                         ones([1 5])*max(yrmat(:)/1000), ...
                                         zvec([1 1 sz(3) sz(3) 1]), ...
                                         'Color', [1 1 1]*0.3, ...
                                         'LineWidth', 1, 'Tag', 'dcline');

            handles.hplane = patch(xvec([1 sz(2) sz(2) 1 1])/1000, ...
                                   ones([1 5])*yvec(imy)/1000, ...
                                   [zvec([1 1])' 1 1 zvec(1)], 'k', ...
                                   'EdgeColor', [1 1 1]*0.3, ...
                                   'FaceAlpha', 0.08);
        else
            % y-z section
            imx = runs.eddy.imx(tind) - xrange(1) + 1;

            handles.hsect = ...
                surface(squeeze(max(xrmat(:)/1000) * ones(size(xrmat(:,imx,:)))), ...
                        squeeze(yrmat(:,imx,:)/1000), ...
                        squeeze(zrmat(:,imx,:)), ...
                        squeeze(mergedye(:,imx,:)), ...
                        'EdgeColor', 'none', 'FaceColor', 'interp', ...
                        'AmbientStrength', 0.8);

            xvec = xrmat(:,imx,1);
            yvec = yrmat(:,imx,1);
            zvec = squeeze(zrmat(end,imx,:));

            handles.hsectoutline = plot3(max(xrmat(:)/1000)*ones([1 5]), ...
                                         yvec([1 sz(1) sz(1) 1 1])/1000, ...
                                         zvec([1 1 sz(3) sz(3) 1]), ...
                                         'Color', [1 1 1]*0.3, ...
                                         'LineWidth', 1, 'Tag', 'dcline');

            handles.hplane = patch(ones([1 5]) * xvec(imx)/1000, ...
                                   yvec([1 sz(1) sz(1) 1 1])/1000, ...
                                   [zvec([1 1])' 1 1 zvec(1)], 'k', ...
                                   'EdgeColor', [1 1 1]*0.3, ...
                                   'FaceAlpha', 0.08);
        end
        colormap(flip(cbrewer('div', 'RdBu', 50)));
        caxis([-1.5 1.5]);
    end

    % shelf-slope water ribbon
    for kk=1:length(cslevel)
        hcsd(kk) = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                    smooth3(csdye),cslevel(kk)), ...
                         'FaceAlpha', 0.8, 'FaceColor', [107 174 214]/255, ...
                         'EdgeColor', 'none', 'AmbientStrength', 0.65);

        % be fancy and draw contours at surface to highlight downwelling
        cm = contourc(xrmat(1,:,end)/1000, yrmat(:,1,end)/1000, csdye(:,:,end), ...
                      [1 1] * cslevel(kk));
        ind = 1;
        indstart = 1;
        indend = 1;
        while indend < size(cm,2)
            NumPoints = cm(2,indend);
            indstart = indend+1;
            indend = indend + NumPoints+1;
            range = indstart:indend-1;

            xind = vecfind(xrmat(1,:,1)/1000, cm(1,range));
            yind = vecfind(yrmat(:,1,1)/1000, cm(2,range));
            clear zvec;
            for xx=1:length(xind)
                zvec(xx) = zrmat(yind(xx), xind(xx), end);
            end
            hcsdsurf{kk}(ind) = plot3(cm(1,range), cm(2,range), ...
                                      zvec, 'Color', runs.shelfSlopeColor('dark'), ...
                                      'LineWidth', 1, 'Tag', 'dcline');
            ind = ind+1;
        end

        if opt.finalize
            isonormals(smooth3(csdye(:,:,:,1)), hcsd(kk));
        end
        % decimate
        reducepatch(hcsd(kk), opt.csdreducepatch, 'verbose');
    end

    %runs.read_zeta(tind);
    %keyboard;
    %ssh = fillnan(runs.zeta(xrange,yrange,tind(1))' .* (eddVolume(:,:,end) > opt.eddthresh)*1e4, 0);
    %ssh = ssh - min(ssh(:)) + max(zrmat(:));
    %hzeta = surf(xrmat(:,:,1)/1000, yrmat(:,:,1)/1000, ssh, 'EdgeColor', 'none');
    if ~opt.nolabels
        titlestr = 'dyes';
        ht = runs.set_title(titlestr,tind(1));
        xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (m)');
    else
        hax.XAxis.Visible = 'off';
        hax.YAxis.Visible = 'off';
        hax.ZAxis.Visible = 'off';
        hax.Box = 'off';
    end
    xlim([min(xrmat(:)) max(xrmat(:))]/1000)
    ylim([min(yrmat(:)) max(yrmat(:))]/1000)
    beautify;
    view(-120,30);

    handles.heddcap = heddcap;
    handles.hedd = hedd;
    handles.hbathy = hbathy;
    handles.hlight = hlight;
    if ~isempty(cslevel)
        handles.hcsd = hcsd;
        handles.hcsdsurf = hcsdsurf;
    end

    % bring bottom up to have less empty space.
    if opt.MoveToZLevel ~= 0
        opt.MoveToZLevel = abs(opt.MoveToZLevel);
        bathy = hbathy.ZData;
        handles.hbathy.ZData(bathy <= -opt.MoveToZLevel) = -opt.MoveToZLevel;
        if ~isempty(cslevel)
            handles.hsectoutline.ZData(handles.hsectoutline.ZData <= -opt.MoveToZLevel) = ...
                -opt.MoveToZLevel;
            handles.hplane.ZData(handles.hplane.ZData <= -opt.MoveToZLevel) = ...
                -opt.MoveToZLevel;
        end
        zlim([-abs(opt.MoveToZLevel) 1.1]);
    end

    % repeat for animation
    for ii=2:4:size(eddye,4)
        pause();
        heddye = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                            bsxfun(@times,eddye(:,:,:,ii) > thresh,mask(:,:,ii)'),eddlevel);
        set(hedd,'Vertices',heddye.vertices,'Faces',heddye.faces, ...
                 'FaceVertexCData',heddye.facevertexcdata);
        set(hedd,'ZData',runs.zeta(:,:,ii));
        for kk=1:length(cslevel)
            hcsdye = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                csdye(:,:,:,ii),cslevel(kk));
            set(hcsd(kk),'Vertices',hcsdye.vertices,'Faces',hcsdye.faces);
        end
        runs.update_title(ht,titlestr,ii);
    end
end
