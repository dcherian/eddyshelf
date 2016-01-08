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
        opt.csdcontours = [rusn.bathy.xsb+5]*1000;
    end
    if ~isfield(opt, 'eddthresh')
        opt.eddthresh = 0.8;
    end

    tind = runs.process_time(tind);
    imx = runs.eddy.imx;

    if length(tind) == 1
        xrange = imx(tind)-100:stride(1):imx(tind)+100;
        yrange = 1:stride(2):size(runs.rgrid.xr,2);
    else
        xrange = 1:stride(1):size(runs.rgrid.xr,1);
        yrange = 1:stride(2):size(runs.rgrid.xr,2);
    end
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

    % chop off extra eddy dye
    if opt.linefilter
        x = opt.x;
        y = opt.y;
        m = diff(y)./diff(x);
        c = y(2) - m * x(2);
        assert(c == (y(1) - m*x(1)));
        mask(yrmat(:,:,1) - m*xrmat(:,:,1) - c > 0) = 0;
    end

    %% make isosurface plot
    eddlevel = zrmat; %opt.eddthresh;
    xsb = runs.bathy.xsb/1000;
    cslevel = opt.csdcontours;
    sbcolors = cbrewer('seq', 'Blues',8);
    sbcolors = sbcolors(end-length(cslevel):end,:);

    if ~exist('hax', 'var') | isempty(hax)
        figure; maximize;
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

    % add eddy
    eddColorMap = runs.eddyeColormap;
    eddVolume = smooth3(bsxfun(@times,eddye(:,:,:,1),mask(:,:,1)));
    hedd = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, eddVolume, opt.eddthresh), ...
                 'EdgeColor', 'none', 'AmbientStrength', 0.65, 'DiffuseStrength', 0.6);
    reducepatch(hedd, opt.eddreducepatch, 'verbose');
    if opt.finalize
        isonormals(eddVolume, hedd);
    end
    hedd.FaceColor = eddColorMap(end,:);
    hedd.FaceAlpha = 1;
    %axis tight;

    heddcap = patch(isocaps(xrmat/1000,yrmat/1000,zrmat, eddVolume, opt.eddthresh), ...
                    'EdgeColor', 'none');
    linkprop([hedd heddcap], {'FaceColor', 'FaceAlpha'});

    % heddfull = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
    %                             smooth3(bsxfun(@times,eddye(:,:,:,1),1 - mask(:,:,1)))), ...
    %                             'EdgeColor', 'none');
    % reducepatch(heddfull, 0.3, 'verbose');
    % heddfull.FaceAlpha = 0.7;
    % heddfull.FaceColor = eddColorMap(end-3,:);

    hlight = camlight;
    hlight.Position(3) = -200;
    hlight.Position(2) = 1400;
    hlight.Position(1) = -738;

    imy = runs.eddy.imy(tind) - 5;
    handles.hysect = ...
        surface(squeeze(xrmat(imy,:,:)/1000), ...
                squeeze(max(ylim) * ones(size(yrmat(imy,:,:)))), ...
                squeeze(zrmat(imy,:,:)), ...
                squeeze(fillnan(smooth3((eddVolume(imy,:,:) >= opt.eddthresh) + ...
                                        -1 * (csdye(imy,:,:) <= cslevel(1))),0)), ...
                'EdgeColor', 'none', 'FaceColor', 'interp', ...
                'AmbientStrength', 0.8);
    colormap(flip(cbrewer('div', 'RdBu', 10)));
    caxis([-1 1.5]);

    sz = size(xrmat);
    xvec = xrmat(imy,:,1);
    yvec = yrmat(imy,:,1);
    zvec = squeeze(zrmat(imy,1,:));

    handles.hysectoutline = plot3(xvec([1 sz(2) sz(2) 1 1])/1000, ones([1 5])*max(ylim), ...
                                  zvec([1 1 sz(3) sz(3) 1]), 'Color', [1 1 1]*0.3, ...
                                  'LineWidth', 1, 'Tag', 'dcline');

    handles.hyplane = patch(xvec([1 sz(2) sz(2) 1 1])/1000, ones([1 5])*yvec(imy)/1000, ...
                            zvec([1 1 sz(3) sz(3) 1]), 'k', 'EdgeColor', [1 1 1]*0.3, ...
                            'FaceAlpha', 0.08);

    for kk=1:length(cslevel)
        hcsd(kk) = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                    smooth3(csdye(:,:,:,1)),cslevel(kk)));
        if opt.finalize
            isonormals(smooth3(csdye(:,:,:,1)), hcsd(kk));
        end
        hcsd(kk).FaceColor = sbcolors(kk,:);
        hcsd(kk).EdgeColor = 'none';
        hcsd(kk).FaceAlpha = 0.5;
        hcsd(kk).AmbientStrength = 0.65;
        reducepatch(hcsd(kk), opt.csdreducepatch, 'verbose');
    end

    %runs.read_zeta(tind);
    %keyboard;
    %ssh = fillnan(runs.zeta(xrange,yrange,tind(1))' .* (eddVolume(:,:,end) > opt.eddthresh)*1e4, 0);
    %ssh = ssh - min(ssh(:)) + max(zrmat(:));
    %hzeta = surf(xrmat(:,:,1)/1000, yrmat(:,:,1)/1000, ssh, 'EdgeColor', 'none');
    titlestr = 'dyes';
    ht = runs.set_title(titlestr,tind(1));
    %view(-104,30);
    %view(-150,66);
    %view(2)
    %view(-50,50); zoom(1.2);
    xlim([min(xrmat(:)) max(xrmat(:))]/1000)
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (m)');
    beautify;

    view(-120,30);

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

    handles.heddcap = heddcap;
    handles.hedd = hedd;
    handles.hcsd = hcsd;
    handles.hbathy = hbathy;
    handles.hlight = hlight;
end
