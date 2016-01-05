function [handles] = animate_3d(runs, tind, hax)
    stride = [1 1 1 1];

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
    x = [300, 375] * 1000;
    y = [300, 0] * 1000;
    m = diff(y)./diff(x);
    c = y(2) - m * x(2);
    assert(c == (y(1) - m*x(1)));
    mask(yrmat(:,:,1) - m*xrmat(:,:,1) - c > 0) = 0;

    %% make isosurface plot
    eddlevel = zrmat; %0.8;
    thresh = 0.8;
    xsb = runs.bathy.xsb/1000;
    cslevel = [xsb+5]*1000;
    sbcolors = cbrewer('seq', 'Blues',8);
    sbcolors = sbcolors(end-length(cslevel)+1:end,:);

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
    hbathy.FaceColor = [1 1 1] * 0.65;
    hbathy.FaceAlpha = 0.5;

    % add eddy
    eddColorMap = runs.eddyeColormap;
    hedd = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                            smooth3(bsxfun(@times,eddye(:,:,:,1),mask(:,:,1))), ...
                            0.8), 'EdgeColor', 'none');
    reducepatch(hedd, 0.5, 'verbose');
    hedd.FaceColor = eddColorMap(end,:);
    hedd.FaceAlpha = 1;
    axis tight;

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

    %hedd = patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors, ...
    %             'FaceColor','interp','EdgeColor','none');
    %colormap(flipud(cbrewer('div', 'RdYlGn', 32))); freezeColors;
    %colorbar;
    %set(hedd,'EdgeColor','none','FaceAlpha',0.5);

    for kk=1:length(cslevel)
        hcsd(kk) = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                    smooth3(csdye(:,:,:,1)),cslevel(kk)));
        hcsd(kk).FaceColor = sbcolors(kk,:);
        hcsd(kk).EdgeColor = 'none';
        hcsd(kk).FaceAlpha = 0.5;
        reducepatch(hcsd(kk),0.1,'verbose');
    end

    %runs.read_zeta(tind);
    %[~,hzeta] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,tind(1)));
    titlestr = 'dyes';
    ht = runs.set_title(titlestr,tind(1));
    %view(-104,30);
    %view(-150,66);
    %view(2)
    %view(-50,50); zoom(1.2);
    xlim([min(xrmat(:)) max(xrmat(:))]/1000)
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (m)');
    beautify;

    view(-120,30); zoom(1.2);

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

    handles.hedd = hedd;
    handles.hcsd = hcsd;
    handles.hbathy = hbathy;
    handles.hlight = hlight;
end
