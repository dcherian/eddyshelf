function [] = animate_3d(runs, tind)
    stride = [1 1 1 1];

    xrmat = repmat(runs.rgrid.xr(1:stride(1):end,1:stride(2):end)', ...
                   [1 1 runs.rgrid.N]);
    yrmat = repmat(runs.rgrid.yr(1:stride(1):end,1:stride(2):end)', ...
                   [1 1 runs.rgrid.N]);
    zrmat = permute(runs.rgrid.zr(1:stride(1):end,1:stride(2):end,:),[2 1 3]);

    tic; csdye = ncread(runs.out_file,runs.csdname); toc;
    tic; eddye = ncread(runs.out_file,runs.eddname); toc;
    csdye = permute(csdye,[2 1 3 4]);
    eddye = permute(eddye,[2 1 3 4]);
    mask = zeros(size(eddye,2),size(eddye,1),size(eddye,4));
    mask(2:end-1,2:end-1,:)=runs.eddy.vor.mask;
    mask = 1 + zeros(size(mask));

    %% make isosurface plot

    eddlevel = zrmat; %0.8;
    thresh = 0.8;
    xsb = runs.bathy.xsb/1000;
    cslevel = [xsb-10 xsb]*1000;
    sbcolors = distinguishable_colors(length(cslevel));

    clf; clear pcsd pedd;
    hold on
    hbathy = surf(runs.rgrid.xr/1000,runs.rgrid.yr/1000,-runs.bathy.h);
    colormap(copper); freezeColors;
    set(hbathy,'FaceColor','Flat','EdgeColor','None');

    ii=1;
    [faces,verts,colors] = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                      bsxfun(@times,eddye(:,:,:,1) > thresh,mask(:,:,1)'),eddlevel);
    pedd = patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors, ...
                 'FaceColor','interp','EdgeColor','none');
    colormap(flipud(cbrewer('div', 'RdYlGn', 32))); freezeColors;
    colorbar; cbfreeze;
    %set(pedd,'EdgeColor','none','FaceAlpha',0.5);
    view(3)

    for kk=1:length(cslevel)
        pcsd(kk) = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                    csdye(:,:,:,ii),cslevel(kk)));
        set(pcsd(kk),'FaceColor',sbcolors(kk,:));
        set(pcsd(kk),'EdgeColor','none');
        %reducepatch(pcsd(kk),0.5,'verbose');
    end
    [~,hedd] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,1));

    titlestr = 'dyes';
    ht = runs.set_title(titlestr,ii);
    %view(-104,30);
    view(-150,66);
    xlim([min(xrmat(:)) max(xrmat(:))]/1000)
    xlabel('X'); ylabel('Y'); zlabel('Z');
    for ii=2:4:size(eddye,4)
        pause();
        heddye = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                            bsxfun(@times,eddye(:,:,:,ii) > thresh,mask(:,:,ii)'),eddlevel);
        set(pedd,'Vertices',heddye.vertices,'Faces',heddye.faces, ...
                 'FaceVertexCData',heddye.facevertexcdata);
        set(hedd,'ZData',runs.zeta(:,:,ii));
        for kk=1:length(cslevel)
            hcsdye = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                                csdye(:,:,:,ii),cslevel(kk));
            set(pcsd(kk),'Vertices',hcsdye.vertices,'Faces',hcsdye.faces);
        end
        runs.update_title(ht,titlestr,ii);
    end
end
