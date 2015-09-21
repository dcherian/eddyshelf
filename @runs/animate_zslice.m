function [] = animate_zslice(runs,varname,depth,tind)

    if ~exist('tind','var'), tind = []; end
    [~,tind,~,nt,stride] = roms_tindices(tind,Inf,length(runs.time));

    % quiver options
    addvelquiver = 1;
    dxi = 10; dyi = 3;
    uref = runs.eddy.V(1)/10; vref = uref;

    runs.video_init(['z' num2str(depth) '-' varname]);
    if strcmpi(varname,'rv') || strcmpi(varname, 'pv')
        grids = [runs.dir '/ocean_vor.nc'];
    else
        grids = runs.rgrid;
    end

    [grd.xax,grd.yax,grd.zax,~] = dc_roms_extract(grids,varname,{},1);

    if addvelquiver
        grdu.xax = repmat(runs.rgrid.x_u',[1 1 runs.rgrid.N]);
        grdu.yax = repmat(runs.rgrid.y_u',[1 1 runs.rgrid.N]);
        grdu.zax = permute(runs.rgrid.z_u,[3 2 1]);

        grdv.xax = repmat(runs.rgrid.x_v',[1 1 runs.rgrid.N]);
        grdv.yax = repmat(runs.rgrid.y_v',[1 1 runs.rgrid.N]);
        grdv.zax = permute(runs.rgrid.z_v,[3 2 1]);
    end

    datain = 0;
    if nt < 20
        tic; disp('Reading data...');
        data = dc_roms_read_data(runs.dir, varname, tind, {}, [], grids);
        if addvelquiver
            u1 = dc_roms_read_data(runs.dir, 'u', tind, {}, [], grids);
            v1 = dc_roms_read_data(runs.dir, 'v', tind, {}, [], grids);
        end

        datain = 1;
        var = nan([size(data,1) size(data,2) nt]);
        toc;
    end

    % read data
    for mmm = 1:nt
        if ~datain
            disp(['reading & interpolating timestep ' num2str(mmm) '/' ...
                  num2str(nt)]);

            data = dc_roms_read_data(runs.dir, varname, mmm, {}, [], ...
                                     grids);
            if strcmpi(varname, 'rho')
                rback = permute(dc_roms_read_data(runs.dir, 'rho', 1, ...
                                          {'x' 1 1}, [], grids), [3 1 2]);
                data = bsxfun(@minus, data, rback);
            end
            if mmm == 1
                var = nan([size(data,1) size(data,2) nt]);
            end
            var(:,:,mmm) = dc_roms_zslice_var(data,depth,grd);
        else
            if strcmpi(varname, 'rho')
                rback = permute(dc_roms_read_data(runs.dir, 'rho', 1, ...
                                          {'x' 1 1}, [], grids), [3 1 2]);
                data = bsxfun(@minus, data, rback);
            end
            disp(['interpolating timestep ' num2str(mmm) '/' num2str(nt)]);
            var(:,:,mmm) = dc_roms_zslice_var(data(:,:,:,mmm),depth,grd);
            if addvelquiver
                u(:,:,mmm) = dc_roms_zslice_var(u1(:,:,:,mmm),depth,grdu);
                v(:,:,mmm) = dc_roms_zslice_var(v1(:,:,:,mmm),depth,grdv);
            end
        end
    end
    clear data u1 v1


    % get on interior RHO points
    u = avg1(u(:,2:end-1,:),1);
    v = avg1(v(2:end-1,:,:),2);
    % decimate for quiver
    u = u(1:dxi:end,1:dyi:end,:);
    v = v(1:dxi:end,1:dyi:end,:);

    if strcmpi(varname, 'pv')
        var = log10(var);
    end

    runs.read_csdsurf;
    runs.read_eddsurf;

    titlestr = [varname ' | z = ' num2str(depth) ' m '];

    % animate
    figure; maximize;
    insertAnnotation([runs.name '.animate_zslice']);
    xax = grd.xax(:,:,1)/1000; yax = grd.yax(:,:,1)/1000; clear grd;
    tt = 1;
    [~,hc] = contourf(xax,yax,var(:,:,tt), 30);
    hc.EdgeColor = 'none';
    hold on
    if strcmpi(varname, 'rho')
        [~,hrho] = contour(xax, yax, var(:,:,tt), ...
                           [1 1]*runs.eddy.drhothresh(1), ...
                           'Color', 'b', 'LineWidth', 2);
    end
    he = runs.plot_rho_contour('contour',tind(1) + tt-1);
    [~,hcsd] = contour(runs.rgrid.x_rho'/1e3, runs.rgrid.y_rho'/1e3, ...
                       runs.csdsurf(:,:,tind(1)+tt-1), runs.csflux.x([1 4 6]), ...
                       'Color', [1 1 1]*0.55, 'LineWidth', 2);
    [~,hedd] = contour(runs.rgrid.x_rho'/1e3, runs.rgrid.y_rho'/1e3, ...
                       runs.eddsurf(:,:,tind(1)+tt-1), [0.9 0.9], ...
                       'Color', 'k', 'LineWidth', 2);
    shading flat;
    ht = runs.set_title(titlestr, tind(1)+tt-1);
    axis image;
    xlim([min(xax(:)) max(xax(:))]);
    ylim([min(yax(:)) max(yax(:))]);
    colorbar;
    if strcmpi(varname, runs.eddname)
        center_colorbar;
        caxis([-1 1]);
    else
        caxis([min(var(:)) max(var(:))]);
    end
    xlabel('X (km)'); ylabel('Y (km)');
    runs.plot_bathy('contour','k');

    hq = quiver(runs.eddy.xr(1:dxi:end,1:dyi:end)/1000, ...
                runs.eddy.yr(1:dxi:end,1:dyi:end)/1000, ...
                u(:,:,tt)./uref, v(:,:,tt)./vref, ...
                'Color', 'k', 'LineWidth', 2);

    for tt=2:nt
        runs.video_update;
        pause(0.1);
        set(hc,'ZData',var(:,:,tt));
        set(hq,'UData',u(:,:,tt));
        set(hq,'VData',v(:,:,tt));
        shading flat
        runs.update_rho_contour(he,tind(1) + tt-1);
        hcsd.ZData = runs.csdsurf(:,:,tind(1) + tt-1);
        runs.update_title(ht, titlestr, tind(1)+tt-1);
        if strcmpi(varname, 'rho')
            hrho.ZData = var(:,:,tt);
        end
    end
    runs.video_write;
end
