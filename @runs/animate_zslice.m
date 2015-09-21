function [] = animate_zslice(runs,varname,depth,tind)

    if ~exist('tind','var'), tind = []; end
    [~,tind,~,nt,stride] = roms_tindices(tind,Inf,length(runs.time));

    runs.video_init(['z' num2str(depth) '-' varname]);
    if strcmp(varname,'vor');
        grids = [runs.dir '/ocean_vor.nc'];
    else
        grids = runs.rgrid;
    end

    [grd.xax,grd.yax,grd.zax,~] = dc_roms_extract(grids,varname,{},1);
    datain = 0;
    if nt < 20
        tic; disp('Reading data...');
        data = dc_roms_read_data(runs.dir, varname, tind, {}, [], ...
                                     grids);
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
            disp(['interpolating timestep ' num2str(mmm) '/' ...
                  num2str(nt)]);
            var(:,:,mmm) = dc_roms_zslice_var(data(:,:,:,mmm),depth,grd);
        end
    end
    clear data

    runs.read_csdsurf;
    runs.read_eddsurf;

    % animate
    figure; maximize;
    insertAnnotation([runs.name '.animate_zslice']);
    xax = grd.xax(:,:,1)/1000; yax=  grd.yax(:,:,1)/1000; clear grd;
    tt = 1;
    [~,hc] = contourf(xax,yax,var(:,:,tt), 10, 'LineWidth', 0);
    hold on
    if strcmpi(varname, 'rho')
        [~,hrho] = contour(xax, yax, var(:,:,tt), ...
                           [1 1]*runs.eddy.drhothresh(1), ...
                           'Color', 'b', 'LineWidth', 2);
    end
    he = runs.plot_rho_contour('contour',tind(1) + tt-1);
    [~,hcsd] = contour(xax, yax, runs.csdsurf(:,:,tind(1)+tt-1), ...
                       runs.csflux.x([1 4 6]), ...
                       'Color', [1 1 1]*0.55, 'LineWidth', 2);
    % [~,hedd] = contour(xax, yax, runs.eddsurf(:,:,tind(1)+tt-1), ...
    %                    [0.9 0.9], 'Color', 'k', 'LineWidth', 2);
    shading flat;
    ht = title([varname ' | z = ' num2str(depth) ' m | t = ' ...
                num2str(runs.time(tind(1)+tt-1)/86400) ' days']);
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

    for tt=2:nt
        runs.video_update;
        pause(0.1);
        set(hc,'ZData',var(:,:,tt));
        shading flat
        runs.update_rho_contour(he,tind(1) + tt-1);
        hcsd.ZData = runs.csdsurf(:,:,tind(1) + tt-1);
        set(ht,'String',[varname ' | z = ' num2str(depth) ' m | t = ' ...
                         num2str(runs.time(tind(1)+tt-1)/86400) 'days']);
        if strcmpi(varname, 'rho')
            hrho.ZData = var(:,:,tt);
        end
    end
    runs.video_write;
end
