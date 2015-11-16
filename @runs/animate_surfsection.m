% animate variables at the surface and cross-section
function [] = animate_surfsection(runs, varname, varname1, t0, ntimes)

    dx = 60; % start from shelfbreak - dx km
    dt = 3;

    if ~exist('ntimes', 'var')
        ntimes = length(runs.time);
    else
        ntimes = t0 + dt * (ntimes - 1);
    end
    if ~exist('t0', 'var'), t0 = 1; end

    if isempty(varname1)
        runs.video_init(['section-' varname]);
    else
        runs.video_init(['section-' varname '-' varname1]);
    end
    makeVideo = runs.makeVideo;

    y0 = runs.eddy.my;
    % y0 = (runs.bathy.xsb+30e3) * ones(size(runs.time));
    csfluxflag = 0;
    if csfluxflag
        isobath = 1;
        y1 = runs.csflux.x(isobath) * ones(size(runs.time));
    else
        y1 = y0 - runs.sgntamp * runs.eddy.vor.dia/4;
        %y1 = (runs.bathy.xsb+15e3) * ones(size(runs.time));
    end
    iy = vecfind(runs.rgrid.y_rho(:,1), y0);
    iy1 = vecfind(runs.rgrid.y_rho(:,1), y1);

    tt = t0;
    x0 = runs.eddy.mx(tt)/1000;
    xr = runs.rgrid.x_rho(1,:)'/1000;
    xvec = xr(2:end-1) - x0;
    zmat = runs.rgrid.z_r;

    varname = runs.process_varname(varname);

    % read data here, so that I don't incur overhead at least for
    % surface fields.

    % cross-shelf dye?
    if strcmpi(varname, runs.csdname);
        runs.read_csdsurf;
        varname = runs.csdname;
    end

    % eddy dye?
    if strcmpi(varname, runs.eddname)
        runs.read_eddsurf;
        varname = runs.eddname;
    end

    if ~exist('varname1', 'var') || isempty(varname1)
        varname1 = varname;
    end

    % if different variables, do same location
    if ~strcmpi(varname1, varname), iy1 = iy; y1 = y0; end

    % process first variable
    v = dc_roms_read_data(runs.dir, varname, ...
                          tt, {runs.bathy.axis iy(tt) iy(tt)}, ...
                          [], runs.rgrid, 'his', 'single');
    v = v(2:end-1,:,:);
    if strcmpi(varname, 'dye_01'), v = v / 1000; end

    % process second variable
    if strcmpi(varname1, 'pv') || strcmpi(varname1, 'rv')
        xvec1 = avg1(xr,1) - x0;
        zmat1 = permute(ncread([runs.dir '/ocean_vor.nc'], ['z_' varname1]), ...
                        [3 2 1]);
        tpv = ncread([runs.dir '/ocean_vor.nc'], 'ocean_time');
        itpv = find_approx(tpv, runs.time(tt));
        v1 = dc_roms_read_data(runs.dir, varname1, itpv, ...
                               {runs.bathy.axis iy1(tt) iy1(tt)}, ...
                               [], runs.rgrid, 'his', 'single');
    else
        xvec1 = xvec;
        zmat1 = runs.rgrid.z_r;
        v1 = dc_roms_read_data(runs.dir, varname1, ...
                               tt, {runs.bathy.axis iy1(tt) iy1(tt)}, ...
                               [], runs.rgrid, 'his', 'single');
        v1 = v1(2:end-1,:,:);
    end
    if strcmpi(varname1, 'dye_01'), v1 = v1 / 1000; end

    % remove background?
    if strcmpi(varname, 'rho') || strcmpi(varname, 'dye_02')
        rback = dc_roms_read_data(runs.dir, varname, 1, {}, [], ...
                                  runs.rgrid, 'his', 'single');
        v = v - squeeze(rback(2:end-1,iy(tt),:));
        %v1 = bsxfun(@minus, v1, v1(1,:));
    end
    if strcmpi(varname1, 'rho') || strcmpi(varname1, 'dye_02')
        rback = dc_roms_read_data(runs.dir, varname1, 1, {}, [], ...
                                  runs.rgrid, 'his', 'single');
        v1 = v1 - squeeze(rback(2:end-1,iy1(tt),:));
        %v1 = bsxfun(@minus, v1, v1(1,:));
    end

    figure
    hax(1) = subplot(2,2,[1 2]);
    runs.makeVideo = 0;
    runs.animate_field(varname, hax, t0, 1);
    runs.makeVideo = makeVideo;
    ylim([runs.bathy.xsb/1000-dx max(ylim)]);
    if strcmpi(varname1, varname)
        liney([y1(tt) y0(tt)]/1000, {'1'; '2'}, 'k');
    else
        liney([y1(tt) y0(tt)]/1000, [], 'k');
    end
    climsurf = caxis;
    if strcmpi(varname, 'csdye') || strcmpi(varname, 'dye_01')
        if runs.sgntamp == 1
            climsurf(1) = runs.bathy.xsb/1000 - dx;
        else
            climsurf(2) = runs.bathy.xsb/1000 + dx;
        end
        climsurf = sort(climsurf);
    end
    caxis(climsurf);

    hax(2) = subplot(2,2,4);
    zvec = zmat(:, iy(tt)+1, 1);
    hplt = pcolor(xvec, zvec, v');
    hold on
    if strcmpi(varname, 'rho')
        clim = caxis;
        hrho = contour(xvec, zvec, v', ...
                       [1 1]*runs.eddy.drhothresh(1), 'k', ...
                       'LineWidth', 2);
        caxis(clim);
    else
        if strcmpi(varname1, 'rho')
            clim = caxis;
            hrho = contour(xvec, zvec, v1', ...
                           [1 1]*runs.eddy.drhothresh(1), 'k', ...
                           'LineWidth', 2);
            caxis(clim);
        end
    end
     if strcmpi(varname, 'rho') || strcmpi(varname, 'dye_02')
        center_colorbar;
    end
    hl = linex([runs.eddy.rhovor.ee(tt) ...
                runs.eddy.rhovor.we(tt) x0]/1000 - x0, [], 'k');
    hly = liney(-1 * runs.eddy.Lgauss(tt), [], 'k');
    liney(-1*runs.bathy.hsb, 'shelfbreak');
    if strcmpi(varname1, varname)
        ht = title([runs.bathy.axis ' = ' num2str(y0(tt)/1000, '%.0f') ' km']);
    else
        title(varname);
    end
    shading interp;
    ylim([-1*min(runs.bathy.h(1,iy1)) 0]);
    colorbar;
    if ~strcmpi(varname, 'rho')
        caxis(climsurf);
    end


    hax(3) = subplot(2,2,3);
    zvec1 = zmat1(:, iy1(tt)+1, 1);
    hplt1 = pcolor(xvec1, zvec1, v1'); hold on;
    hl1 = linex([runs.eddy.rhovor.ee(tt) ...
                 runs.eddy.rhovor.we(tt) x0]/1000 - x0, [], 'k');
    hly1 = liney(-1 * runs.eddy.Lgauss(tt), [], 'k');
    if strcmpi(varname1, 'rho')
        clim1 = caxis;
        hrho = contour(xvec, zvec, v1', ...
                       [1 1]*runs.eddy.drhothresh(1), 'k', ...
                       'LineWidth', 2);
        caxis(clim1);
    end
    if strcmpi(varname1, 'rho') || strcmpi(varname1, 'dye_02')
        center_colorbar;
    end
    liney(-1*runs.bathy.hsb, 'shelfbreak');
    if strcmpi(varname1, varname)
        ht1 = title([runs.bathy.axis ' = ' num2str(y1(tt)/1000, '%.0f') ' km']);
    else
        title(varname1);
    end
    shading interp;
    ylim([-1*min(runs.bathy.h(1,iy)) 0]);
    colorbar;
    if strcmpi(varname1, varname) && ~strcmpi(varname, 'rho')
        caxis(climsurf);
    end

    % linkaxes(hax, 'x');
    xlim([-1 1]* 200);
    ylim([-300 0]);
    linkaxes(hax([3 2]), 'xy');

    if ntimes > 1
        runs.video_update();

        for tt=t0+1:dt:ntimes
            axes(hax(1))
            cla(hax(1))
            runs.makeVideo = 0;
            runs.animate_field(varname, hax, tt, 1);
            runs.makeVideo = makeVideo;
            if strcmpi(varname1, varname)
                liney([y1(tt) y0(tt)]/1000, {'1'; '2'}, 'k');
            else
                liney([y1(tt) y0(tt)]/1000, [], 'k');
            end
            ylim([runs.bathy.xsb/1000-dx max(ylim)])
            caxis(climsurf);

            v = dc_roms_read_data(runs.dir, varname, ...
                                  tt, {runs.bathy.axis iy(tt) iy(tt)}, ...
                                  [], runs.rgrid, 'his', 'single');
            v1 = dc_roms_read_data(runs.dir, varname1, ...
                                   tt, {runs.bathy.axis iy1(tt) iy1(tt)}, ...
                                   [], runs.rgrid, 'his', 'single');
            v = v(2:end-1,:,:);
            v1 = v1(2:end-1,:,:);

            if strcmpi(varname, 'rho') || strcmpi(varname, 'dye_02')
                v = v - squeeze(rback(2:end-1,iy1(tt),:));
            end

            if strcmpi(varname1, 'rho') || strcmpi(varname1, 'dye_02')
                v1 = v1 - squeeze(rback(2:end-1,iy1(tt),:));
            end

            if strcmpi(varname, 'dye_01'), v = v / 1000; end
            if strcmpi(varname1, 'dye_01'), v1 = v1 / 1000; end

            x0 = runs.eddy.mx(tt)/1000;
            xvec = xr(2:end-1) - x0;
            if strcmpi(varname1, 'pv') || strcmpi(varname1, 'rv')
                xvec1 = avg1(xr,1) - x0;
            else
                xvec1 = xvec;
            end
            zvec = zmat(:, iy(tt)+1, 1);
            zvec1 = zmat1(:, iy1(tt)+1, 1);

            hplt1.CData = v1';
            hplt1.XData = xvec1;
            hplt1.YData = zvec1;
            hl1{1}.XData = [1 1] * runs.eddy.rhovor.ee(tt)/1000 - x0;
            hl1{2}.XData = [1 1] * runs.eddy.rhovor.we(tt)/1000 - x0;
            hly1.YData = [1 1] * -1 * runs.eddy.Lgauss(tt);

            hplt.CData = v';
            hplt.XData = xvec;
            hplt.YData = zvec;
            hl{1}.XData = [1 1] * runs.eddy.rhovor.ee(tt)/1000 - x0;
            hl{2}.XData = [1 1] * runs.eddy.rhovor.we(tt)/1000 - x0;
            hly.YData = [1 1] * -1 * runs.eddy.Lgauss(tt);

            ht.String = [runs.bathy.axis ' = ' num2str(y0(tt)/1000, '%.0f') ' km'];
            if strcmpi(varname, varname1)
                ht1.String = [runs.bathy.axis ' = ' num2str(y1(tt)/1000, '%.0f') ' km'];
            end
            runs.video_update();
            pause(1);
        end
        runs.video_write();
    end
end