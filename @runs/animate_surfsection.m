% animate variables at the surface and cross-section
function [] = animate_surfsection(runs, varname, t0, ntimes)

    if ~exist('ntimes', 'var'), ntimes = length(runs.time); end
    if ~exist('t0', 'var'), t0 = 1; end

    dt = 4;

    runs.video_init(['surfsection-' varname]);
    makeVideo = runs.makeVideo;

    iy = vecfind(runs.rgrid.y_rho(:,1), runs.eddy.my);
    tt = t0;

    % read data here, so that I don't incur overhead at least for
    % surface fields.
    % cross-shelf dye?
    if strcmpi(varname, runs.csdname) | strcmpi(varname, 'csdye')
        runs.read_csdsurf;

        % process cross-shelf dye
        var = dc_roms_read_data(runs.dir, runs.csdname, ...
                                  tt, {runs.bathy.axis iy(tt) iy(tt)}, ...
                                  [], runs.rgrid, 'his', 'single')/1000;
    end

    % eddy dye?
    if strcmpi(varname, runs.eddname) | strcmpi(varname, 'eddye')
        runs.read_eddsurf;

        var = dc_roms_read_data(runs.dir, runs.eddname, ...
                                  tt, {runs.bathy.axis iy(tt) iy(tt)}, ...
                                  [], runs.rgrid, 'his', 'single');
    end
    var = var(2:end-1,:,:);

    figure
    hax(1) = subplot(211);
    runs.makeVideo = 0;
    runs.animate_field(varname, hax, t0, 1);
    runs.makeVideo = makeVideo;
    ylim([runs.bathy.xsb/1000-20 max(ylim)])
    liney(runs.eddy.my(tt)/1000, [], 'k');

    hax(2) = subplot(212);
    L = runs.eddy.rhovor.dia(tt)/2;
    xvec = runs.rgrid.x_rho(1,2:end-1);
    zvec = runs.rgrid.z_r(:, iy(tt)+1, 1);
    hplt = pcolor(xvec/1000, zvec, var' - runs.eddy.my(tt)/1000);
    hl = linex([runs.eddy.rhovor.ee(tt) runs.eddy.rhovor.we(tt)]/1000, ...
               [], 'k');
    shading interp;
    ylim([min(runs.rgrid.z_r(:)) 0]);
    colorbar;

    linkaxes(hax, 'x');

    if ntimes > 1
        runs.video_update();

        for tt=t0+1:dt:ntimes
            axes(hax(1))
            cla(hax(1))
            runs.makeVideo = 0;
            runs.animate_field(varname, hax, tt, 1);
            runs.makeVideo = makeVideo;
            liney(runs.eddy.my(tt)/1000, [], 'k');

            axes(hax(2));
            zvec = runs.rgrid.z_r(:, iy(tt)+1, 1);

            % cross-shelf dye?
            if strcmpi(varname, runs.csdname) | strcmpi(varname, 'csdye')
                % process cross-shelf dye
                var = dc_roms_read_data(runs.dir, runs.csdname, ...
                                        tt, {runs.bathy.axis iy(tt) iy(tt)}, ...
                                        [], runs.rgrid, 'his', 'single')/1000;
            end

            % eddy dye?
            if strcmpi(varname, runs.eddname) | strcmpi(varname, 'eddye')
                var = dc_roms_read_data(runs.dir, runs.eddname, ...
                                        tt, {runs.bathy.axis iy(tt) iy(tt)}, ...
                                        [], runs.rgrid, 'his', 'single');
            end
            var = var(2:end-1,:,:);

            hplt.CData = var';
            hplt.YData = zvec;
            hl{1}.XData = [1 1] * runs.eddy.vor.ee(tt)/1000;
            hl{2}.XData = [1 1] * runs.eddy.vor.we(tt)/1000;

            runs.video_update();
            pause(1);
        end
        runs.video_write();
    end
end