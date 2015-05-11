% animate variables at the surface

function [] = animate_surf(runs, varname, hax, t0, ntimes)

    warning('USE ANIMATE_FIELD INSTEAD');

    if ~exist('hax', 'var'), hax = []; end
    if ~exist('ntimes', 'var'), ntimes = length(runs.time); end
    if ~exist('t0', 'var'), t0 = 1; end

    dt = 4;

    runs.video_init(['surf-' varname]);

    % cross-shelf dye?
    if strcmpi(varname, runs.csdname) | strcmpi(varname, 'csdye')
        if isempty(runs.csdsurf) | isnan(runs.zeta(:,:,t0))
            if ntimes == 1
                tindices = t0;
            else
                tindices = [];
            end
            runs.csdsurf = dc_roms_read_data(runs.dir, runs.csdname, tindices, {'z' ...
                                runs.rgrid.N runs.rgrid.N}, [], runs.rgrid, ...
                                             'his');
        end

        vname = 'runs.csdsurf';
        varname = 'cross-shelf dye conc.';
    end

    % eddy dye?
    if strcmpi(varname, runs.eddname) | strcmpi(varname, 'eddye')
        if isempty(runs.eddsurf)
            if ntimes == 1
                tindices = t0;
            else
                tindices = [];
            end

            runs.eddsurf(:,:,tt) = dc_roms_read_data(runs.dir, runs.eddname, tindices, ...
                                                     {'z' runs.rgrid.N ...
                                runs.rgrid.N}, [], runs.rgrid, 'his');
        end

        vname = 'runs.eddsurf';
        varname = 'eddy dye conc.';
    end

    titlestr = ['Surface ' varname];

    xr = runs.rgrid.x_rho'/1000;
    yr = runs.rgrid.y_rho'/1000;
    tt = t0;

    if isempty(hax), figure; else axis(hax); end
    eval(['hpc = pcolorcen(xr, yr, ' vname '(:,:,tt));']);
    hold on; colorbar; freezeColors;

    he = runs.plot_eddy_contour('contour', tt);
    hbathy = runs.plot_bathy('contour', 'k');
    ht = runs.set_title(titlestr, tt);

    linex(runs.asflux.x/1000);

    if ntimes > 1
        runs.video_update();
        for tt=t0+1:dt:ntimes
            eval(['set(hpc, ''CData'', ' vname '(:,:,tt));']);
            runs.update_eddy_contour(he, tt);
            runs.update_title(ht, titlestr, tt);
            runs.video_update();
            pause(1);
        end
        runs.video_write();
    end
end