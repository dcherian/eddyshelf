%        [] = hovmoeller(runs, varname, axname, loc, iz, hax)
function [] = hovmoeller(runs, varname, axname, loc, iz, hax)
    if ~exist('loc', 'var') || isempty(loc)
        loc = [-30 -20 -10 0] * 1000 + runs.bathy.xsb;
        locu = [200 250 300 350]*1000;
    end

    if ~exist('iz', 'var') || isempty(iz)
        iz = runs.rgrid.N;
    end

    if ~exist('hax', 'var')
        figure;
        hax = gca;
    end

    handles.hax = hax; maximize;

    insertAnnotation([runs.name '.hovmoeller']);

    if strcmpi(loc, 'sb')
        loc = num2str(runs.bathy.xsb);
    end

    if strcmpi(varname, 'zeta') | strcmpi(varname, 'ubar') ...
            | strcmpi(varname, 'vbar')
        vol = {axname loc loc};
    else
        vol = {axname loc loc; 'z' iz iz};
    end

    [var, xax, yax] = dc_roms_read_data(runs, varname, [], vol);

    % plot (xvec,t) for yvec = constant
    if axname == 'y'
        xvec = xax;
        yvec = yax;
        if runs.bathy.axis == 'x'
            labx = 'X - X_{sb} (km)';
        else
            labx = 'X (km)';
        end
    else
        xvec = yax;
        yvec = xax;
        if runs.bathy.axis == 'y'
            labx = 'Y - Y_{sb} (km)';
        else
            labx = 'Y (km)';
        end
    end

    ind = find_approx(yvec, str2double(loc), 1);

    [xmat, tmat] = ndgrid(xvec/1000, runs.time/86400);
    if axname ~= runs.bathy.axis
        x0 = runs.bathy.xsb/1000;
    else
        x0 = runs.eddy.mx/1000;
        if x0 ~= 0
            if runs.bathy.axis == 'x'
                labx = 'Y - Y_{eddy} (km)';
            else
                labx = 'X - X_{eddy} (km)';
            end
        end
    end
    xmat = bsxfun(@minus, xmat, x0);

    % step topography phase speed estimate
    cp = runs.params.phys.f0 * runs.bathy.xsb * ...
         (runs.bathy.hsl - runs.bathy.hsb)./runs.bathy.hsl;

    axes(hax);
    if axname == 'y'
        handles.hplt = pcolorcen(xmat, tmat, var);
    else
        handles.hplt = pcolorcen(tmat, xmat, var);
    end

    handles.hcb = colorbar;
    if ~strcmpi(varname, 'rho'), center_colorbar; end
    hold on
    if axname == 'y'
        handles.hedd(1) = plot(runs.eddy.mx/1000 - x0, runs.eddy.t, 'k');
        handles.hedd(2) = plot(runs.eddy.vor.ee/1000 - x0, runs.eddy.t, 'k--');
        handles.hedd(3) = plot(runs.eddy.vor.we/1000 - x0, runs.eddy.t, 'k--');
        xlabel(labx);
        if runs.bathy.axis == 'x', linex(0, 'shelfbreak'); end
        ylabel('Time (days)');
    else
        handles.hedd(1) = plot(runs.eddy.t, runs.eddy.my/1000 - x0, 'k');
        handles.hedd(2) = plot(runs.eddy.t, runs.eddy.vor.ne/1000 - x0, 'k--');
        handles.hedd(3) = plot(runs.eddy.t, runs.eddy.vor.se/1000 - x0, 'k--');
        ylabel(labx);
        xlabel('Time (days)');
        if runs.bathy.axis == 'y', liney(0, 'shelfbreak'); end
    end
    handles.htitle = title([runs.name ' |  ' varname ' | ' axname ' = ' ...
                        num2str(str2double(loc)/1000, '%d') ' km']);
    legend(handles.hedd, ...
           {'Eddy center', 'Western, eastern edges'}, 'Location', 'SouthWest');
    beautify;

    % wave - estimate
    % limx = xlim;
    % xvec = limx(1):10:limx(2);
    % tvec = 0.5 + 1./0.02 .* (xvec * 1000)./runs.eddy.tscale;
    % plot(xvec, tvec);
end
