function [] = plot_xzsection(runs, loc, day, debug_flux)

    if ~exist('debug_flux', 'var')
        debug_flux = 0; % debug flux parameterization
    end

    if ~exist('day', 'var'), day = []; end
    if ischar(day)
        day = str2double(day);
        tindex = vecfind(runs.time, day*86400);
    else
        tindex = day;
    end

    if ~exist('loc', 'var')
        loc = runs.eddy.my(tindex);
    end
    if ischar(loc)
        loc = str2double(loc);
    else
        if loc <= length(runs.csflux.x)
            % should be integer isobath
            isobath = loc;
            loc = runs.csflux.x(isobath);
            if isempty(tindex)
                [~,tindex] = runs.calc_maxflux(...
                    runs.csflux.off.slope(:,isobath,isobath));;
            end
        else
            if runs.bathy.axis == 'y'
                loc = runs.rgrid.y_rho(loc,1);
            else
                loc = runs.rgrid.x_rho(1,loc);
            end
        end
    end
    ix = find_approx(runs.rgrid.y_rho(:,1), loc);

    L = runs.eddy.rhovor.dia(tindex)/2;

    % copied from csfluxes.m
    if runs.csvelname == 'v'
        bathyax = 2;
        xvec = runs.rgrid.x_rho(1,2:end-1) - runs.eddy.mx(tindex);
        zvec = runs.rgrid.z_r(:, ix, 1);
    else
        bathyax = 1;
        xvec = (runs.rgrid.y_rho(2:end-1,1) - runs.eddy.my(tindex))';
        zvec = runs.rgrid.z_r(:, 1, ix);
    end

    csvel = squeeze(avg1( ...
        dc_roms_read_data(runs.dir, runs.csvelname, tindex, ...
                          {runs.bathy.axis ix-1 ix}, [], runs.rgrid, ...
                          'his', 'single'), bathyax));
    csvel = csvel(2:end-1,:,:,:);

    % process cross-shelf dye
    csdye = dc_roms_read_data(runs.dir, runs.csdname, ...
                              tindex, {runs.bathy.axis ix+1 ix+1}, ...
                              [], runs.rgrid, 'his', 'single');
    csdye = csdye(2:end-1,:,:);

    eddye = dc_roms_read_data(runs.dir, runs.eddname, ...
                              tindex, {runs.bathy.axis ix+1 ix+1}, ...
                              [], runs.rgrid, 'his', 'single');
    eddye = eddye(2:end-1,:,:);

    rho = dc_roms_read_data(runs.dir, 'rho', ...
                            tindex, {runs.bathy.axis ix+1 ix+1}, ...
                            [], runs.rgrid, 'his', 'single');
    rho = rho(2:end-1,:,:);

    rback = dc_roms_read_data(runs.dir, 'rho', ...
                              1, {runs.bathy.axis ix+1 ix+1}, ...
                              [], runs.rgrid, 'his', 'single');
    % rho = rho - rback(2:end-1,:,:);

    mask = fillnan(bsxfun(@times, csdye < runs.csflux.x(isobath), ...
                          runs.csflux.offmask(:,tindex,isobath)),0);

    % profile I am assuming
    [videal, idmask] = runs.makeStreamerSection(isobath);

    tind = tindex; % FOR PARAMETERISATION
                   % syms x z;
    a = 2; % 2 for gaussian
           % V0 = runs.eddy.V(tind) * 2.33;
    R = runs.csflux.R;
    L = median(runs.eddy.rhovor.dia(1:tind))/2;
    % Lz = runs.eddy.Lgauss(tindex);
    % H = runs.csflux.h(isobath);
    yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
    y0oL =  R/L * (1 - yoR); % y0/L - used in derivation
    xfrac = -sqrt(1 - y0oL^a);

    if debug_flux
        % diagnosed flux
        flux = runs.csflux.off.slope(tindex,isobath,isobath)
        % should agree with above
        vtrue = trapz(xvec, trapz(zvec, repnan(csvel .* mask,0), 2), 1)

        % idealized velocity with real mask
        vitruemask = trapz(xvec, trapz(zvec, repnan(videal .* mask, 0), 2), 1)

        % idealized parameterization
        vest = trapz(xvec, trapz(zvec, repnan(videal .* idmask, 0), 2), 1)
    end

    figure;
    insertAnnotation([runs.name '.plot_xzsection']);
    maximize;

    ax(1) = subplot(221);
    [~,h1] = contourf(xvec/1000, zvec, csvel', 20);
    h1.EdgeColor = 'none';
    hold on
    contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    runs.add_timelabel(tindex);
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title('Cross-shelf velocity (m/s)');
    %linex(xfrac);
    liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
          {'vertical scale'; 'h_{sb}'});
    caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
    beautify;
    if exist('isobath', 'var') & (isobath == 1)
        ylim([min(zvec(:)) 0]);
    end

    xsb = runs.bathy.xsb;
    ax(2) = subplot(224);
    pcolorcen(xvec/1000, zvec, runs.sgntamp*(csdye'-xsb)/1000); % .* mask
    hold on; shading interp
    contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    hcb = colorbar;
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title('Cross-shelf dye (km)');
    runs.add_timelabel(tindex);
    linkaxes(ax, 'xy');
    %linex(xfrac);
    liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
          {'vertical scale'; 'h_{sb}'});

    % muck with colorbar
    cmin = round(min(runs.sgntamp*(csdye(:) - xsb)/1000));
    cmax = round(max(runs.sgntamp*(csdye(:) - xsb)/1000));
    hcb.TickLabelsMode = 'auto';
    hcb.TickDirection = 'out';
    hcb.Limits = [cmin cmax];
    hcb.Ticks = [cmin cmin/2 0 round((runs.bathy.xsl + xsb)/2000 - xsb/1000) ...
                 round(runs.params.eddy.cy/1000 - xsb/1000)];
    % hcb.TickLabels{1} = ['Shelfbreak - ' num2str(-1*hcb.Ticks(1)) ' km'];
    hcb.TickLabels{2} = 'Shelf Water';
    % hcb.TickLabels{3} = 'Shelfbreak';
    hcb.TickLabels{4} = 'Slope Water';
    hcb.TickLabels{end} = 'Eddy Water';
    beautify;

    ax(3) = subplot(223);
    contour(xvec/1000, zvec, rho', 30, 'k'); % .* mask
    clim = caxis;
    hold on; shading interp;
    contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    caxis(clim);
    colorbar;
    xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
    title('\rho (kg/m^3)');
    runs.add_timelabel(tindex);
    linkaxes(ax, 'xy');
    %linex(xfrac);
    beautify;

    if ~exist('isobath', 'var')
        ax(4) = subplot(224);
        pcolorcen(xvec/1000, zvec, eddye'); % .* mask
        hold on; shading interp;
        contour(xvec/1000, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
        colorbar;
        xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
        title(['Eddy dye | ' runs.name]);
        runs.add_timelabel(tindex);
        linkaxes(ax, 'xy');
        % linex(xfrac);
        liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
              {'vertical scale'; 'h_{sb}'});
        caxis([0 1]);
        beautify;
    else
        subplot(222);
        plot(abs(trapz(xvec, repnan(csvel.*mask,0), 1)), zvec);
        title(runs.name);
        ylabel('Z (m)');
        xlabel('|Transport| (m^2/s)');
        liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
              {'vertical scale'; 'h_{sb}'});
        ylim(ax(1).YLim);
        beautify;
    end

    if debug_flux
        figure;
        insertAnnotation([runs.name '.plot_xzsection']);
        maximize;

        ax(1) = subplot(2,3,[1 2]);
        [~,h1] = contourf(xvec/1000, zvec, csvel', 20);
        hold on
        contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
        runs.add_timelabel(tindex);
        xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
        title(['Cross-shelf velocity (m/s) | ' runs.name]);
        linex(xfrac*L/1000, 'xfrac');
        liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
              {'vertical scale'; 'h_{sb}'});
        caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
        beautify;

        ax(2) = subplot(2,3,[4 5]);
        [~,h1] = contourf(xvec/1000, zvec, videal', 20);
        hold on
        contour(xvec/1000, zvec, repnan(idmask',0), [1 1], 'k', 'LineWidth', 2);
        runs.add_timelabel(tindex);
        xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
        title('Idealized');
        linex(xfrac*L/1000, 'xfrac');
        liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
              {'vertical scale'; 'h_{sb}'});
        caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
        beautify;

        subplot(2,3,[3 6]);
        plot(trapz(xvec, repnan(csvel.*mask,0), 1), zvec);
        hold on;
        plot(trapz(xvec, videal.*idmask, 1), zvec);
        legend('Actual', 'Ideal', 'Location', 'SouthEast');
        beautify;

        linkaxes(ax, 'xy');
    end
end
