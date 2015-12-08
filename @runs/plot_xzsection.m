function [handles] = plot_xzsection(runs, loc, day, debug_flux)

    dxi = 5; dzi = 3;

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
                [~,tindex] = runs.calc_maxflux( ...
                    runs.recalculateFlux(2*runs.bathy.hsb, isobath, isobath));
            end
        else
            if runs.bathy.axis == 'y'
                loc = runs.rgrid.y_rho(loc,1);
            else
                loc = runs.rgrid.x_rho(1,loc);
            end
        end
    end

    if runs.bathy.axis == 'y'
        sx1 = runs.spng.sx1;
        sx2 = runs.spng.sx2;
        ix = find_approx(runs.rgrid.y_rho(:,1), loc);
        bathyax = 2; asax = 1;
        xvec = runs.rgrid.x_rho(1,2:end-1) - runs.eddy.mx(tindex);
        zvec = runs.rgrid.z_r(:, ix, 1);
    else
        sx1 = runs.spng.sy1;
        sx2 = runs.spng.sy2;
        ix = find_approx(runs.rgrid.x_rho(1,:), loc);
        bathyax = 1; asax = 2;
        xvec = (runs.rgrid.y_rho(2:end-1,1) - runs.eddy.my(tindex))';
        zvec = runs.rgrid.z_r(:, 1, ix);
    end

    xvec = xvec(sx1:sx2);
    L = runs.eddy.rhovor.dia(tindex)/2;
    xsb = runs.bathy.xsb;

    csvel = squeeze(avg1( ...
        dc_roms_read_data(runs.dir, runs.csvelname, tindex, ...
                          {runs.bathy.axis ix-1 ix}, [], runs.rgrid, ...
                          'his', 'single'), bathyax));
    csvel = csvel(sx1:sx2,:);

    %asvel = squeeze(avg1( ...
    %    dc_roms_read_data(runs.dir, runs.asvelname, tindex, ...
    %                      {runs.bathy.axis ix ix}, [], runs.rgrid, ...
    %                      'his', 'single'), asax));
    %asvel = asvel(sx1:sx2,:);

    w = squeeze(avg1( ...
        dc_roms_read_data(runs.dir, 'w', tindex, ...
                          {runs.bathy.axis ix ix}, [], runs.rgrid, ...
                          'his', 'single'), 2));
    w = w(sx1:sx2,:);

    % process cross-shelf dye
    csdye = dc_roms_read_data(runs.dir, runs.csdname, ...
                              tindex, {runs.bathy.axis ix+1 ix+1}, ...
                              [], runs.rgrid, 'his', 'single');
    csdye = csdye(sx1:sx2,:);

    eddye = dc_roms_read_data(runs.dir, runs.eddname, ...
                              tindex, {runs.bathy.axis ix+1 ix+1}, ...
                              [], runs.rgrid, 'his', 'single');
    eddye = eddye(sx1:sx2,:);

    rho = dc_roms_read_data(runs.dir, 'rho', ...
                            tindex, {runs.bathy.axis ix+1 ix+1}, ...
                            [], runs.rgrid, 'his', 'single');
    rho = rho(sx1:sx2,:);

    rback = dc_roms_read_data(runs.dir, 'rho', ...
                              1, {runs.bathy.axis ix+1 ix+1}, ...
                              [], runs.rgrid, 'his', 'single');
    % rho = rho - rback(sx1:sx2,:);

    if isobath == 1
        mask = ones(size(csvel));
    else
        mask = fillnan(bsxfun(@times, csdye < runs.csflux.x(isobath), ...
                              runs.csflux.offmask(sx1:sx2,tindex,isobath)),0);
    end

    % profile I am assuming
    [videal, idmask, xmask, zmask] = runs.makeStreamerSection(isobath);
    videal = videal(sx1:sx2,:);
    idmask = idmask(sx1:sx2,:);
    xmask = xmask(sx1:sx2);

    % tind = tindex; % FOR PARAMETERISATION
    %                % syms x z;
    % a = 2; % 2 for gaussian
    %        % V0 = runs.eddy.V(tind) * 2.33;
    % R = runs.csflux.R;
    % L = median(runs.eddy.rhovor.dia(1:tind))/2;
    % % Lz = runs.eddy.Lgauss(tindex);
    % % H = runs.csflux.h(isobath);
    % yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
    % y0oL =  R/L * (1 - yoR); % y0/L - used in derivation
    % xfrac = -sqrt(1 - y0oL^a);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.hax(1) = subplot(221);
    [~,handles.hcsvel] = contourf(xvec/1000, zvec, csvel', 20);
    handles.hcsvel.EdgeColor = 'none';
    hold on
    [~,handles.hmask(1)] = ...
        contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    handles.htitle(1) = title('Cross-shelf velocity (m/s)');
    %linex(xfrac);
    common(runs, tindex);
    caxis([-1 1] * max(abs(csvel(:))));
    handles.hcb(1) = center_colorbar;
    handles.htime = runs.add_timelabel(tindex);
    handles.htime.Position = [0.68 0.13];
    if exist('isobath', 'var')
        handles.htime.String = [{handles.htime.String;  ...
                            ['y/R = ' ...
                            num2str(runs.csflux.ndloc(isobath), '%.2f')]}];
    end
    beautify;
    if exist('isobath', 'var') & (isobath == 1)
        ylim([min(zvec(:)) 0]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.hax(2) = subplot(224);
    handles.hcsdye = pcolorcen(xvec/1000, zvec, runs.sgntamp*(csdye'-xsb)/1000); % .* mask
    hold on; shading interp
    [~,handles.hmask(4)] = ...
        contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    [~,hrho(1)] = contour(xvec/1000, zvec, rho', 30, 'k'); % .* mask
    handles.hcb(4) = colorbar;
    title('Cross-shelf dye (km)');
    linkaxes(handles.hax, 'xy');
    common(runs, tindex);

    % muck with colorbar
    cmin = round(min(runs.sgntamp*(csdye(:) - xsb)/1000));
    cmax = round(max(runs.sgntamp*(csdye(:) - xsb)/1000));
    handles.hcb(4).TickLabelsMode = 'auto';
    handles.hcb(4).TickDirection = 'out';
    handles.hcb(4).Limits = [cmin cmax];
    handles.hcb(4).Ticks = sort([cmin/2 cmin 0 ...
                        round((runs.bathy.xsl + xsb)/2000 - xsb/1000) ...
                        round(runs.params.eddy.cy/1000 - xsb/1000)]);
    % handles.hcb(4).TickLabels{1} = ['Shelfbreak - ' num2str(-1*handles.hcb(4).Ticks(1)) ' km'];
    handles.hcb(4).TickLabels{2} = 'Shelf Water';
    % handles.hcb(4).TickLabels{3} = 'Shelfbreak';
    handles.hcb(4).TickLabels{4} = 'Slope Water';
    handles.hcb(4).TickLabels{end} = 'Eddy Water';
    beautify;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.hax(3) = subplot(223);
    [~,hw] = contourf(xvec/1000, zvec, w', 6);
    hw.EdgeColor = 'none';
    clim = caxis;
    hold on;
    [~,hrho(2)] = contour(xvec/1000, zvec, rho', 30, 'k'); % .* mask
    [~,handles.hmask(3)] = ...
        contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    %handles.hquiv = quiver(xvec(1:dxi:end)/1000, zvec(1:dzi:end), ...
    %                       asvel(1:dxi:end, 1:dzi:end)'/1000, ...
    %                       w(1:dxi:end, 1:dzi:end)');
    caxis(clim); center_colorbar;
    common(runs, tindex);
    title('\rho (kg/m^3)');
    linkaxes(handles.hax, 'xy');
    beautify;

    hrho(2).LevelList = hrho(1).LevelList;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('isobath', 'var')
        handles.hax(4) = subplot(224);
        pcolorcen(xvec/1000, zvec, eddye'); % .* mask
        hold on; shading interp;
        contour(xvec/1000, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
        colorbar;
        title(['Eddy dye | ' runs.name]);
        linkaxes(handles.hax, 'xy');
        common(runs, tindex);
        caxis([0 1]);
        beautify;
    else
        handles.hax(2) = subplot(222);
        plot(abs(trapz(xvec, repnan(csvel.*mask,0), 1)), zvec);
        title(runs.name);
        common(runs, tindex);
        xlabel('|Transport| (m^2/s)');
        ylim(handles.hax(1).YLim);
        pbaspect([1.3 1 1]);
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
        common(runs, tindex);
        title(['Cross-shelf velocity (m/s) | ' runs.name]);
        % linex(xfrac*L/1000, 'xfrac');
        caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
        beautify;

        ax(2) = subplot(2,3,[4 5]);
        [~,h1] = contourf(xmask/1000, zmask, videal', 20);
        hold on
        contour(xmask/1000, zmask, repnan(idmask',0), [1 1], 'k', 'LineWidth', 2);
        runs.add_timelabel(tindex);
        common(runs, tindex);
        title('Idealized');
        % linex(xfrac*L/1000, 'xfrac');
        caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
        beautify;

        subplot(2,3,[3 6]);
        plot(trapz(xvec, repnan(csvel.*mask,0), 1), zvec);
        hold on;
        plot(trapz(xmask, videal.*idmask, 1), zmask);
        legend('Actual', 'Ideal', 'Location', 'SouthEast');
        beautify;

        linkaxes(ax, 'xy');
    end

    handles.hrho = hrho;
end

function [] = common(obj, tindex)
    xlabel('X - X_{eddy} (km)'); ylabel('Z (m)');
    liney(-1 * [obj.eddy.Lgauss(tindex) obj.bathy.hsb], ...
          {'vertical scale'; 'h_{sb}'});
    beautify;
end