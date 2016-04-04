% plot eddye - y-z cross-sections to compare against diagnosed vertical
% scale
function [handles] = PlotSingleYZSection(runs, varname, days, loc, hax)

    if ~exist('hax', 'var')
        hf = figure; hax = gca;
    else
        axes(hax);
    end

    plot_pbot = 0;
    if plot_pbot
        pb = load([runs.dir '/pbot.mat'], 'pbot', 'xvec', 'yvec');
    end

    [~,~,tind] = runs.locate_resistance;

    tindices = runs.process_time(days)
    varname = runs.process_varname(varname);
    days = runs.eddy.t(tindices);

    drho = []; ed = []; axall = [];

    if runs.bathy.axis == 'y'
        cen = runs.eddy.mx;
        bathyax = 'x';
        yz = repmat(runs.rgrid.y_rho(:,1), [1 runs.rgrid.N]) / 1000;
        zmat = runs.rgrid.z_r(:,:,1)';
        tback = double(squeeze(ncread(runs.out_file, 'rho', [1 1 1 1], ...
                                      [1 Inf Inf 1])));
    else
        cen = runs.eddy.my;
        bathyax = 'y';
        yz = repmat(runs.rgrid.x_rho(1,:)', [1 runs.rgrid.N]) / ...
             1000;
        zmat = squeeze(runs.rgrid.z_r(:,1,:))';
        tback = double(squeeze(ncread(runs.out_file, 'rho', [1 1 1 1], ...
                                      [Inf 1 Inf 1])));
    end

    if ~exist('loc', 'var') || isempty(loc)
        loc = cen(tindices);
    end

    locstr = num2str(loc);

    if plot_pbot
        iloc = find_approx(pb.xvec, str2double(loc),1);
    end

    if strcmpi(varname, 'rhoanom')
        var = dc_roms_read_data(runs.dir, 'rho', tindices, ...
                                {bathyax locstr locstr}, [], runs.rgrid, 'his');
        var = bsxfun(@minus, var, tback);
        varname = '\rho''';
    else
        var = dc_roms_read_data(runs.dir, varname, tindices, ...
                                {bathyax locstr locstr}, [], runs.rgrid, 'his');
    end

    if strcmpi(varname, runs.csdname)
        var = var/1000;
    end

    handles.hvar = contourf(yz, zmat, var, 20, 'EdgeColor', 'none');
    handles.hlz = liney(-1 * runs.eddy.Lgauss(tindices));
    if strcmpi(varname, runs.eddname)
        colormap(cbrewer('seq', 'Greys', 32));
        caxis([0 1]);
    end
    handles.htitle = title([runs.name ' | ' varname]);
    handles = common(runs, yz, zmat, drho, var, days, locstr, tindices);
end

function handles = common(obj, yz, zmat, drho, ed, days, loc, tindices)
% do common tasks
    drawnow;
    clim = caxis; limx = xlim; limy = ylim;

    colorbar;

    if ~isempty(ed)
        [~,handles.hed] = contour(yz, zmat, ed, 1, 'r', 'LineWidth', 2);
    end
    if ~isempty(drho)
        [~,handles.hdrho] = contour(yz, zmat, drho, [1 1]* obj.eddy.drhothreshssh(1), ...
                                    'Color', [44 162 95]/256, 'LineWidth', 2);
    end
    caxis(clim);

    % vertical scale
    liney(-1 * obj.eddy.Lgauss(tindices), [], 'k');

    linex(obj.eddy.my(tindices)/1000, [], 'k');

    % patch bathymetry
    if obj.bathy.loc == 'h'
        % northern coast
        patch(([obj.rgrid.y_rho(:,1); max(obj.rgrid.y_rho(:,1))])./1000, ...
              -1*[obj.rgrid.h(:,1); max(obj.rgrid.h(:))], 'k');
        textx = 0.15;
    else
        if obj.bathy.axis == 'y'
            % southern coast
            patch(([obj.rgrid.y_rho(:,1); min(obj.rgrid.y_rho(:,1))])./1000, ...
                  -1*[min(obj.rgrid.h(:)); obj.rgrid.h(:,1)], 'k');
            textx = 0.85;
        else
            % western coast
            patch(([obj.rgrid.x_rho(1,:)'; min(obj.rgrid.x_rho(1,:))])./1000, ...
                  -1*[min(obj.rgrid.h(:)); obj.rgrid.h(1,:)'], 'k');
            textx = 0.85;
        end
    end

    text(textx*limx(2), -1 * obj.eddy.Lgauss(tindices), ...
         {'vertical','scale'}, 'VerticalAlignment', 'Bottom', ...
         'HorizontalAlignment','Center');

    handles.htlabel = ....
        text(0.05 , 0.05, {['t = ' num2str(days) ' days']; ...
                        ['x = ' num2str(str2double(loc)/1000, '%3.0f') ' km']}, ...
             'Units', 'normalized', 'Color', 'w');

    ylabel('Z (m)');
    xlabel([upper(obj.bathy.axis) '(km)']);
    beautify([15 15 18]);
end