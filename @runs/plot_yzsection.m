% plot eddye - y-z cross-sections to compare against diagnosed vertical
% scale
function [] = plot_yzsection(runs, days, loc)

    plot_pbot = 0;
    if plot_pbot
        pb = load([runs.dir '/pbot.mat'], 'pbot', 'xvec', 'yvec');
    end

    [~,~,tind] = runs.locate_resistance;
    % hack for when I'm trying to provide non-dimensional times
    %if all(days < 1)
    %tindices = vecfind(runs.ndtime, days);
    %tindices = tind;
        %else
        %tindices = vecfind(runs.time/86400, days)
        %end

    tindices = runs.process_time(days);
    days = runs.eddy.t(tindices);

    nt = length(tindices);

    hf1 = figure; maximize();% - eddye
    %hf2 = figure; maximize();% - rho
    % hf3 = figure; maximize();% - zdye
    %hf4 = figure; maximize();% - u
    %hf5 = figure; maximize();% - v
    %hf6 = figure; maximize();% - csdye
    %zdback = double(squeeze(ncread(runs.out_file, runs.zdname, ...
    %                               [1 1 1 1], [1 Inf Inf
    %                               1])));
    %zback = runs.rgrid.z_r(:,:,1)';

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

    for ii=1:nt
        locstr = num2str(loc(ii));

        if plot_pbot
            iloc = find_approx(pb.xvec, str2double(loc),1);
        end

        ed = dc_roms_read_data(runs.dir, runs.eddname, tindices(ii), ...
                               {bathyax locstr locstr}, [], runs.rgrid, 'his');

        if exist('hf1', 'var')
            figure(hf1);
            ax1(ii) = subplot(1, nt, ii);
            axall = [axall ax1(ii)];
            contourf(yz, zmat, ed, [0.1:0.1:1], 'EdgeColor', 'none');
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            colormap(cbrewer('seq', 'Greys', 32));
            caxis([0 1]);
            title('Eddy dye');

            common(runs, hf1, yz, zmat, drho, ed, ii, days, locstr, tindices);
        end

        if exist('hf2', 'var')
            figure(hf2);
            temp = dc_roms_read_data(runs.dir, 'rho', tindices(ii), ...
                                     {bathyax locstr locstr}, [], runs.rgrid, 'his');

            ax2(ii) = subplot(1, nt, ii);
            axall = [axall ax2(ii)];
            drho = bsxfun(@minus, temp, tback);
            [cc,hh] = contourf(yz, zmat, temp, 40, 'EdgeColor', 'none');
            % colormap(flipud(cbrewer('seq','Blues',12)));
            %           hh.EdgeColor = 'none';

            shading flat; hold on;

            hcbar = colorbar;
            hcbar.Label.String = '$$\rho - \bar\rho(z)$$';
            hcbar.Label.Interpreter = 'latex';

            if ii == 1
                clim = caxis; %[-0.0553 -0.0021]; %caxis;
            end
            caxis(clim);

            if plot_pbot
                pbvec = pb.pbot(iloc,:,tindices(ii));
                pbvec = pbvec ./ max(abs(pbvec(:)));
                z0 = mean(ylim);
                plot(pb.yvec/1000, z0 + 100*pbvec)
                liney(z0, 'pbot anom = 0','k');
            end

            common(runs, hf2, yz, zmat, drho, ed, ii, days, locstr, tindices);
        end

        if exist('hf3', 'var')
            figure(hf3);
            zd = dc_roms_read_data(runs.dir, runs.zdname, tindices(ii), ...
                                  {bathyax locstr locstr}, [], runs.rgrid, 'his');

            ax3(ii) = subplot(1, nt, ii);
            axall = [axall ax3(ii)];
            contourf(yz, zmat, zd-zback);
            shading flat; hold on
            caxis( [-1 1] * max(abs(zd(:)-zback(:))) );
            center_colorbar;

            common(runs, hf3, yz, zmat, drho, ed, ii, days, locstr, tindices);
        end

        if exist('hf4', 'var')
            figure(hf4);
            u = dc_roms_read_data(runs.dir, 'u', tindices(ii), ...
                                  {bathyax locstr locstr}, [], ...
                                  runs.rgrid, 'his');

            ax4(ii) = subplot(1, nt, ii);
            axall = [axall ax4(ii)];
            if runs.bathy.axis == 'y'
                contourf(yz, zmat, u, 20, 'EdgeColor', 'none');
                hold on;
                contour(yz, zmat, u, [0 0], 'Color', 'k', 'LineWidth', 2);
            else
                contourf(yz(2:end-1,:), zmat(2:end-1,:), avg1(u,1), 20, ...
                         'EdgeColor', 'none');
                hold on;
                contour(yz(2:end-1,:), zmat(2:end-1,:), u, ...
                        [0 0], 'Color', 'k', 'LineWidth', 2);
            end
            shading flat;

            if plot_pbot
                pbvec = pb.pbot(iloc,:,tindices(ii));
                pbvec = pbvec ./ max(abs(pbvec(:)));
                z0 = mean(ylim);
                plot(pb.yvec/1000, z0 + 100*pbvec)
                liney(z0, 'pbot anom = 0','k');
            end

            caxis([-1 1] * max(abs(u(:)))); center_colorbar;
            common(runs, hf4, yz, zmat, drho, ed, ii, days, locstr, tindices);
        end

        if exist('hf5', 'var')
            figure(hf5);
            v = dc_roms_read_data(runs.dir, 'v', tindices(ii), ...
                                  {bathyax locstr locstr}, [], runs.rgrid, 'his');
            ax5(ii) = subplot(1, nt, ii);
            axall = [axall ax5(ii)];
            if runs.bathy.axis == 'y'
                contourf(yz(2:end-1,:), zmat(2:end-1,:), avg1(v,1), 20, ...
                         'EdgeColor', 'none');
                hold on
                contour(yz(2:end-1,:), zmat(2:end-1,:), avg1(v,1), ...
                         [0 0], 'LineWidth', 2, 'Color', 'k');
            else
                contourf(yz, zmat, v, 'EdgeColor', 'none');
                hold on
                contour(yz, zmat, v, [0 0], 'k', 'LineWidth', 2);
            end
            shading flat;
            caxis( [-1 1] * max(abs(v(:)))); center_colorbar;
            common(runs, hf5, yz, zmat, drho, ed, ii, days, locstr, tindices);
        end

        if exist('hf6', 'var')
            [csdye, ~, yax, zax, ~] = ...
                dc_roms_read_data(runs.dir, runs.csdname, tindices(ii), ...
                                  {bathyax locstr locstr}, [], runs.rgrid, 'his');

            figure(hf6);
            ax6(ii) = subplot(1, nt, ii);
            axall = [axall ax6(ii)];
            contourf(yax/1000, zax, csdye/1000, 40, 'EdgeColor', 'none');
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));

            common(runs, hf6, yz, zmat, drho, ed, ii, days, locstr, tindices);
        end

    end

    if exist('hf1')
        figure(hf1)
        % suplabel('eddy dye', 't');
        insertAnnotation([runs.name '.plot_yzsection']);
    end

    if exist('hf2', 'var')
        figure(hf2)
        %[~,ht] = suplabel(['\rho anomaly | (black, grey) contours = (eddye, \rho ' ...
        %                   'threshold)'], 't');
        insertAnnotation([runs.name '.plot_yzsection']);
        if ~isempty(findall(gcf, 'type', 'colorbar'))
            hcbar = findall(gcf,'type','colorbar');
            for ii = 1:length(hcbar)
                hcbar(ii).Label.Rotation = 90;
            end
        end
    end

    if exist('hf3', 'var')
        figure(hf3)
        suplabel('z-dye - z-level', 't');
        insertAnnotation([runs.name '.plot_yzsection']);
    end

    if exist('hf4', 'var')
        figure(hf4)
        suplabel('u - along-shore', 't');
        insertAnnotation([runs.name '.plot_yzsection']);
    end

    if exist('hf5', 'var')
        figure(hf5)
        suplabel('v - cross-shore', 't');
        insertAnnotation([runs.name '.plot_yzsection']);
    end

    if exist('hf6', 'var')
        figure(hf6)
        suplabel('Cross-shelf dye', 't');
        insertAnnotation([runs.name '.plot_yzsection']);
    end

    linkaxes(axall, 'xy');
end

function common(obj, hf, yz, zmat, drho, ed, ii, days, loc, tindices)
% do common tasks
    figure(hf); drawnow;
    clim = caxis; limx = xlim; limy = ylim;

    % colorbar on last panel only
    if ii == length(days)
        colorbar;
    else
        colorbar('hide');
    end

    if ~isempty(ed)
        contour(yz, zmat, ed, 1, 'r', 'LineWidth', 2);
    end
    if ~isempty(drho)
        contour(yz, zmat, drho, [1 1]* obj.eddy.drhothreshssh(1), ...
                'Color', [44 162 95]/256, 'LineWidth', 2);
    end
    caxis(clim);

    % vertical scale
    liney(-1 * obj.eddy.Lgauss(tindices(ii)), [], 'k');

    linex(obj.eddy.my(tindices(ii))/1000, [], 'k');

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
        end
    end

    text(textx*limx(2), -1 * obj.eddy.Lgauss(tindices(ii)), ...
         {'vertical','scale'}, 'VerticalAlignment', 'Bottom', ...
         'HorizontalAlignment','Center');

    figure(hf);
    text(0.05 , 0.05, {['t = ' num2str(days(ii)) ' days']; ...
                       ['x = ' num2str(str2double(loc)/1000, '%3.0f') ' km']}, ...
         'Units', 'normalized', 'Color', 'w');

    if ii == 1, ylabel('Z (m)'); end
    xlabel([upper(obj.bathy.axis) '(km)']);
    beautify([15 15 18]);
end