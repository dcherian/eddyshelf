% plot eddye - y-z cross-sections to compare against diagnosed vertical
% scale
function [] = plot_yzsection(runs, days)

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

    tindices = days;
    days = runs.eddy.t(days);

    nt = length(tindices);

    hf1 = figure; maximize();% - eddye
    hf2 = figure; maximize();% - rho
    %hf3 = figure; maximize();% - zdye
    hf4 = figure; maximize();% - u
    hf5 = figure; maximize();% - v
    %zdback = double(squeeze(ncread(runs.out_file, runs.zdname, ...
    %                               [1 1 1 1], [1 Inf Inf
    %                               1])));
    %zback = runs.rgrid.z_r(:,:,1)';

    drho = []; ed = [];

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

    for ii=1:nt
        loc = num2str(cen(tindices(ii)))
        loc  = '395000';
        if plot_pbot
            iloc = find_approx(pb.xvec, str2double(loc),1);
        end

        ed = dc_roms_read_data(runs.dir, runs.eddname, tindices(ii), ...
                               {bathyax loc loc}, [], runs.rgrid, 'his');

        if exist('hf1', 'var')
            figure(hf1);
            ax1(ii) = subplot(1, nt, ii);
            contourf(yz, zmat, ed, [0.1:0.1:1]);
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            center_colorbar;

            common(runs, hf1, yz, zmat, drho, ed, ii, days, tindices);
        end

        if exist('hf2', 'var')
            figure(hf2);
            temp = dc_roms_read_data(runs.dir, 'rho', tindices(ii), ...
                                     {bathyax loc loc}, [], runs.rgrid, 'his');

            ax2(ii) = subplot(1, nt, ii);
            drho = bsxfun(@minus, temp, tback);
            [cc,hh] = contourf(yz, zmat, temp, 40);
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

            common(runs, hf2, yz, zmat, drho, ed, ii, days, tindices);
        end

        if exist('hf3', 'var')
            figure(hf3);
            zd = dc_roms_read_data(runs.dir, runs.zdname, tindices(ii), ...
                                  {bathyax loc loc}, [], runs.rgrid, 'his');

            ax3(ii) = subplot(1, nt, ii);
            contourf(yz, zmat, zd-zback);
            shading flat; hold on
            caxis( [-1 1] * max(abs(zd(:)-zback(:))) );
            center_colorbar;

            common(runs, hf3, yz, zmat, drho, ed, ii, days, tindices);
        end

        if exist('hf4', 'var')
            figure(hf4);
            u = dc_roms_read_data(runs.dir, 'u', tindices(ii), ...
                                  {bathyax loc loc}, [], ...
                                  runs.rgrid, 'his');

            ax4(ii) = subplot(1, nt, ii);
            if runs.bathy.axis == 'y'
                contourf(yz, zmat, u, 20);
                hold on;
                contour(yz, zmat, u, [0 0], 'Color', 'k', 'LineWidth', 2);
            else
                contourf(yz(2:end-1,:), zmat(2:end-1,:), avg1(u,1), 20);
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
            common(runs, hf4, yz, zmat, drho, ed, ii, days, tindices);
        end

        if exist('hf5', 'var')
            figure(hf5);
            v = dc_roms_read_data(runs.dir, 'v', tindices(ii), ...
                                  {bathyax loc loc}, [], runs.rgrid, 'his');
            ax5(ii) = subplot(1, nt, ii);
            if runs.bathy.axis == 'y'
                contourf(yz(2:end-1,:), zmat(2:end-1,:), avg1(v,1), 20);
                hold on
                contour(yz(2:end-1,:), zmat(2:end-1,:), avg1(v,1), ...
                         [0 0], 'LineWidth', 2, 'Color', 'k');
            else
                contourf(yz, zmat, v);
                hold on
                contour(yz, zmat, v, [0 0], 'k', 'LineWidth', 2);
            end
            shading flat;
            caxis( [-1 1] * max(abs(v(:)))); center_colorbar;
            common(runs, hf5, yz, zmat, drho, ed, ii, days, tindices);
        end
    end

    if exist('hf1')
        figure(hf1)
        suplabel('eddy dye', 't');
        spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
        linkaxes(ax1, 'xy');
        insertAnnotation([runs.name '.plot_eddye']);
    end

    if exist('hf2', 'var')
        figure(hf2)
        %[~,ht] = suplabel(['\rho anomaly | (black, grey) contours = (eddye, \rho ' ...
        %                   'threshold)'], 't');
        %set(ht, 'FontSize', 20);
        linkaxes(ax2, 'xy');
        insertAnnotation([runs.name '.plot_eddye']);
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
        spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
        linkaxes(ax3, 'xy');
        insertAnnotation([runs.name '.plot_eddye']);
    end

    if exist('hf4', 'var')
        figure(hf4)
        suplabel('u - along-shore', 't');
        %spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
        linkaxes(ax4, 'xy');
        insertAnnotation([runs.name '.plot_eddye']);
    end

    if exist('hf5', 'var')
        figure(hf5)
        suplabel('v - cross-shore', 't');
        spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
        linkaxes(ax5, 'xy');
        insertAnnotation([runs.name '.plot_eddye']);
    end
end

function common(obj, hf, yz, zmat, drho, ed, ii, days, tindices)
% do common tasks
    figure(hf); drawnow;
    colorbar; clim = caxis; limx = xlim; limy = ylim;
    if ~isempty(ed)
        contour(yz, zmat, ed, 1, 'r', 'LineWidth', 2);
    end
    if ~isempty(drho)
        contour(yz, zmat, drho, [1 1]* obj.eddy.drhothreshssh(1), ...
                'Color', [1 1 1]*0, 'LineWidth', 2);
    end
    caxis(clim);

    % vertical scale
    liney(-1 * obj.eddy.Lgauss(tindices(ii)), [], 'k');
    text(0.85*limx(2), -1 * obj.eddy.Lgauss(tindices(ii)), ...
         {'vertical','scale'}, 'VerticalAlignment', 'Bottom', ...
         'HorizontalAlignment','Center');

    % patch bathymetry
    patch(([obj.rgrid.y_rho(:,1); min(obj.rgrid.y_rho(:,1))])./1000, ...
          -1*[min(obj.rgrid.h(:)); obj.rgrid.h(:,1)]./1, 'k');

    figure(hf);
    text(0.05 , 0.05, ['t = ' num2str(days(ii)) ' days'], ...
         'Units', 'normalized', 'Color', 'w');
    title(['day' num2str(days(ii))]);
    if ii == 1, ylabel('Z (m)'); end
    xlabel([upper(obj.bathy.axis) '(km)']);
    beautify([15 15 18]);

end