
% plot eddye - y-z cross-sections to compare against diagnosed vertical
% scale
function [] = plot_eddye(runs, days)

    % hack for when I'm trying to provide non-dimensional times
    if all(days < 10)
        tindices = vecfind(runs.time./runs.tscale, days);
    else
        tindices = vecfind(runs.time/86400, days)
    end

    nt = length(tindices);
    yz = repmat(runs.rgrid.y_rho(:,1), [1 runs.rgrid.N]) / 1000;

    %hf1 = figure; maximize();% - eddye
    hf2 = figure; maximize();% - rho
    %hf3 = figure; maximize(); - zdye
    %hf4 = figure; maximize(); - u
    %hf5 = figure; maximize(); - v

    tback = double(squeeze(ncread(runs.out_file, 'rho', [1 1 1 1], ...
                                  [1 Inf Inf 1])));
    %zdback = double(squeeze(ncread(runs.out_file, runs.zdname, ...
    %                               [1 1 1 1], [1 Inf Inf
    %                               1])));
    zback = runs.rgrid.z_r(:,:,1)';
    for ii=1:nt

        loc = num2str(runs.eddy.mx(tindices(ii)));

        ed = dc_roms_read_data(runs.dir, runs.eddname, tindices(ii), ...
                               {'x' loc loc}, [], runs.rgrid, 'his');

        if exist('hf1', 'var')
            figure(hf1);
            ax1(ii) = subplot(1, nt, ii);
            contour(yz, runs.rgrid.z_r(:,:,1)', ed, [0.1:0.1:1]);
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            colorbar; caxis([0 1]);
            colormap(flipud(colormap('bone')))
            title(['day' num2str(days(ii))]);
            beautify;
        end

        if exist('hf2', 'var')
            figure(hf2);
            temp = dc_roms_read_data(runs.dir, 'rho', tindices(ii), ...
                                     {'x' loc loc}, [], runs.rgrid, 'his');

            ax2(ii) = subplot(1, nt, ii);
            drho = bsxfun(@minus, temp, tback);
            contourf(yz, runs.rgrid.z_r(:,:,1)', drho, 20);
            shading flat;
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            colorbar;
            if ii == 1
                clim = caxis; %[-0.0553 -0.0021]; %caxis;
                ylabel('Z (m)');
            end
            caxis(clim);
            hold on;
            contour(yz, runs.rgrid.z_r(:,:,1)', ed, 1, 'k', ...
                    'LineWidth', 2);
            contour(yz, runs.rgrid.z_r(:,:,1)', drho, ...
                    [runs.eddy.drhothresh(1) runs.eddy.drhothreshssh(1)], ...
                    'Color', [1 1 1]*0.3, 'LineWidth', 2);
            caxis(clim);
            title(['day ' num2str(days(ii))]);
            xlabel('Y (km)');
            %axis square
            beautify([15 15 18]);
        end

        if exist('hf3', 'var')
            figure(hf3);
            zd = dc_roms_read_data(runs.dir, runs.zdname, tindices(ii), ...
                                  {'x' loc loc}, [], runs.rgrid, 'his');

            ax3(ii) = subplot(1, nt, ii);
            contourf(yz/1000, runs.rgrid.z_r(:,:,1)', zd-zback);
            shading flat;
            hold on
            contour(yz/1000, runs.rgrid.z_r(:,:,1)', ed, 1, 'k', ...
                    'LineWidth', 2);
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            colorbar;
            caxis( [-1 1] * max(abs(zd(:)-zback(:))) );
            title(['day' num2str(days(ii))]);
            beautify;
        end

        if exist('hf4', 'var')
            figure(hf4);
            u = dc_roms_read_data(runs.dir, 'u', tindices(ii), ...
                                  {'x' num2str(runs.eddy.cx(tindices(ii))) ...
                                num2str(runs.eddy.cx(tindices(ii)))}, [], ...
                                  runs.rgrid, 'avg');

            ax4(ii) = subplot(1, nt, ii);
            contourf(yz/1000, runs.rgrid.z_r(:,:,1)', u);
            shading flat;
            hold on
            contour(yz/1000, runs.rgrid.z_r(:,:,1)', ed, 1, 'k', ...
                    'LineWidth', 2);
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            colorbar;
            caxis( [-1 1] * max(abs(u(:))));
            title(['day' num2str(days(ii))]);
        end

        if exist('hf5', 'var')
            figure(hf5);
            v = dc_roms_read_data(runs.dir, 'v', tindices(ii), ...
                                  {'x' loc loc}, [], runs.rgrid, 'his');
            ax5(ii) = subplot(1, nt, ii);
            contourf(yz(2:end-1,:)/1000, runs.rgrid.z_r(:,2:end-1,1)', avg1(v,1));
            shading flat;
            hold on
            contour(yz/1000, runs.rgrid.z_r(:,:,1)', ed, 1, 'k', ...
                    'LineWidth', 2);
            liney(-1 * runs.eddy.Lgauss(tindices(ii)));
            colorbar;
            caxis( [-1 1] * max(abs(v(:))));
            title(['day' num2str(days(ii))]);
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
        %suplabel('Density anomaly', 't');
        [~,ht] = suplabel(['\rho anomaly | (black, grey) contours = (eddye, \rho ' ...
                           'threshold) | grey line = Gaussian fit vertical scale'], ...
                          't');
        set(ht, 'FontSize', 20);
        %spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
        linkaxes(ax2, 'xy');
        insertAnnotation([runs.name '.plot_eddye']);
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
        spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
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
