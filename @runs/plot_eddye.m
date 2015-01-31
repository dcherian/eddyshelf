
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

    hf1 = figure; maximize();% - eddye
    hf2 = figure; maximize();% - rho
    %hf3 = figure; maximize();
    %hf4 = figure; maximize();
    %hf5 = figure; maximize();

    tback = double(squeeze(ncread(runs.out_file, 'rho', [1 1 1 1], ...
                                  [1 Inf Inf 1])));
    %zdback = double(squeeze(ncread(runs.out_file, runs.zdname, ...
    %                               [1 1 1 1], [1 Inf Inf
    %                               1])));
    zback = runs.rgrid.z_r(:,:,1)';
    for ii=1:nt
        figure(hf1);
        ed = dc_roms_read_data(runs.dir, runs.eddname, tindices(ii), ...
                               {'x' num2str(runs.eddy.cx(tindices(ii))) ...
                            num2str(runs.eddy.cx(tindices(ii)))}, [], ...
                               runs.rgrid, 'his');

        ax1(ii) = subplot(1, nt, ii);
        contour(yz, runs.rgrid.z_r(:,:,1)', ed, [0.1:0.1:1]);
        liney(-1 * runs.eddy.Lgauss(tindices(ii)));
        colorbar;
        colormap(flipud(colormap('bone')))
        title(['day' num2str(days(ii))]);
        beautify;

        figure(hf2);
        temp = dc_roms_read_data(runs.dir, 'rho', tindices(ii), ...
                                 {'x' num2str(runs.eddy.cx(tindices(ii))) ...
                            num2str(runs.eddy.cx(tindices(ii)))}, [], ...
                                 runs.rgrid, 'his');

        ax2(ii) = subplot(1, nt, ii);
        contourf(yz, runs.rgrid.z_r(:,:,1)', temp, 20);
        shading flat;
        liney(-1 * runs.eddy.Lgauss(tindices(ii)));
        colorbar;
        clim = caxis;
        %caxis([-0.05 0.05]); % [-1 1] * max(abs(temp(:))) );
        hold on;
        contour(yz, runs.rgrid.z_r(:,:,1)', ed, 1, 'k', ...
                'LineWidth', 2);
        caxis(clim);
        title(['day' num2str(days(ii))]);
        axis square
        beautify;

        %{
        figure(hf3);
        zd = dc_roms_read_data(runs.dir, runs.zdname, tindices(ii), ...
        {'x' num2str(runs.eddy.cx(tindices(ii))) ...
        num2str(runs.eddy.cx(tindices(ii)))}, [], ...
        runs.rgrid, 'avg');

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

        figure(hf5);
        v = dc_roms_read_data(runs.dir, 'v', tindices(ii), ...
                              {'x' num2str(runs.eddy.cx(tindices(ii))) ...
                            num2str(runs.eddy.cx(tindices(ii)))}, [], ...
                              runs.rgrid, 'avg');

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
        %}
    end

    figure(hf1)
    suplabel('eddy dye', 't');
    spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
    linkaxes(ax1, 'xy');
    insertAnnotation([runs.name '.plot_eddye']);

    figure(hf2)
    suplabel('temp anomaly', 't');
    spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));

    %figure(hf3)
    %suplabel('z-dye - z-level', 't');
    %spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));

    spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
    linkaxes(ax2, 'xy');
    insertAnnotation([runs.name '.plot_eddye']);

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
