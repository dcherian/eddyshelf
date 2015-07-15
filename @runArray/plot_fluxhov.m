% hovmoeller diagrams of cross-shelfbreak flux
function [] = plot_fluxhov(runArray, axis)

    if ~exist('axis', 'var') || isempty(axis)
        axis = 'z';
    end

    cmap = flipud(cbrewer('div','RdBu',24));

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);
        run = runArray.array(ii);

        figure;
        hold on
        insertAnnotation('runArray.plot_fluxhov')

        if axis == 'x'
            matrix = run.csflux.shelfxt(:,:,1);
            tmat = repmat(run.csflux.time/86400, [size(matrix, 1) 1]);
            xmat = repmat(run.rgrid.x_rho(1,2:end-1)'/1000, ...
                          [1 size(tmat, 2)]);
        else
            matrix = run.csflux.shelfzt(:,:,1);
            tmat = repmat(run.csflux.time/run.tscale, [size(matrix, 1) 1]);
            xmat = repmat(run.rgrid.z_r(:,run.bathy.isb,1), ...
                          [1 size(tmat, 2)]);
            xmat = xmat./max(abs(xmat(:)));
        end

        pcolorcen(xmat, tmat, matrix);
        colormap(cmap);
        colorbar; center_colorbar;
        title(run.name);

        if axis == 'x'
            plot(run.eddy.vor.cx/1000, run.eddy.t, 'k');
            plot(run.eddy.vor.ee/1000, run.eddy.t, 'k');
            plot(run.eddy.vor.we/1000, run.eddy.t, 'k');
            ylabel('Time (days)');
            xlabel('X (km)');
        else
            xlim([0 3]);
            ylim([-1 0]);
            ylabel('z/H_sb');
            xlabel('Time (days)');
        end
        cblabel('Normalized transport');
    end
end
