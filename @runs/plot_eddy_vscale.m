% check eddy vertical scale estimations
function [] = plot_eddy_vscale(runs)

    c = hypot(runs.eddy.cvx, runs.eddy.cvy) / 86.4;
    c = nanmean(c(1:50));

    % read in initial velocity field
    u0 = ncread(runs.out_file, 'u', [1 1 1 1], [Inf Inf Inf ...
                        1]);
    v0 = ncread(runs.out_file, 'v', [1 1 1 1], [Inf Inf Inf ...
                        1]);
    U = hypot(avg1(u0(:,2:end-1,:), 1), avg1(v0(2:end-1,:,:), ...
                                             2));

    % volume of eddy that satisfies U/c criterion
    dV = bsxfun(@times, runs.rgrid.dV(2:end-1, 2:end-1,:) ...
                .* (U > c), runs.eddy.vor.mask(:,:,1));
    runs.eddy.Ucvol = nansum(dV(:));

    if ~isfield(runs.eddy, 'zT') || isempty(runs.eddy.zT)
        for ii=1:size(runs.eddy.T, 1)
            ix = vecfind(runs.rgrid.x_rho(1,:), ...
                         runs.eddy.vor.cx(ii));
            iy = vecfind(runs.rgrid.y_rho(:,1), ...
                         runs.eddy.vor.cy(ii));
            runs.eddy.zT(ii,:) = squeeze(runs.rgrid.z_r(:, iy, ix))';
        end
        runs.eddy.tmat = repmat(runs.time', [1 size(runs.eddy.T, ...
                                                    2)]);
        eddy = runs.eddy;
        save([runs.dir '/eddytrack.mat'], 'eddy');
    end

    figure;
    if isfield(runs.eddy, 'dyecen')
        subplot(211)
        pcolorcen(runs.eddy.tmat, runs.eddy.zT, ...
                  runs.eddy.dyecen);
        %   colormap(flipud(colormap('bone')));
        caxis([0 1]);
        colorbar;
        hold all
        plot(runs.time/86400, -1*runs.eddy.Lz2, 'c');
        plot(runs.time/86400, -1*runs.eddy.Lgauss, 'm');

        tcen = find_approx(runs.eddy.my, runs.bathy.xsl, ...
                           1);
        tse = find_approx(runs.eddy.se, runs.bathy.xsl, 1);
        linex([tse tcen], [], 'r');
        title(['Eddy dye profiles | ' runs.name]);
        subplot(212)
    end
    pcolorcen(runs.eddy.tmat, runs.eddy.zT, ...
              runs.eddy.T./max(runs.eddy.T(1,:)));
    %        colormap(flipud(colormap('bone')));
    caxis([0 1]);
    hold all
    plot(runs.time/86400, -1*runs.eddy.Lz2, 'c');
    plot(runs.time/86400, -1*runs.eddy.Lgauss, 'm');
    ylabel('Z (m)'); xlabel('Time (days)');
    title('Scaled temp anomaly at eddy center');
    colorbar;
    legend('Scaled temp anomaly', 'sine fit', 'Gaussian fit', 'Location', ...
           'SouthEast');
    contour(runs.eddy.tmat, runs.eddy.zT, runs.eddy.T, [0], ...
            'LineWidth', 2,'Color', 'k');
end
