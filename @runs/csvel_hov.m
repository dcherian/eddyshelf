function [] = csvel_hov(runs, iz, loc)
    if ~exist('loc', 'var') || isempty(loc)
        loc = [-30 -20 -10 0] * 1000 + runs.bathy.xsb;
        locu = [200 250 300 350]*1000;
    end

    flag_v = 1; % csvel
    flag_z = 0; % zeta
    flag_u = 0; % csvel y-z

    if ~exist('iz', 'var') || isempty(iz)
        iz = runs.rgrid.N;
    end

    if iz == runs.rgrid.N
        if isempty(runs.vsurf)
            runs.read_velsurf;
        end
        u = runs.usurf;
        v = runs.vsurf;
        titlestr = 'surface';
    else
        if isempty(runs.vbot)
            runs.read_velbot;
        end
        u = runs.ubot;
        v = runs.vbot;
        titlestr = 'bottom';
    end

    % step topography phase speed estimate
    cp = runs.params.phys.f0 * runs.bathy.xsb * ...
         (runs.bathy.hsl - runs.bathy.hsb)./runs.bathy.hsl;

    if flag_v
        figure;
        for ii=1:length(loc)
            ax(ii) = subplot(2,2,ii);

            ind = find_approx(runs.rgrid.y_v(:,1), loc(ii), 1);

            xmat = repmat(runs.rgrid.x_v(1,:)', [1 size(runs.vsurf, 3)]);
            tmat = repmat(runs.time, [size(runs.vsurf, 1) 1]);
            pcolorcen(xmat/1000, tmat./runs.eddy.tscale, squeeze(runs.vsurf(:,ind,:)));
            colorbar; center_colorbar;
            xlabel('X (km)');
            ylabel('Time (non-dimensional)');
            title([titlestr ' v (m/s) at y = ' num2str(loc(ii)) ' km']);

            hold on
            plot(runs.eddy.vor.cx/1000, runs.eddy.t*86400 / runs.eddy.tscale);
            plot(runs.eddy.vor.ee/1000, runs.eddy.t*86400 / runs.eddy.tscale);
            plot(runs.eddy.vor.we/1000, runs.eddy.t*86400 / runs.eddy.tscale);

            limx = xlim;
            xvec = limx(1):10:limx(2);
            tvec = 0.5 + 1./0.02 .* (xvec * 1000)./runs.eddy.tscale;
            plot(xvec, tvec);

            beautify;
        end

        suplabel(runs.name, 't');
        linkaxes(ax, 'xy');
    end

    %%%%%%%%%%%%%% ζ
    % along-shelfbreak pressure gradient at shelfbreak
    if flag_z
        %zx = bsxfun(@rdivide, squeeze(diff(runs.zeta(:, runs.bathy.isb, :), ...
        %                          1, 1)), diff(runs.rgrid.x_rho(1,:)', ...
        %1, 1));
        loc2 = [-20 0 10 20] + runs.bathy.isb;
        runs.read_zeta;
        figure;
        for ii=1:length(loc2)
            ax(ii) = subplot(2,2,ii);
            pcolorcen(xmat/1000, tmat/runs.eddy.tscale, ...
                      squeeze(runs.zeta(:, loc2(ii), :)));
            colorbar; center_colorbar;
            hold on
            plot(runs.eddy.vor.cx/1000, runs.eddy.t*86400 / runs.eddy.tscale);
            plot(runs.eddy.vor.ee/1000, runs.eddy.t*86400 / runs.eddy.tscale);
            plot(runs.eddy.vor.we/1000, runs.eddy.t*86400 / ...
                 runs.eddy.tscale);
            xlabel('X (km)');
            ylabel('Time (non-dimensional)');
            title(['zeta (m) at y = ' ...
                   num2str(runs.rgrid.y_rho(loc2(ii),1))]);
            beautify;
        end

        %suplabel(runs.name, 't');
        linkaxes(ax, 'xy');
    end

    %%%%%%%%%%%%%% u
    if flag_u
        figure;
        for ii=1:length(locu)
            ax(ii) = subplot(2,2,ii);

            ind = find_approx(runs.rgrid.x_v(1,:), locu(ii), 1);

            ymat = repmat(runs.rgrid.y_u(:,1), [1 size(runs.usurf, 3)]);
            tmat = repmat(runs.time, [size(runs.usurf, 2) 1]);
            contourf(tmat./86400, ymat./1000, squeeze(runs.usurf(ind,:,:)));
            colorbar; center_colorbar;
            ylabel('Y (km)');
            xlabel('Time (non-dimensional)');
            ylim([0 runs.bathy.xsb/1000]);
            title([titlestr ' u (m/s) at x = ' num2str(locu(ii)/1000) ' km']);

            beautify;
        end
    end
end
