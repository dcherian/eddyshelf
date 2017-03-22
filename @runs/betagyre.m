function [] = betagyre(runs, depth)

    if ~exist('depth', 'var'), depth = 0; end

    use_cartesian = 1;
    use_psi = 1;

    runs.eddy.betagyre = nan(size(runs.eddy.mask));

    if ~use_psi
        runs.read_zeta;
        titlestr = 'SSH';
    else
        titlestr = '\psi';
    end
    if use_cartesian
        titlestr = [titlestr ' | cartesian'];
    else
        titlestr = [titlestr ' | polar'];
    end

    tic;
    for tt=1:1:size(runs.eddy.betagyre,3)
        if use_psi
            [var, iynan] = runs.streamfunction(tt, depth);
        else
            iynan = 1;
            var = runs.zeta(2:end-1,2:end-1,tt);
        end

        % var = var - nanmin(nanmin(var, [], 1), [], 2);

        xr = runs.rgrid.x_rho(1, 2:end-1)';
        yr = runs.rgrid.y_rho(2:end-1, 1)';
        ix = find_approx(xr, runs.eddy.mx(tt));
        iy = find_approx(yr(iynan:end), runs.eddy.my(tt));
        xvec = xr - runs.eddy.mx(tt);
        yvec = yr(iynan:end) - runs.eddy.my(tt);

        [x0, y0, v0x(tt), v0y(tt), exitx, exity] = ...
            FindCenter(xvec, yvec, var(:,iynan:end), ix, iy);

        if exitx == 0 | exity == 0
            disp('fit didn''t work. skipping')
            continue;
        end

        xref = runs.rgrid.x_rho(2:end-1, 2:end-1)' - runs.eddy.mx(tt) - x0;
        yref = runs.rgrid.y_rho(2:end-1, 2:end-1)' - runs.eddy.my(tt) - y0;

        if use_cartesian
            % Use LaCasce (1996) method of scaling down the symmetric
            % initial condition, relocating it to the center
            % (detected using fits), and then subtracting that from
            % the total signal.
            if tt == 1
                xref0 = xref;
                yref0 = yref;
                var0 = var;
            end

            vamp(tt) = mean([v0x(tt) v0y(tt)]);
            vmax(tt) = nanmax(var(:));

            var0int = interp2(xref0', yref0', var0', ...
                              xref, yref, 'spline');
            BetaGyre = var - var0int * vamp(tt)/vamp(1);

            if tt == 1
                BetaGyre0 = BetaGyre;
            end

        else
            % Determine radially-invariant signal
            % Doesn't work so well because of scatteredInterpolant.
            % Not sure how to do better.
            [thmat,rmat] = cart2pol(xref, yref);

            % zetaback = Interp2Cart(thmat, rmat, zeta, xref, yref);
            % zetaback and zeta agree really really well

            F = scatteredInterpolant(thmat(:), rmat(:), var(:));
            F.Method = 'natural';
            F.ExtrapolationMethod = 'none';
            [thi, ri] = ndgrid(linspace(min(thmat(:)), max(thmat(:)), 300), ...
                               linspace(0, max(rmat(:)), 300));
            varpol = F(thi, ri);
            % zetagrid = griddata(thmat(:), rmat(:), zeta(:), thi, ri, 'cubic');
            % zetagridback = Interp2Cart(thi, ri, zetagrid, xref, yref);
            % zetapolback and zeta have discrepancies.
            % zetapolback = Interp2Cart(thi, ri, zetapol, xref, yref);
            % dzeta = zetapolback - runs.zeta(:,:,tt);
            varsympol = repmat(nanmean(varpol, 1), [size(thi, 1) 1]);
            varsym = Interp2Cart(thi, ri, varsympol, xref, yref);
            varsym(varsym < 1e-3 * nanmax(abs(varsym(:)))) = NaN;
            % varsym(isnan(varsym)) = 0;

            BetaGyre = var - varsym;
        end

        runs.eddy.betagyre(:,:,tt) = BetaGyre;

        clf;
        pcolorcen(xref/1000, yref/1000, BetaGyre);
        center_colorbar; shading interp;
        L = 120;
        xlim([-1 1]*L); ylim([-1 1]*L);
        hold on;
        % contour(xref/1000, yref/1000, ...
        % runs.zeta(2:end-1, 2:end-1, tt), 20, 'k');
        contour(xref/1000, yref/1000, var, 20, 'k');
        plot(runs.eddy.mx/1000 - runs.eddy.mx(tt)/1000 - x0/1000, ...
             runs.eddy.my/1000 - runs.eddy.my(tt)/1000 - y0/1000, ...
             'b-');
        axis square;
        betatimescale = 1/(runs.params.phys.beta*runs.params.eddy.dia/2)/86400;
        title([titlestr ' | t = ' num2str(runs.eddy.t(tt)) ' days' ...
               ' = ' num2str(runs.eddy.t(tt)/betatimescale, '%.2f') ' 1/(\beta L)']);
        keyboard;
    end
    toc;
end

function [x0, y0, v0x, v0y, exitx, exity] = FindCenter(xvec, yvec, var, ix, iy)

    vxvec = double(var(:,iy))';
    vyvec = double(var(ix,:));

    % limit the region so that the fit works.
    dxy = 100;
    rangex = max(1, ix-dxy):min(ix+dxy, length(vxvec));
    rangey = max(1, iy-dxy):min(iy+dxy, length(vyvec));

    % fit used to find "actual" center that is not on a grid-point
    [v0x, ~, x0, ~, exitx] = ...
        gauss_fit(xvec(rangex),  vxvec(rangex), 0);
    [v0y, ~, y0, ~, exity] = ...
        gauss_fit(yvec(rangey), vyvec(rangey), 0);
end

function [out] = Interp2Cart(theta, r, in, x, y)
    xin = r.*cos(theta);
    yin = r.*sin(theta);
    [~,I,~] = unique([xin(:) yin(:)], 'first', 'rows');
    F = scatteredInterpolant(xin(I), yin(I), in(I));
    F.Method = 'natural';
    F.ExtrapolationMethod = 'none';
    out = F(x, y);
end