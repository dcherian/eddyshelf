% Returns [supply jet width, error on supplpy jet width, eddy penetration width]
function [supply, errsupp, eddyonshelf] = SupplyJetEddyonShelf(runs)

    debug = 0;
    xsb = runs.bathy.xsb;

    %% eddy penetration on shelf
    dz = diff(runs.supply.zmat, 1, 2);

    if debug
        figure;
        contourf(runs.supply.ymat/1000, runs.supply.zmat, ...
                 runs.supply.csdmean/1000, 40, 'EdgeColor', 'none');
        shading flat;
        hold on;
        contour(runs.supply.ymat/1000, runs.supply.zmat, ...
                runs.supply.csdmean/1000, ...
                [1 1]*runs.bathy.xsb/1000, 'k', 'LineWidth', 2);
        contour(runs.supply.ymat/1000, runs.supply.zmat, ...
                runs.supply.eddmean, ...
                [1 1]*0.2, 'g', 'LineWidth', 2);
    end

    csdint = sum(dz .* avg1(runs.supply.csdmean,2), 2)./sum(dz,2);
    eddint = sum(dz .* avg1(runs.supply.eddmean,2), 2)./sum(dz,2);
    eddcumvol = cumtrapz(runs.supply.ymat(:,1), ...
                         sum(dz .* avg1(runs.supply.eddmean,2), 2));
    eddcumvol = eddcumvol./max(eddcumvol);

    % different ways of finding dye front
    if debug
        N = 120;
        [yi, zi] = meshgrid(runs.supply.ymat(:, 1), ...
                            linspace(double(min(runs.supply.zmat(:))), 0, N));

        F = scatteredInterpolant(runs.supply.ymat(:), ...
                                 double(runs.supply.zmat(:)), ...
                                 runs.supply.csdmean(:), 'linear', 'none');
        csdi = F(yi, zi);
        ind = find(csdi(ceil(1/2*N),:) <= (runs.bathy.xsb), 1, 'last');
        liney(zi(ceil(1/2*N), 1));
        linex(runs.supply.ymat(ind, 1)/1000);
    end

    ind = find(csdint >= (runs.bathy.xsb), 1, 'first');
    %ind = find(runs.supply.eddmean(:,end) >= 0.7, 1, 'first');
    %ind = find(eddint >= 0.3, 1, 'first');
    %ind = find(eddcumvol >= 0.15, 1, 'first');
    eddyonshelf = (xsb - runs.supply.ymat(ind, 1))/1000;

    %% shelf water envelope

    if runs.bathy.axis == 'y'
        [~,stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1), 0.01);
        start = 1;
    else
        [~,stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1), 0.3);
        start = 1;
    end

    if strcmpi(runs.name, 'ew-8392')
        % biased by decrease near the end
        % stop = stop - 30;
    elseif strcmpi(runs.name, 'ew-8384') ...
            | strcmpi(runs.name, 'ew-8385') ...
            | strcmpi(runs.name, 'ew-8381') ...
            | strcmpi(runs.name, 'ew-8341')
        stop = length(runs.csflux.time);
    end

    env = xsb/1000 - runs.csflux.off.slopewater.envelope(start:stop,1)/1000;
    time = runs.csflux.time(start:stop)/86400;
    tvec = time(~isnan(env))/86400;
    env = env(~isnan(env));
    menv = mean(env);

    [y0,T,t0,y1,conf,fitobj] = tanh_fit(tvec, env, debug);
    if debug, title(runs.name); end
    supply = y0 + y1;

    % env time series is weird.
    % I think the fits are underestimating
    if strcmpi(runs.name, 'ew-8341') ...
            | strcmpi(runs.name, 'ew-8351-2')
        supply = supply + 3;
    end
    if strcmpi(runs.name, 'ew-8041')
         supply = supply + 2;
    end
    errsupp = hypot(conf(1,1)-y0, conf(1,4)-y1);

    % use transport binned by on-shelf origin
    % bins = avg1(run.csflux.off.bins{1})/1000;
    % %bins = bins - max(bins);
    % trans = run.csflux.off.slopewater.itrans{1}/1e10;
    % [y0, X, x0, y1, ~, conf] = gauss_fit(bins, trans, 0);
    % supply = abs(X);
    % errsupp = abs(abs(conf(1,2)) - abs(X));

end