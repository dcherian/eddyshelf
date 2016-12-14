function [] = balance_torques(runs)

    S = runs.rgrid;

    zrmat = runs.rgrid.z_r;
    f = runs.rgrid.f';

    extract_params(runs);

    % cylindrical co-ordinates
    r0 = eddy.dia/2; % This is required to account for exponential decay
    [th,r] = cart2pol((S.x_rho-eddy.cx),(S.y_rho-eddy.cy));
    rnorm  = r./r0; % normalized radius
    eddy.ix = find_approx(S.x_rho(:,1),eddy.cx,1);
    eddy.iy = find_approx(S.y_rho(1,:),eddy.cy,1);

    eddy.xyprof = nan(size(rnorm));
    % eddy temp. field - xy profile - Reference: Katsman et al. (2003)
    if flags.solidbody_katsman % solid body rotation
        % rnorm = rnorm .* 2;
        % normalized rm = rm_tilde in Katsman et al. (2003)
        % values determined by matching first and second derivatives of
        % temp. profile at rm - see paper for details
        rmnorm = nthroot( (eddy.a-2)/(eddy.a-1) , eddy.a);
        gamma = -2 * (eddy.a-2)./eddy.a ./ (rmnorm)^2;

        exponent = (eddy.a - 1)/eddy.a .* (rnorm.^(eddy.a) - rmnorm.^(eddy.a));
        eddy.xyprof = (gamma/2 * rnorm.^2 + 1) .* (rnorm <= rmnorm) ...
                       + (gamma/2 *rmnorm^2 + 1) .* exp( -1 * exponent ) .* ...
                                                           (rnorm > rmnorm);
    end
    eddy.xyprof = eddy.xyprof./max(eddy.xyprof(:));
    eddy.xyprof = repmat(eddy.xyprof', [1 1 S.N]);

    eddy.temp = eddy.tamp .* eddy.xyprof .* exp(- (zrmat./ ...
                                                   eddy.depth).^2);

    % Radial
    if max(~isnan(eddy.temp(:))) % && flags.use_radial

        % Calculate azimuthal velocity shear (r d(theta)/dt)_z using geostrophic balance
        if flags.solidbody_katsman
            dTdr = gamma * rnorm ./ r0 .* (rnorm <= rmnorm) ...
                    + (gamma/2 .* rmnorm^2 + 1) .*  (-(eddy.a-1) ./ r0 .* rnorm.^(eddy.a-1)) ...
                               .* exp(-exponent) .* (rnorm > rmnorm);
        end
        dTdr = dTdr';

        % loop through points and integrate
        sz = size(eddy.xyprof);
        itz = nan(sz);
        eddy.zeta = nan(sz(1:2));
        tic;
        for ii=1:size(eddy.xyprof, 1)
            for jj=1:size(eddy.xyprof, 2)
                eddy.zeta(ii,jj) = phys.TCOEF * trapz(zrmat(:,jj,ii), ...
                                                      eddy.temp(ii,jj,:) ...
                                                      - eddy.temp(eddy.ix,eddy.iy,:),3);

                eddy.tz = exp(-(zrmat(:,jj,ii) ./ eddy.depth).^2);
                itz(ii,jj,:) = cumtrapz(zrmat(:,jj,ii), eddy.tz, 1);
            end
        end
        toc;
        eddy.zeta = eddy.zeta - min(eddy.zeta(:));

        % azimuthal velocity = r d(theta)/dt
        rut = zeros(size(eddy.xyprof));
        rut = eddy.tamp * (phys.TCOEF*phys.g) .* bsxfun(@times, itz, dTdr./f);

        % solve quadratic for vel. if gradient wind balance
        vgeo = rut;
        if flags.use_gradient
            rfb2 = r'.*f ./ 2;
            sdisc = sqrt(1 + bsxfun(@times,vgeo,2./rfb2));% sqrt(discriminant)
            if isreal(sdisc) % gradient wind doesn't always work with anticyclones
                rut = bsxfun(@times,(-1 + sdisc), rfb2);
                disp('Using gradient wind balance.');
            else
                % cyclostrophic balance doesn't work yet
                error(['gradient wind calculated complex v! - ' ...
                         'Ro > 0.25']);
            end
        end

        eddy.u = -1 * bsxfun(@times,rut, sin(th'));
        eddy.v =      bsxfun(@times,rut, cos(th'));
end