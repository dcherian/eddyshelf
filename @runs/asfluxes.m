% This functions calculates downstream mass and energy fluxes that
% make up the "eddy leakage"
% saves output in runs.asflux;

function asfluxes(runs)
    ticstart = tic;

    do_energy = 1;

    locations = '';

    % locations - grid indices - RHO points
    % find "resistance" location, then add 2 rossby radii
    [xx,yy,tind] = locate_resistance(runs);
    loc1 = find_approx(runs.rgrid.x_rho(1,:), ...
                       xx + 2 * runs.rrdeep);

    % find sponge edges
    sz = size(runs.sponge);
    sx1 = find(runs.sponge(1:sz(1)/2,sz(2)/2) == 0, 1, 'first');
    sx2 = sz(1)/2 + find(runs.sponge(sz(1)/2:end,sz(2)/2) == 1, 1, ...
                         'first') - 2;
    sy1 = find(runs.sponge(sz(1)/2, :) == 1, 1, 'first') - 1;

    locs = [loc1 sx1 sx2];
    locations = [locations 'resistance | '];
    locations = [locations 'sponge | '];
    %locs = loc1;
    ax = 'x';

    assert(runs.bathy.axis == 'y');

    % mask based on eddy contour's northern edge
    % permute to (y,z,t) like eddmask later.
    %nedge = permute(bsxfun(@lt, runs.rgrid.y_rho(:,1), runs.eddy.vor.ne), ...
    %                [1 3 2]);

    % read data
    % (y,z,t, location)
    % preallocate for speed

    tind = [1 length(runs.time)];
    tind = [300 400];

    u = nan([sz(2) runs.rgrid.N tind(2)-tind(1)+1 length(locs)]);
    rho = u;
    eddmask = u;
    if do_energy
        v = nan([sz(2)-1 runs.rgrid.N tind(2)-tind(1)+1 length(locs)]);
    end

    for ii=1:length(locs)
        u(:,:,:,ii) = squeeze(avg1(dc_roms_read_data(runs.dir, ...
                                                     'u', tind, ...
                                                     {ax locs(ii) locs(ii)+1}, ...
                                                     [], runs.rgrid, ...
                                                     'his', 'single'),1));
        rho(:,:,:,ii) = dc_roms_read_data(runs.dir, 'rho', tind, {ax ...
                            locs(ii) locs(ii)}, [], runs.rgrid, ...
                                          'his', 'single');

        % mask based on dye
        eddmask(:,:,:,ii) = (dc_roms_read_data(runs.dir, runs.eddname, tind, {ax ...
                            locs(ii) locs(ii)}, [], runs.rgrid, ...
                                          'his', 'single') > runs.eddy_thresh);

        %eddmask(:,:,:,ii) = bsxfun(@and, eddmask(:,:,:,ii), ...
        %                           nedge);

        if do_energy
            v(:,:,:,ii) = dc_roms_read_data(runs.dir, 'v', tind, {ax ...
                                locs(ii) locs(ii)}, [], runs.rgrid, ...
                                            'his', 'single');
            %rback(:,:,1,ii) = dc_roms_read_data(runs.dir, 'rho', [1 1], {ax ...
            %                    locs(ii) locs(ii)}, [], runs.rgrid, ...
            %                                    'his', 'single');
        end
    end

    rho0 = runs.params.phys.rho0;
    % full rho field
    rho = rho + rho0;

    % calculate mass flux ρu
    rflux = rho .* u .* eddmask;

    % calculate energy flux
    if do_energy
        % at interior RHO points
        ke = 1/2 * (u(2:end-1,:,:,:).^2 + avg1(v,1).^2) .* ...
             (rho(2:end-1,:,:,:));
        pe = - runs.params.phys.g * bsxfun(@times, rho, ...
                                           runs.rgrid.z_r(:,:,1)');

        % calculate pu term
        % dz = (y,z) > 0
        dzmat = permute(diff(runs.rgrid.z_w(:,:,1), 1, 1), [2 1]);
        % p = ∫_ζ^z ρg dz (no negative sign)
        % p = (y,z,t)
        pres = 9.81 * bsxfun(@times, rho, dzmat);
        % flip because I'm integrating from surface to (z)
        pres = flipdim(pres, 2);
        % integrate
        pres = cumsum(pres, 2);
        % flip back
        pres = flipdim(pres, 2);

        % energy flux (y, z, t, location)
        peflux = u(2:end-1,:,:,:) .* (pe(2:end-1,:,:,:));
        keflux = u(2:end-1,:,:,:) .* ke + ...
                 u(2:end-1,:,:,:) .* pres(2:end-1,:,:,:);

        %peflux = bsxfun(@times, u(2:end-1,:,:,:) .* (pe(2:end-1,:,:,:)), ...
        %                nedge(2:end-1,:,:));
        %keflux = bsxfun(@times, u(2:end-1,:,:,:) .* (ke), nedge(2:end-1,:,:));

        peflux_edd = peflux .* eddmask(2:end-1,:,:,:);
        keflux_edd = keflux .* eddmask(2:end-1,:,:,:);
    end

    % integrate flux
    disp('Looping over locations - integrating energy / mass fluxes');
    % preallocate for speed
    sz2 = size(rho);
    sz2 = [sz2 1];
    sz2(1) = sz2(1) - 2; % 2:end-1 for y
    irfluxyt = nan([sz2(1) sz2(3) sz2(4)]);
    ikefluxyt = irfluxyt;
    ipefluxyt = irfluxyt;
    ikefluxyt_edd = irfluxyt;
    ipefluxyt_edd = irfluxyt;

    irflux = nan([sz2(3) sz2(4)]);
    ikeflux = irflux;
    ipeflux = irflux;
    ikeflux_edd = irflux;
    ipeflux_edd = irflux;

    tic;
    for ii=1:length(locs)
        % grid vectors for integration.
        zvec = runs.rgrid.z_r(:,2:end-1,locs(ii))';
        yvec = runs.rgrid.y_rho(2:end-1,locs(ii));

        % integrate vertically
        for jj=1:size(yvec,1)
            irfluxyt(jj,:,ii) = squeeze(trapz(zvec(jj,:), ...
                                              rflux(jj,:,:,ii), 2));
            if do_energy
                % integrated flux (y, t, location)
                ikefluxyt(jj,:,ii) = squeeze(trapz(zvec(jj,:), ...
                                                   keflux(jj,:,:,ii), 2));
                % integrated flux (y, t, location)
                ipefluxyt(jj,:,ii) = squeeze(trapz(zvec(jj,:), ...
                                                   peflux(jj,:,:,ii), ...
                                                   2));
                ikefluxyt_edd(jj,:,ii) = squeeze(trapz(zvec(jj,:), ...
                                                   keflux_edd(jj,:,:,ii), 2));
                % integrated flux (y, t, location)
                ipefluxyt_edd(jj,:,ii) = squeeze(trapz(zvec(jj,:), ...
                                                   peflux_edd(jj,:,:,ii), ...
                                                   2));
            end
        end

        irflux(:,ii) = squeeze(trapz(yvec, irfluxyt(:,:,ii), 1));

        % integrate in y
        if do_energy
            % integrated flux (t, location)
            ikeflux(:,ii) = squeeze(trapz(yvec, ikefluxyt(:,:,ii), ...
                                          1));
            ipeflux(:,ii) = squeeze(trapz(yvec, ipefluxyt(:,:,ii), ...
                                          1));

            ikeflux_edd(:,ii) = squeeze(trapz(yvec, ikefluxyt_edd(:,:,ii), ...
                                          1));
            ipeflux_edd(:,ii) = squeeze(trapz(yvec, ipefluxyt_edd(:,:,ii), ...
                                          1));
        end
    end
    toc;

    % save data
    runs.asflux.comment = ['ix = indices | x = x-locations (m) ' ...
                        '| ik(p)eflux(t, locations) - integrated ' ...
                        'KE/PE flux | (r)flux - mass fluxes | **' ...
                        'fluxyt (y, t, locations) - depth-integrated fluxes'];
    runs.asflux.hash = githash;
    runs.asflux.time = runs.time;
    runs.asflux.locations = locations;
    runs.asflux.ix = locs;
    runs.asflux.x = runs.rgrid.x_rho(1,locs);
    runs.asflux.eddy.irflux = irflux;
    runs.asflux.eddy.irfluxyt = irfluxyt;

    if do_energy
        runs.asflux.ikeflux = ikeflux;
        runs.asflux.ipeflux = ipeflux;
        runs.asflux.ikefluxyt = ikefluxyt;
        runs.asflux.ipefluxyt = ipefluxyt;

        % eddy water
        runs.asflux.eddy.ikeflux = ikeflux_edd;
        runs.asflux.eddy.ipeflux = ipeflux_edd;
        runs.asflux.eddy.ikefluxyt = ikefluxyt_edd;
        runs.asflux.eddy.ipefluxyt = ipefluxyt_edd;
    end

    asflux = runs.asflux;
    save([runs.dir '/fluxes.mat'], 'asflux', '-append');

    disp('Completed leakage calculations');
    toc(ticstart);
end
