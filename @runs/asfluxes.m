% This functions calculates downstream mass and energy fluxes that
% make up the "eddy leakage"
% saves output in runs.asflux;

function asfluxes(runs)
    ticstart = tic;

    do_energy = 0;

    % locations - grid indices - RHO points
    % find location at ndtime == 1, then add 2 rossby radii
    loc1 = find_approx(runs.rgrid.x_rho(1,:), ...
                       runs.eddy.mx(find_approx( ...
                           runs.eddy.t*86400/runs.tscale, 1)) + ...
                       2 * runs.rrdeep);
    locs = [loc1 runs.params.grid.ixn runs.params.grid.ixp]

    ax = 'x';

    % read data
    % (y,z,t, location)
    for ii=1:length(locs)
        u(:,:,:,ii) = squeeze(avg1(dc_roms_read_data(runs.dir, ...
                                                     'u', [], ...
                                                     {ax locs(ii) locs(ii)+1}, ...
                                                     [], runs.rgrid, ...
                                                     'his', 'single'),1));
        rho(:,:,:,ii) = dc_roms_read_data(runs.dir, 'rho', [], {ax ...
                            locs(ii) locs(ii)}, [], runs.rgrid, ...
                                          'his', 'single');
        if do_energy
            v(:,:,:,ii) = dc_roms_read_data(runs.dir, 'v', [], {ax ...
                                locs(ii) locs(ii)}, [], runs.rgrid, ...
                                            'his', 'single');
            rback(:,:,1,ii) = dc_roms_read_data(runs.dir, 'rho', [1 1], {ax ...
                                locs(ii) locs(ii)}, [], runs.rgrid, ...
                                                'his', 'single');
        end
    end

    % calculate mass flux ρu
    rflux = rho .* u;

    % calculate energy flux
    if do_energy
        % at interior RHO points
        ke = 1/2 * (u(2:end-1,:,:,:).^2 + avg1(v,1).^2);
        pe = - runs.params.phys.g * ...
             bsxfun(@times, bsxfun(@minus, rho, rback), ...
                    runs.rgrid.z_r(:,:,1)');

        % energy flux (y, z, t, location)
        peflux = u(2:end-1,:,:,:) .* (pe(2:end-1,:,:,:));
        keflux = u(2:end-1,:,:,:) .* (ke);
    end

    % integrate flux
    disp('Looping over locations - integrating energy flux');
    tic;
    for ii=1:length(locs)
        % grid vectors for integration.
        zvec = runs.rgrid.z_r(:,2:end-1,locs(ii))';
        yvec = runs.rgrid.y_rho(2:end-1,locs(ii));

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
            end
        end

        irflux(:,ii) = squeeze(trapz(yvec, irfluxyt(:,:,ii), 1));

        if do_energy
            % integrated flux (t, location)
            ikeflux(:,ii) = squeeze(trapz(yvec, ikefluxyt(:,:,ii), ...
                                          1));
            ipeflux(:,ii) = squeeze(trapz(yvec, ipefluxyt(:,:,ii), ...
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
    runs.asflux.ix = locs;
    runs.asflux.x = runs.rgrid.x_rho(1,locs);
    runs.asflux.irflux = irflux;
    runs.asflux.irfluxyt = irfluxyt;

    if do_energy
        runs.asflux.ikeflux = ikeflux;
        runs.asflux.ipeflux = ipeflux;
        runs.asflux.ikefluxyt = ikefluxyt;
        runs.asflux.ipefluxyt = ipefluxyt;
    end

    asflux = runs.asflux;
    save([runs.dir '/fluxes.mat'], 'asflux', '-append');

    disp('Completed leakage calculations');
    toc(ticstart);
end
