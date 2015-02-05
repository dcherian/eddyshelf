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

    locs = [sx1 sx2 2 sz(1)-1 ceil(sx2+(sz(1)-sx2)/2)];
    %locations = [locations 'resistance | '];
    locations = [locations 'sponge | '];
    locations = [locations 'domain edge | '];
    %locs = loc1;
    ax = 'x';
    ax2 ='y';

    % y (ax2) limits - between shelfbreak & edge of sponge
    sy1 = runs.bathy.isb;
    sy2 = find(runs.sponge(sz(1)/2, :) == 1, 1, 'first') - 1;
    locations = [locations 'between shelfbreak & edge of sponge in ' ...
                        'y'];

    runs.asflux = [];

    assert(runs.bathy.axis == 'y');

    % mask based on eddy contour's northern edge
    % permute to (y,z,t) like eddmask later.
    %nedge = permute(bsxfun(@lt, runs.rgrid.y_rho(:,1), runs.eddy.vor.ne), ...
    %                [1 3 2]);

    % read data
    % (y,z,t, location)
    % preallocate for speed
    tind = [1 length(runs.eddy.t)];

    u = nan([sy2-sy1+1 runs.rgrid.N tind(2)-tind(1)+1 length(locs)]);
    rho = u;
    eddmask = u;
    if do_energy
        v = nan([sy2-sy1+1 runs.rgrid.N tind(2)-tind(1)+1 length(locs)]);
    end

    % loop over locations - read data
    for ii=1:length(locs)
        % works for rho points
        volumer = {ax locs(ii) locs(ii); ...
                   ax2 sy1 sy2};
        % for u points
        volumeu = {ax locs(ii)-1 locs(ii); ...
                  ax2 sy1 sy2};
        % for v points
        volumev = {ax locs(ii) locs(ii); ...
                  ax2 sy1-1 sy2};

        u(:,:,:,ii) = squeeze(avg1(dc_roms_read_data(runs.dir, ...
                                                     'u', tind, ...
                                                     volumeu, ...
                                                     [], runs.rgrid, ...
                                                     'his', 'single'),1));
        rho(:,:,:,ii) = dc_roms_read_data(runs.dir, 'rho', tind, ...
                                          volumer, [], runs.rgrid, ...
                                          'his', 'single');

        % mask based on dye
        eddmask(:,:,:,ii) = (dc_roms_read_data(runs.dir, runs.eddname, ...
                                               tind, volumer, [], runs.rgrid, ...
                                          'his', 'single') > runs.eddy_thresh);

        %eddmask(:,:,:,ii) = bsxfun(@and, eddmask(:,:,:,ii), ...
        %                           nedge);

        if do_energy
            v(:,:,:,ii) = avg1(dc_roms_read_data(runs.dir, 'v', tind, ...
                                            volumev, [], runs.rgrid, ...
                                            'his', 'single'), 1);
            %rback(:,:,1,ii) = dc_roms_read_data(runs.dir, 'rho', [1 1], {ax ...
            %                    locs(ii) locs(ii)}, [], runs.rgrid, ...
            %                                    'his', 'single');
        end
    end

    % full rho field
    rho = rho + 1000;

    % calculate mass flux ρu
    rflux = rho .* u .* eddmask;

    % calculate energy flux
    if do_energy
        % at interior RHO points
        ke = 1/2 * (u.^2 + v.^2) .* rho;
        pe = - runs.params.phys.g * bsxfun(@times, rho, ...
                                           runs.rgrid.z_r(:,sy1:sy2,1)');

        % calculate pu term
        % dz = (y,z) > 0
        dzmat = permute(diff(runs.rgrid.z_w(:,sy1:sy2,1), 1, 1), [2 1]);
        % p = ∫_ζ^z ρg dz (no negative sign)
        % p = (y,z,t)
        pres = 9.81 * bsxfun(@times, rho, dzmat);
        % flip because I'm integrating from surface to (z)
        pres = flip(pres, 2);
        % integrate
        pres = cumsum(pres, 2);
        % flip back
        pres = flip(pres, 2);

        % energy flux (y, z, t, location)
        peflux = u .* pe;
        keflux = u .* ke + u .* pres;

        %peflux = bsxfun(@times, u(2:end-1,:,:,:) .* (pe(2:end-1,:,:,:)), ...
        %                nedge(2:end-1,:,:));
        %keflux = bsxfun(@times, u(2:end-1,:,:,:) .* (ke), nedge(2:end-1,:,:));

        peflux_edd = peflux .* eddmask;
        keflux_edd = keflux .* eddmask;
    end

    % integrate flux
    disp('Looping over locations - integrating energy / mass fluxes');
    % preallocate for speed
    sz2 = size(rho);
    sz2 = [sz2 1];
    sz2(1) = sz2(1);
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
    % loop over locations
    for ii=1:length(locs)
        % grid vectors for integration.
        zvec = runs.rgrid.z_r(:,sy1:sy2,locs(ii))';
        yvec = runs.rgrid.y_rho(sy1:sy2,locs(ii));

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

        % integrate in y
        irflux(:,ii) = squeeze(trapz(yvec, irfluxyt(:,:,ii), 1));

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

            % correct for sy1 = isb by default
            isl = runs.bathy.isl - sy1 + 1;
            % these ranges are correct! even if you don't think so.
            % ikeflux = deep.ikeflux + topo.ikeflux
            ideep = [isl:length(yvec)];
            itopo = [1:isl];

            % deep water
            runs.asflux.deep.ikeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                ikefluxyt(ideep,:,ii), 1));
            runs.asflux.deep.ipeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                ipefluxyt(ideep,:,ii), 1));
            runs.asflux.deep.eddy.ikeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                ikefluxyt_edd(ideep,:,ii), 1));
            runs.asflux.deep.eddy.ipeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                ipefluxyt_edd(ideep,:,ii), 1));

            % over shelf-slope
            runs.asflux.topo.ikeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                ikefluxyt(itopo,:,ii), 1));
            runs.asflux.topo.ipeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                ipefluxyt(itopo,:,ii), 1));
            runs.asflux.topo.eddy.ikeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                ikefluxyt_edd(itopo,:,ii), 1));
            runs.asflux.topo.eddy.ipeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                ipefluxyt_edd(itopo,:,ii), 1));
        end
    end
    toc;

    % save data
    runs.asflux.comment = ['ix = indices | x = x-locations (m) ' ...
                        '| ik(p)eflux(t, locations) - integrated ' ...
                        'KE/PE flux | (r)flux - mass fluxes | **' ...
                        'fluxyt (y, t, locations) - depth-integrated fluxes' ...
                        ' | .eddy. = eddye mask applied | .deep. = integrated' ...
                        ' over deep only | .topo. = integrated over shelf-slope'];
    runs.asflux.hash = githash([mfilename('fullpath') '.m']);
    runs.asflux.time = runs.time;
    runs.asflux.locations = locations;
    runs.asflux.ix = locs;
    runs.asflux.iy = [sy1 sy2];
    runs.asflux.yvec = yvec;
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