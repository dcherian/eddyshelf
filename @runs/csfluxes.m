% quantify cross-shelfbreak and along-shelfbreak fluxes of a whole
% bunch of stuff:
function [] = csfluxes(runs, ftype)

% Things I want to calculate fluxes of:
% 1. RV & PV
% 2. Shelf water dye & eddy dye
    ticstart = tic;

    disp('===================')
    disp([runs.name '.csfluxes']);
    disp('===================')

    reset = 1; % 1 means redo all isobaths

    % Use history or avg files?
    if ~exist('ftype', 'var') || isempty(ftype), ftype = 'his'; end

    precision = 'single';

    % calculate energy fluxes?
    do_energy = 0;

    % debug?
    debug = 0;

    % do PV fluxes?
    dopv = 0;

    vorname = [runs.dir '/ocean_vor.nc'];

    tstart = 1;
    % revind = runs.eddy.trevind;

    % quantities for area-averaging if needed
    h = runs.bathy.h(2:end-1,2:end-1);

    % sort out isobaths across which to calculate transports
    Y = max(runs.rgrid.y_rho(:));

    % is cyclone? if so, -1; else 0
    sgntamp = runs.sgntamp;

    xsb = runs.bathy.xsb; xsl = runs.bathy.xsl;

    % Non-dimensional isobath = (y_{isobath} - y_{sb})/(Eddy center
    % location)

    if reset
        runs.csflux = [];
        if runs.bathy.axis == 'y'
            [~,R,restind] = runs.locate_resistance;
        else
            [R,~,restind] = runs.locate_resistance;
        end
        if isempty(R)
            error(['locate_resistance didn''t return anything for ' runs.name]);
        end
        R = abs(R - xsb);
        runs.csflux.ndloc = linspace(0,2,13);
        runs.csflux.R = R;
        loc = runs.csflux.ndloc * R * sgntamp + xsb;
        if runs.params.flags.flat_bottom
            % remove the wall location
            runs.csflux.ndloc(1) = [];
            loc(1) = [];
        end

        if debug
            runs.animate_field('zeta',[],restind,1);
            if runs.bathy.axis == 'y'
                liney(loc/1000, [], 'r');
            else
                linex(loc/1000, [], 'r');
            end
        end
    else
        R = runs.csflux.R;
        loc = runs.csflux.ndloc * R * sgntamp + xsb;
    end

    if runs.params.bathy.axis == 'x'
        csvelid = 'u';
        asvelid = 'v';
        bathyax = 1;
        xvec = runs.rgrid.yr(1,:);
        indices = vecfind(runs.rgrid.x_rho(1,:),loc);

        % append sponge edge to existing values
        sy2 = runs.spng.sx2;
    else
        csvelid = 'v';
        asvelid = 'u';
        bathyax = 2;
        indices = vecfind(runs.rgrid.y_rho(:,1),loc);
        xvec = runs.rgrid.xr(:,1);

        % append sponge edge to existing values
        sy2 = runs.spng.sy2;
    end

    sz = size(runs.sponge);

    if reset
        % save locations
        runs.csflux.x = loc;
        % save indices for locations - w.r.t INTERIOR RHO POINTS
        runs.csflux.ix = indices;
        % save isobath values
        runs.csflux.h = ceil(runs.bathy.h(1,runs.csflux.ix));

        runs.csflux.comment = ['shelf / slope / eddy= which water mass am I ', ...
                            ' targeting? |\n (x,ix,h) = (location, index, depth) ' ...
                            'at which I''m calculating transport |\n ' ...
                            'Slope(kk) = water onshore of loc(kk)'];
    end

    % how much of the time vector should I read?
    if isfield(runs.eddy, 'tend') & (runs.eddy.tend ~= 0)
        tinf = runs.eddy.tend;
    else
        tinf = length(runs.time);
    end

    time = runs.time(1:tinf);
    runs.csflux.time = time;
    t0 = 1;

    % size for initialization
    nloc = length(loc);
    szfull = size(runs.bathy.h);
    szflux = [tinf nloc nloc];
    if bathyax == 1
        szfluxxt = [szfull(2)-2 tinf nloc];
    else
        szfluxxt = [szfull(1)-2 tinf nloc];
    end

    if reset
        % initialize - mass flux variables
        runs.csflux.off.slope = nan(szflux);
        runs.csflux.off.eddy = nan(szflux);

        runs.csflux.off.itrans.slope = nan(szflux);
        runs.csflux.off.itrans.eddy = nan(szflux);

        runs.csflux.on.slope = nan(szflux);
        runs.csflux.on.eddy = nan(szflux);

        runs.csflux.on.itrans.slope = nan(szflux);
        runs.csflux.on.itrans.eddy = nan(szflux);

        runs.csflux.off.bins = cell([1 nloc]);

        runs.csflux.off.slopewater.vmax = nan([tinf nloc]);

        runs.csflux.off.slopewater.trans = cell([1 nloc nloc]);
        runs.csflux.on.slopewater.vtrans = cell([1 nloc nloc]);
        runs.csflux.off.slopewater.itrans = cell([1 nloc nloc]);
        runs.csflux.off.slopewater.envelope = nan([tinf nloc]);

        runs.csflux.on.slopewater.trans = cell([1 nloc nloc]);
        runs.csflux.on.slopewater.vtrans = cell([1 nloc nloc]);
        runs.csflux.on.slopewater.itrans = cell([1 nloc nloc]);
        runs.csflux.on.slopewater.envelope = nan([tinf nloc]);

        runs.csflux.slopext  = nan(szfluxxt);
        runs.csflux.eddyxt   = nan(szfluxxt);
        runs.csflux.offmask = logical(zeros(szfluxxt));
        runs.csflux.onmask = logical(zeros(szfluxxt));

        runs.csflux.off.slopezt = nan([runs.rgrid.N tinf nloc nloc]);
        runs.csflux.off.eddyzt = nan([runs.rgrid.N tinf nloc nloc]);
        runs.csflux.on.slopezt = nan([runs.rgrid.N tinf nloc nloc]);
        runs.csflux.on.eddyzt = nan([runs.rgrid.N tinf nloc nloc]);

        runs.csflux.off.slopewater.vertitrans = nan([runs.rgrid.N nloc]);
        runs.csflux.off.eddywater.vertitrans = nan([runs.rgrid.N nloc]);
        runs.csflux.on.slopewater.vertitrans = nan([runs.rgrid.N nloc]);
        runs.csflux.on.eddywater.vertitrans = nan([runs.rgrid.N nloc]);

        % initialize - pv fluxes
        if exist(vorname, 'file') && dopv == 1
            pvtime = ncread(vorname, 'ocean_time');
            if isequal(pvtime', time) || isequal(pvtime, time)
                dopv = 1;
                runs.csflux.off.pv = nan(szflux);
                runs.csflux.off.rv = nan(szflux);
                runs.csflux.on.pv = nan(szflux);
                runs.csflux.on.rv = nan(szflux);
            end
        end

        % initialize - energy fluxes
        if do_energy
            % total
            runs.csflux.ikeflux = nan(szflux);
            runs.csflux.ipeflux = nan(szflux);
            runs.csflux.ikefluxxt = nan(szfluxxt);
            runs.csflux.ipefluxxt = nan(szfluxxt);

            % masked by eddye
            runs.csflux.eddy.ikeflux = nan(szflux);
            runs.csflux.eddy.ipeflux = nan(szflux);
            runs.csflux.eddy.ikefluxxt = nan(szfluxxt);
            runs.csflux.eddy.ipefluxxt = nan(szfluxxt);
        end
    end

    % sponge mask
    if runs.bathy.axis == 'y'
        spongemask = ~runs.sponge(2:end-1,sz(2)/2);
    else
        spongemask = ~runs.sponge(sz(1)/2,2:end-1)';
    end

    % for integrated transport diagnostics
    dt = [time(2)-time(1) diff(time)'];
    if bathyax == 2
        dx = 1./runs.rgrid.pm(1,2:end-1)';
    else
        dx = 1./runs.rgrid.pn(2:end-1,1);
    end

    % loop over all isobaths
    for kk=1:length(indices)
        disp(['Doing isobath ' num2str(kk) '/', ...
              num2str(length(indices))]);

        volr = {runs.bathy.axis runs.csflux.ix(kk) runs.csflux.ix(kk)};
        volv = {runs.bathy.axis runs.csflux.ix(kk)-1 runs.csflux.ix(kk)};

        % read along-shore section of cross-shore vel.
        % dimensions = (x/y , z , t )
        % average carefully to get values at RHO point
        csvel = squeeze(avg1(dc_roms_read_data(runs.dir, csvelid, [t0 tinf], ...
                                               volv, [], runs.rgrid, ftype, ...
                                               precision), bathyax));
        csvel = csvel(2:end-1,:,:,:) * sgntamp;

        % process cross-shelf dye
        csdye = dc_roms_read_data(runs.dir, runs.csdname, [t0 tinf], ...
                                  volr, [], runs.rgrid, ftype, 'single');
        csdye = csdye(2:end-1,:,:);

        % read eddye & define eddy water mass
        eddymask = dc_roms_read_data(runs.dir, runs.eddname, [t0 tinf], ...
                                     volr, [], runs.rgrid, ftype, 'single') > runs.eddy_thresh;
        eddymask = eddymask(2:end-1,:,:);

        % offshore and onshore masks
        if bathyax == 2
            cxi = runs.eddy.mx(tstart:end);
        else
            cxi = runs.eddy.my(tstart:end);
        end

        if bathyax == 2
            offmask = bsxfun(@times, ...
                             bsxfun(@lt, runs.eddy.xr(:,1), cxi(t0:tinf)), ...
                             spongemask);
        else
            offmask = bsxfun(@times, ...
                             bsxfun(@gt, runs.eddy.yr(1,:)', cxi(t0:tinf)), ...
                             spongemask);
        end
        onmask = bsxfun(@times, 1 - offmask, spongemask);

        offmask = logical(offmask); onmask = logical(onmask);

        if kk ~= 1
            % New approach:
            % +++++++++++++
            % Go to first zero-crossing (i.e, do only offshore velocities).
            % This is required only right at the shelfbreak. For ndloc = 0.17, this is not an
            % issue. There, the eddy center based mask works well.
            %
            % Older approach:
            % +++++++++++++++
            % cxi at shelfbreak = eastern edge because export generally occurs west of the eastern edge
            % dimensions - (along-shore) x time
            % The problem is that this is sensitive to contour detection.

            csvonmask = 1;
            csvoffmask = 1;
        else
            csvonmask = 1;
            csvoffmask = 1;
            offmask = repmat(spongemask, [1 size(runs.csflux.offmask, 2)]);
            onmask = repmat(spongemask, [1 size(runs.csflux.offmask, 2)]);

            flag_found = 1;
            % % define a mask based on first zero-crossing at each depth
            % csvoffmask = csvel > 0;
            % % pick out region with eddy center
            % for tt=1:size(csvel, 3)
            %     flag_found(tt) = 0;
            %     regions = bwconncomp(csvoffmask(:,:,tt), 8);
            %     ix = find_approx(xvec, cxi(tt), 1);
            %     iwe = find_approx(xvec, runs.eddy.rhovor.we(tt), 1);

            %     % find region with eddy center
            %     for zz=1:regions.NumObjects
            %         maskreg = zeros(regions.ImageSize);
            %         maskreg(regions.PixelIdxList{zz}) = 1;

            %         if any(maskreg(ix,:) == 1)
            %             csvoffmask(:,:,tt) = maskreg;
            %             flag_found(tt) = 1;
            %             break;
            %         end
            %     end

            %     % sometimes eddy center is a few grid points off,
            %     % then look for midpoint of western edge and center
            %     if ~flag_found(tt)
            %         for zz=1:regions.NumObjects
            %             maskreg = zeros(regions.ImageSize);
            %             maskreg(regions.PixelIdxList{zz}) = 1;

            %             if any(maskreg(ceil((ix+iwe)/2),:) == 1)
            %                 csvoffmask(:,:,tt) = maskreg;
            %                 flag_found(tt) = 1;
            %                 break;
            %             end
            %         end
            %     end

            %     if ~flag_found(tt)
            %         for zz=1:regions.NumObjects
            %             maskreg = zeros(regions.ImageSize);
            %             maskreg(regions.PixelIdxList{zz}) = 1;

            %             if any(maskreg(ix-5,:) == 1)
            %                 csvoffmask(:,:,tt) = maskreg;
            %                 flag_found(tt) = 1;
            %                 break;
            %             end
            %         end
            %     end
            % end

            % csvelbar = squeeze(avg1(dc_roms_read_data(runs.dir, [csvelid 'bar'], [t0 tinf], ...
            %                                           volv, [], runs.rgrid, ftype, ...
            %                                           precision), bathyax));
            % csvelbar = csvelbar(2:end-1,:,:,:) * sgntamp;

            % % apply other refinements
            % csvoffmask = bsxfun(@and, csvoffmask, permute(csvelbar > 0, [1 3 2]));
            % csvoffmask = bsxfun(@and, csvoffmask, spongemask);
            % csvoffmask = bsxfun(@or, csvoffmask, permute(offmask, [1 3 2]));
            % clear maskreg

            % csvonmask = ~csvoffmask;
            % offmask = 1;
            % onmask  = 1;

            % save here because the variables get over-written later.
            runs.csflux.csvoffmask = csvoffmask;
            runs.csflux.csvonmask = csvonmask;
        end

        runs.csflux.offmask(:,:,kk) = logical(offmask);
        runs.csflux.onmask(:,:,kk) = logical(onmask);

        if bathyax == 2
            zvec = runs.rgrid.z_r(:, runs.csflux.ix(kk), 1);
        else
            zvec = runs.rgrid.z_r(:, 1, runs.csflux.ix(kk));
        end
        runs.csflux.vertbins(:,kk) = zvec;

        runs.csflux.eddyxt(:,:,kk) = squeeze(trapz(zvec, eddymask .* csvel,2));
        runs.csflux.off.eddyxt(:,:,kk) = squeeze(trapz(zvec, eddymask .* csvel .* csvoffmask,2));
        runs.csflux.on.eddyxt(:,:,kk)  = squeeze(trapz(zvec, eddymask .* csvel .* csvonmask ,2));

        runs.csflux.off.eddyzt(:,:,kk) = ...
            squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(eddymask .* csvel .* csvoffmask), ...
                                                 permute(offmask, [1 3 2])), dx),1));
        runs.csflux.on.eddyzt(:,:,kk) = ...
            squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(eddymask .* csvel .* csvonmask), ...
                                                 permute(onmask, [1 3 2])), dx),1));

        % transports as a function of z only
        runs.csflux.off.eddywater.vertitrans(:,kk) = ...
            nansum(bsxfun(@times, runs.csflux.off.eddyzt(:,:,kk), dt), 2);
        runs.csflux.on.eddywater.vertitrans(:,kk) = ...
            nansum(bsxfun(@times, runs.csflux.on.eddyzt(:,:,kk), dt), 2);

        % calculate fluxes
        if runs.csflux.ix(kk) ~= sy2
            runs.csflux.off.eddy(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.off.eddyxt(:,:,kk) .* offmask, dx),1))';

            runs.csflux.on.eddy(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.on.eddyxt(:,:,kk) .* onmask, dx),1))';

            [runs.csflux.off.itrans.eddy(:,kk), ...
             runs.csflux.off.avgflux.eddy(kk)] = ...
                runs.integrate_flux(time, runs.csflux.off.eddy(:,kk));

            [runs.csflux.on.itrans.eddy(:,kk), ...
             runs.csflux.on.avgflux.eddy(kk)] = ...
                runs.integrate_flux(time, runs.csflux.on.eddy(:,kk));
        end

        if do_energy
            disp('Calculating energy...');
            error('CHECK INDICES OF READING');
            % read velocity and calculate pressure etc.
            asvel = avg1(dc_roms_read_data(runs.dir, asvelid, [t0 tinf], ...
                                           {runs.bathy.axis runs.csflux.ix(kk) ...
                                runs.csflux.ix(kk)}, [], runs.rgrid, ftype, 'single'), 1);

            % read velocity and calculate pressure etc.
            rho = dc_roms_read_data(runs.dir, 'rho', [t0 tinf], ...
                                      {runs.bathy.axis runs.csflux.ix(kk) ...
                                runs.csflux.ix(kk)}, [], runs.rgrid, ftype, 'single');
            rho = rho(2:end-1,:,:);

            ke = 1/2 * rho .* (csvel.^2 + asvel.^2);

            pe = - runs.params.phys.g * bsxfun(@times, rho, ...
                                               squeeze(runs.rgrid.z_r(:,runs.csflux.ix(kk),2:end-1))');

            % calculate pu term
            % dz = (y,z) > 0
            dzmat = permute(diff(squeeze(...
                runs.rgrid.z_w(:, runs.csflux.ix(kk),2:end-1)), 1, 1), [2 1]);

            % p = ∫_ζ^z ρg dz (no negative sign)
            % p = (y,z,t)
            pres = 9.81 * bsxfun(@times, rho, dzmat);
            % flip because I'm integrating from surface to (z)
            pres = flip(pres, 2);
            % integrate
            pres = cumsum(pres, 2);
            % flip back
            pres = flip(pres, 2);

            % energy flux (y, z, t, location) - mask by sponge
            peflux = bsxfun(@times, csvel .* pe, spongemask);
            keflux = bsxfun(@times, csvel .* ke + csvel .* pres, spongemask);
        end

        % check velocity plots
        debug = 0;
        if debug
            figure;
            indices = [60 100 140 170];
            for zzz = 1:length(indices)
                subplot(2,2,zzz);
                pind = indices(zzz);
                %pcolorcen( bsxfun(@times, squeeze(csvel(:,:,:,pind)), ...
                %                  offmask(:,pind))');
                pcolorcen(squeeze(csvel(:,:,:,pind))');
                center_colorbar; limc = caxis;
                hold on;
                contour(squeeze(shelfmask(:,:,:,pind))', [1 1], 'k', ...
                        'LineWidth', 2);
                %contour(squeeze(slopemask(:,:,:,pind))', [1 1], 'r', ...
                %        'LineWidth', 2);
                contour(squeeze(eddymask(:,:,:,pind))', [1 1], 'b', ...
                        'LineWidth', 2);
                caxis(limc);
                title(['Day ' num2str(time(pind)/86400)]);
            end
        end

        % Locate source of slope water being transported
        % i.e., save depth integrated transport = fn{loc}(time,bin)
        tic;
        disp('Binning horizontally...');
        if (bathyax == 2) && (sgntamp == -1)
            bins = loc(kk):3000:max(runs.rgrid.y_rho(:));
        else
            bins = flip(loc(kk):-3000:1000);
        end
        binmat = repmat(bins, [tinf 1]);
        runs.csflux.off.bins{kk} = bins;
        for mmm = 1:length(bins)-1
            % calculate transport for each water mass bin
            binmask = (csdye < bins(mmm + 1)) & (csdye >= bins(mmm));
            bintransoff = squeeze(csvel .* binmask .* csvoffmask);
            bintranson = squeeze(csvel .* binmask .* csvonmask);
            runs.csflux.off.slopewater.trans{kk}(t0:tinf, mmm) ...
                = squeeze(trapz(zvec, nansum(bsxfun(@times, ...
                                                    bsxfun(@times, bintransoff, ...
                                                           permute(offmask, [1 3 2])), dx),1),2));
            runs.csflux.on.slopewater.trans{kk}(t0:tinf, mmm) ...
                = squeeze(trapz(zvec, nansum(bsxfun(@times, ...
                                                    bsxfun(@times, bintranson, ...
                                                           permute(onmask, [1 3 2])), dx),1),2));

            % fn (z,t,bin)
            runs.csflux.off.slopewater.vtrans{kk}(:,t0:tinf, mmm) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, bintransoff, ...
                                                     permute(offmask, [1 3 2])), dx),1));
            runs.csflux.on.slopewater.vtrans{kk}(:,t0:tinf, mmm) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, bintranson, ...
                                                     permute(onmask, [1 3 2])), dx),1));
        end

        % save envelope for across-shelfbreak only = f(time,loc)
        % integrate trans over bins
        if bathyax == 2 && sgntamp == -1
            runs.csflux.off.slopewater.envelope(:,kk) = ...
                nanmax(avg1(binmat,2) .* fillnan( ...
                    runs.csflux.off.slopewater.trans{kk} > 0, 0), [], 2);
            runs.csflux.on.slopewater.envelope(:,kk) = ...
                nanmax(avg1(binmat,2) .* fillnan( ...
                    runs.csflux.on.slopewater.trans{kk} > 0, 0), [], 2);
        else
            runs.csflux.off.slopewater.envelope(:,kk) = ...
                nanmin(avg1(binmat,2) .* fillnan( ...
                    runs.csflux.off.slopewater.trans{kk} > 0, 0), [], 2);
            runs.csflux.on.slopewater.envelope(:,kk) = ...
                nanmin(avg1(binmat,2) .* fillnan( ...
                    runs.csflux.on.slopewater.trans{kk} > 0, 0), [], 2);
        end

        % integrated transport = fn(y-origin bins, loc)
        % integrate trans over time
        runs.csflux.off.slopewater.itrans{kk} = squeeze(nansum( ...
            bsxfun(@times, runs.csflux.off.slopewater.trans{kk}, ...
                   dt'), 1));
        runs.csflux.on.slopewater.itrans{kk} = squeeze(nansum( ...
            bsxfun(@times, runs.csflux.on.slopewater.trans{kk}, dt'), 1));
        toc;

        for ll=1:kk
            disp(['Parititioning slope water: ' num2str(ll) '/' num2str(kk)]);
            if (bathyax == 2) && (sgntamp == -1)
                % cyclones moving northwards
                slopemask = (csdye > loc(ll));
            else
                % everything else
                slopemask = (csdye < loc(ll));
            end

            % transports - integrate in z - f(x,t)
            disp('Integrating in z...');
            runs.csflux.slopext(:,:,kk,ll) = squeeze(trapz(zvec, slopemask .* csvel,2));
            runs.csflux.off.slopext(:,:,kk,ll) = squeeze(trapz(zvec, ...
                                                              slopemask .* csvel .* csvoffmask,2));
            runs.csflux.on.slopext(:,:,kk,ll) = squeeze(trapz(zvec, ...
                                                              slopemask .* csvel .* csvonmask,2));

            if do_energy
                runs.csflux.ikefluxxt(:,:,kk,ll) = squeeze(trapz(zvec, keflux,2));
                runs.csflux.ipefluxxt(:,:,kk,ll) = squeeze(trapz(zvec, peflux,2));

                runs.csflux.eddy.ikefluxxt(:,:,kk) = squeeze(trapz(zvec, eddymask .* keflux,2));
                runs.csflux.eddy.ipefluxxt(:,:,kk) = squeeze(trapz(zvec, eddymask .* peflux,2));
            end

            if ll == kk
                temp = bsxfun(@times, squeeze(slopemask .* csvel .* csvoffmask), ...
                              permute(offmask, [1 3 2]));
                runs.csflux.off.slopewater.vmax(:,kk) = squeeze(nanmax(nanmax(abs(temp), [], 1), [], 2));
            end

            % calculate transport as fn of vertical depth - west of
            % eddy only - f(z,t)
            disp('Binning vertically...');
            runs.csflux.off.slopezt(:,:,kk,ll) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(slopemask .* csvel .* csvoffmask), ...
                                                     permute(offmask, [1 3 2])), dx),1));
            runs.csflux.on.slopezt(:,:,kk,ll) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(slopemask .* csvel .* csvonmask), ...
                                                     permute(onmask, [1 3 2])), dx),1));

            % transports as a function of z only
            runs.csflux.off.slopewater.vertitrans(:,kk,ll) = ...
                nansum(bsxfun(@times, runs.csflux.off.slopezt(:,:,kk,ll), dt), 2);
            runs.csflux.on.slopewater.vertitrans(:,kk,ll) = ...
                nansum(bsxfun(@times, runs.csflux.on.slopezt(:,:,kk,ll), dt), 2);

            % don't need on-off paritition for edge of northern sponge
            disp('Calculating fluxes...');
            if runs.csflux.ix(kk) ~= sy2
                % west of center - f(x,t)
                runs.csflux.off.slope(t0:tinf,kk,ll) = squeeze(nansum( ...
                    bsxfun(@times, runs.csflux.off.slopext(:,:,kk,ll) .* offmask, dx),1))';

                % save average flux and itrans (west of center)
                [runs.csflux.off.itrans.slope(:,kk,ll), ...
                 runs.csflux.off.avgflux.slope(kk,ll)] = ...
                    runs.integrate_flux(time, runs.csflux.off.slope(:,kk,ll));

                % east of center
                runs.csflux.on.slope(t0:tinf,kk,ll) = squeeze(nansum( ...
                    bsxfun(@times, runs.csflux.on.slopext(:,:,kk,ll) .* onmask, dx),1))';

                % save average flux and itrans (east of center)
                [runs.csflux.on.itrans.slope(:,kk,ll), ...
                 runs.csflux.on.avgflux.slope(kk,ll)] = ...
                    runs.integrate_flux(time, runs.csflux.on.slope(:,kk,ll));
            end
        end

        % process pv
        if dopv
            % both at interior RHO points
            error('CHECK INDICES OF READING');
            if bathy.ax == 2
                start = [1 runs.csflux.ix(kk) 1 t0];
                count = [Inf 1 Inf Inf];
            else
                start = [runs.csflux.ix(kk) 1 1 t0];
                count = [1 Inf Inf Inf];
                error('PV fluxes not tested for NS isobaths');
            end
            pv = ncread(vorname,'pv',start,count);
            rv = avg1(avg1(ncread(vorname,'rv',start,count+[0 1 0 0]),1),2);

            % get vorticity fluxes
            % first, depth integrated
            pvcsflux = squeeze(trapz(avg1(runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1),1), ...
                                     pv .* avg1(csvel,3),3));
            rvcsflux = squeeze(trapz(avg1(runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1),1), ...
                                     rv .* avg1(csvel,3),3));

            % now east & west of eddy center
            runs.csflux.off.pv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, pvcsflux .* offmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';
            runs.csflux.on.pv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, pvcsflux .* onmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.off.rv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, rvcsflux .* offmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';
            runs.csflux.on.rv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, rvcsflux .* onmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';
        end
    end

    % process energy for all locations
    if do_energy
        xvec = runs.rgrid.x_rho(1,2:end-1);
        runs.csflux.ikeflux = squeeze(trapz(xvec, ...
                                            runs.csflux.ikefluxxt, 1));
        runs.csflux.ipeflux = squeeze(trapz(xvec, ...
                                            runs.csflux.ipefluxxt, 1));

        runs.csflux.eddy.ikeflux = squeeze(trapz(xvec, ...
                                                 runs.csflux.eddy.ikefluxxt, 1));
        runs.csflux.eddy.ipeflux = squeeze(trapz(xvec, ...
                                                 runs.csflux.eddy.ipefluxxt, 1));
    end

    % save fluxes
    runs.csflux.flag_found = flag_found;

    hash = githash([mfilename('fullpath') '.m']);
    runs.csflux.hash = hash;

    csflux = runs.csflux;
    asflux = runs.asflux;

    if strcmpi(precision, 'single')
        runs.csflux.slopext = single(runs.csflux.slopext);
        runs.csflux.eddyxt = single(runs.csflux.eddyxt);
        runs.csflux.off.slopext = single(runs.csflux.off.slopext);
        runs.csflux.off.eddyxt = single(runs.csflux.off.eddyxt);
        runs.csflux.on.slopext = single(runs.csflux.on.slopext);
        runs.csflux.on.eddyxt = single(runs.csflux.on.eddyxt);
        runs.csflux.off.slopezt = single(runs.csflux.off.slopezt);
        runs.csflux.off.eddyzt = single(runs.csflux.off.eddyzt);
        runs.csflux.on.slopezt = single(runs.csflux.on.slopezt);
        runs.csflux.on.eddyzt = single(runs.csflux.on.eddyzt);
    end

    disp('Saving data...');
    save([runs.dir '/fluxes.mat'], 'csflux', 'asflux');
    toc(ticstart);
end