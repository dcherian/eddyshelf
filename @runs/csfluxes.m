% quantify cross-shelfbreak and along-shelfbreak fluxes of a whole
% bunch of stuff:
function [] = csfluxes(runs, ftype)

% Things I want to calculate fluxes of:
% 1. RV & PV
% 2. Shelf water dye & eddy dye
    ticstart = tic;

    % Use history or avg files?
    if ~exist('ftype', 'var') || isempty(ftype), ftype = 'his'; end

    % calculate energy fluxes?
    do_energy = 0;

    % do PV fluxes?
    dopv = 0;

    vorname = [runs.dir '/ocean_vor.nc'];
    % need some kind of initial time instant - decided by streamer mask
    % now
    runs.csflux = [];

    tstart = 1;%find(repnan(runs.streamer.time,0) ==
               %0,1,'last') + 1;
    revind = runs.eddy.trevind;

    % quantities for area-averaging if needed
    h = runs.bathy.h(2:end-1,2:end-1);

    % sort out isobaths across which to calculate transports
    Y = max(runs.rgrid.y_rho(:));

    % Non-dimensional isobath = (y_{isobath} - y_{sb})/(Eddy center
    % location)
    xsb = runs.bathy.xsb; xsl = runs.bathy.xsl;
    if runs.bathy.axis == 'y'
        [~,R,restind] = runs.locate_resistance;
    else
        [R,~,restind] = runs.locate_resistance;
    end
    if isempty(R)
        error(['locate_resistance didn''t return anything for ' runs.name]);
    end

    R = R - xsb;
    runs.csflux.ndloc = linspace(0,2,13);
    runs.csflux.R = R;
    loc = runs.csflux.ndloc* R + xsb;
    runs.csflux.ndloc(loc > xsl) = [];
    loc(loc > xsl) = [];
    runs.animate_field('zeta',[],restind,1);
    if runs.bathy.axis == 'y'
        liney(loc/1000, [], 'r');
    else
        linex(loc/1000, [], 'r');
    end

    if runs.params.bathy.axis == 'x'
        csvelid = 'u';
        asvelid = 'v';
        bathyax = 1;
        indices = vecfind(runs.eddy.xr(:,1),loc);

        % append sponge edge to existing values
        sy2 = runs.spng.sx2;
        % indices = [indices sy2];
        % loc = [loc runs.eddy.xr(sy2, 1)];
    else
        csvelid = 'v';
        asvelid = 'u';
        bathyax = 2;
        indices = vecfind(runs.eddy.yr(1,:),loc);

        % append sponge edge to existing values
        sy2 = runs.spng.sy2;
        % indices = [indices sy2];
        % loc = [loc runs.eddy.yr(1, sy2)];
    end

    sz = size(runs.sponge);

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

    % how much of the time vector should I read?
    if isfield(runs.eddy, 'tend')
        tinf = runs.eddy.tend;
    else
        tinf = length(runs.time);
    end

    time = runs.time(1:tinf);
    t0 = 1;

    % interpolate center locations
    % Don't think this is required anymore. I don't use average
    % files at all.
    % if strcmpi(ftype, 'his')
    %     time = dc_roms_read_data(runs.dir, 'ocean_time', [1 tinf], {}, [], ...
    %                              [], 'his');
    %     if length(time) ~= length(runs.eddy.t(1:runs.eddy.tend))
    %         t0 = find_approx(time, runs.time(tstart), 1);
    %         if bathyax == 2
    %             cxi = interp1(runs.eddy.t(tstart:end)*86400, runs.eddy.vor.ee(tstart:end), ...
    %                           time(t0:end));
    %         else
    %             cxi = interp1(runs.eddy.t(tstart:end)*86400, runs.eddy.vor.se(tstart:end), ...
    %                           time(t0:end));
    %         end
    %     else
    %         t0 = tstart;
    %         if bathyax == 2
    %             cxi = runs.eddy.vor.ee(tstart:end);
    %         else
    %             cxi = runs.eddy.vor.se(tstart:end);
    %         end
    %     end
    % else
    %     if strcmpi(ftype, 'avg')
    %         t0 = tstart;
    %         time = runs.time;
    %         if bathyax == 2
    %             cxi = runs.eddy.vor.ee(t0:end);
    %         else
    %             cxi = runs.eddy.vor.se(t0:end);
    %         end
    %     end
    % end

    % size for initialization
    nloc = length(loc);
    szfull = size(runs.bathy.h);
    szflux = [tinf nloc nloc];
    if bathyax == 1
        szfluxxt = [szfull(2)-2 tinf nloc];
    else
        szfluxxt = [szfull(1)-2 tinf nloc];
    end

    % initialize - mass flux variables
    runs.csflux.west.slope = nan(szflux);
    runs.csflux.west.eddy = nan(szflux);

    runs.csflux.west.itrans.slope = nan(szflux);
    runs.csflux.west.itrans.eddy = nan(szflux);

    runs.csflux.east.slope = nan(szflux);
    runs.csflux.east.eddy = nan(szflux);

    runs.csflux.east.itrans.slope = nan(szflux);
    runs.csflux.east.itrans.eddy = nan(szflux);

    runs.csflux.west.bins = cell([1 nloc]);

    runs.csflux.west.slopewater.vmax = nan([tinf nloc]);

    runs.csflux.west.slopewater.trans = cell([1 nloc nloc]);
    runs.csflux.east.slopewater.vtrans = cell([1 nloc nloc]);
    runs.csflux.west.slopewater.itrans = cell([1 nloc nloc]);
    runs.csflux.west.slopewater.envelope = nan([tinf nloc]);

    runs.csflux.east.slopewater.trans = cell([1 nloc nloc]);
    runs.csflux.east.slopewater.vtrans = cell([1 nloc nloc]);
    runs.csflux.east.slopewater.itrans = cell([1 nloc nloc]);
    runs.csflux.east.slopewater.envelope = nan([tinf nloc]);

    runs.csflux.slopext  = nan(szfluxxt);
    runs.csflux.eddyxt   = nan(szfluxxt);
    runs.csflux.westmask = nan(szfluxxt);
    runs.csflux.eastmask = nan(szfluxxt);

    runs.csflux.west.slopezt = nan([runs.rgrid.N tinf nloc nloc]);
    runs.csflux.west.eddyzt = nan([runs.rgrid.N tinf nloc nloc]);
    runs.csflux.east.slopezt = nan([runs.rgrid.N tinf nloc nloc]);
    runs.csflux.east.eddyzt = nan([runs.rgrid.N tinf nloc nloc]);

    runs.csflux.west.slopewater.vertitrans = nan([runs.rgrid.N nloc]);
    runs.csflux.west.eddywater.vertitrans = nan([runs.rgrid.N nloc]);
    runs.csflux.east.slopewater.vertitrans = nan([runs.rgrid.N nloc]);
    runs.csflux.east.eddywater.vertitrans = nan([runs.rgrid.N nloc]);

    % initialize - pv fluxes
    if exist(vorname, 'file') && dopv == 1
        pvtime = ncread(vorname, 'ocean_time');
        if isequal(pvtime', time) || isequal(pvtime, time)
            dopv = 1;
            runs.csflux.west.pv = nan(szflux);
            runs.csflux.west.rv = nan(szflux);
            runs.csflux.east.pv = nan(szflux);
            runs.csflux.east.rv = nan(szflux);
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

    % sponge mask
    if runs.bathy.axis == 'y'
        spongemask = ~runs.sponge(2:end-1,1);
    else
        spongemask = ~runs.sponge(1,2:end-1)';
    end

    % for integrated transport diagnostics
    dt = [time(2)-time(1) diff(time)];
    if bathyax == 2
        dx = 1./runs.rgrid.pm(1,2:end-1)';
    else
        dx = 1./runs.rgrid.pn(2:end-1,1);
    end

    % loop over all isobaths
    for kk=1:length(indices)
        disp(['Doing isobath ' num2str(kk) '/', ...
              num2str(length(indices))]);

        if kk == 1
            % cxi here = eastern edge because export occurs west of the eastern edge
            % dimensions - (along-shore) x time
            % The problem is that this is sensitive to contour
            % detection.
            if bathyax == 2
                cxi = runs.eddy.vor.ee(tstart:end);
            else
                cxi = runs.eddy.vor.se(tstart:end);
            end
        else
            if bathyax == 2
                cxi = runs.eddy.mx(tstart:end);
            else
                cxi = runs.eddy.my(tstart:end);
            end
        end

        % east and west (w.r.t cxi) masks
        if bathyax == 2
            westmask = bsxfun(@times, ...
                              bsxfun(@lt, runs.eddy.xr(:,1), cxi(t0:tinf)), ...
                              spongemask);
        else
            westmask = bsxfun(@times, ...
                              bsxfun(@gt, runs.eddy.yr(1,:)', cxi(t0:tinf)), ...
                              spongemask);
        end
        eastmask = bsxfun(@times, 1 - westmask, spongemask);

        runs.csflux.westmask(:,:,kk) = westmask;
        runs.csflux.eastmask(:,:,kk) = eastmask;

        % read along-shore section of cross-shore vel.
        % dimensions = (x/y , z , t )
        % average carefully to get values at RHO point
        csvel = squeeze(avg1(dc_roms_read_data(runs.dir, csvelid, ...
                                       [t0 tinf], {runs.bathy.axis runs.csflux.ix(kk)-1 ...
                            runs.csflux.ix(kk)}, [], runs.rgrid, ftype, ...
                                       'single'), bathyax));
        csvel = csvel(2:end-1,:,:,:);

        % process cross-shelf dye
        csdye = dc_roms_read_data(runs.dir, runs.csdname, ...
                                  [t0 tinf], {runs.bathy.axis runs.csflux.ix(kk)+1 ...
                            runs.csflux.ix(kk)+1}, [], runs.rgrid, ...
                                  ftype, 'single');
        csdye = csdye(2:end-1,:,:);

        % read eddye
        eddye = dc_roms_read_data(runs.dir, runs.eddname, ...
                                  [t0 tinf], {runs.bathy.axis runs.csflux.ix(kk)+1 ...
                            runs.csflux.ix(kk)+1}, [], runs.rgrid, ...
                                  ftype, 'single');
        eddye = eddye(2:end-1,:,:);

        % define eddy water mass
        eddymask = eddye > runs.eddy_thresh;

        if bathyax == 2
            zvec = runs.rgrid.z_r(:, runs.csflux.ix(kk)+1, 1);
        else
            zvec = runs.rgrid.z_r(:, 1, runs.csflux.ix(kk)+1);
        end
        runs.csflux.vertbins(:,kk) = ...
            runs.rgrid.z_r(:, runs.csflux.ix(kk)+1, 1);

        runs.csflux.eddyxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                   eddymask .* csvel,2));
        runs.csflux.west.eddyzt(:,:,kk) = ...
            squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(eddymask .* csvel), ...
                                                 permute(westmask, ...
                                                         [1 3 2])), dx),1));
        runs.csflux.east.eddyzt(:,:,kk) = ...
            squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(eddymask .* csvel), ...
                                                 permute(eastmask, ...
                                                         [1 3 2])), dx),1));
        % transports as a function of z only
        runs.csflux.west.eddywater.vertitrans(:,kk) = ...
            nansum(bsxfun(@times, runs.csflux.west.eddyzt(:,:,kk), dt), 2);
        runs.csflux.east.eddywater.vertitrans(:,kk) = ...
            nansum(bsxfun(@times, runs.csflux.east.eddyzt(:,:,kk), dt), 2);

        % calculate fluxes
        if runs.csflux.ix(kk) ~= sy2
            runs.csflux.west.eddy(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.eddyxt(:,:,kk) .* westmask, ...
                       dx),1))';

            runs.csflux.east.eddy(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.eddyxt(:,:,kk) .* eastmask, ...
                       dx),1))';

            [runs.csflux.west.itrans.eddy(:,kk), ...
             runs.csflux.west.avgflux.eddy(kk)] = ...
                runs.integrate_flux(time, runs.csflux.west.eddy(:,kk));

            [runs.csflux.east.itrans.eddy(:,kk), ...
             runs.csflux.east.avgflux.eddy(kk)] = ...
                runs.integrate_flux(time, runs.csflux.east.eddy(:,kk));
        end

        if do_energy
            disp('Calculating energy...');
            % read velocity and calculate pressure etc.
            asvel = avg1(dc_roms_read_data(runs.dir, asvelid, [t0 tinf], ...
                                      {runs.bathy.axis runs.csflux.ix(kk) ...
                                runs.csflux.ix(kk)}, ...
                                      [], runs.rgrid, ftype, ...
                                      'single'), 1);

            % read velocity and calculate pressure etc.
            rho = dc_roms_read_data(runs.dir, 'rho', [t0 tinf], ...
                                      {runs.bathy.axis runs.csflux.ix(kk) ...
                                runs.csflux.ix(kk)}, ...
                                      [], runs.rgrid, ftype, ...
                                      'single');
            rho = rho(2:end-1,:,:);

            ke = 1/2 * rho .* (csvel.^2 + asvel.^2);

            pe = - runs.params.phys.g * bsxfun(@times, rho, ...
                                               squeeze(runs.rgrid ...
                                                       .z_r(:,runs.csflux.ix(kk),2:end-1))');

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
                %                  westmask(:,pind))');
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
        bins = flip(loc(kk):-3000:1000);
        binmat = repmat(bins, [tinf 1]);
        runs.csflux.west.bins{kk} = bins;
        for mmm = 1:length(bins)-1
            % calculate transport for each water mass bin
            binmask = (csdye < bins(mmm + 1)) & (csdye >= bins(mmm));
            bintrans = squeeze(csvel .* binmask);
            runs.csflux.west.slopewater.trans{kk}(t0:tinf, mmm) ...
                = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1), ...
                    nansum(bsxfun(@times, bsxfun(@times, bintrans, ...
                                                 permute(westmask, ...
                                                         [1 3 2])), dx),1),2));
            runs.csflux.east.slopewater.trans{kk}(t0:tinf, mmm) ...
                = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1), ...
                    nansum(bsxfun(@times, bsxfun(@times, bintrans, ...
                                                 permute(eastmask, ...
                                                         [1 3 2])), dx),1),2));

            % fn (z,t,bin)
            runs.csflux.west.slopewater.vtrans{kk}(:,t0:tinf, mmm) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, bintrans, ...
                                                     permute(westmask, ...
                                                             [1 3 2])), dx),1));
            runs.csflux.east.slopewater.vtrans{kk}(:,t0:tinf, mmm) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, bintrans, ...
                                                     permute(eastmask, ...
                                                             [1 3 2])), dx),1));
        end

        % save envelope for across-shelfbreak only = f(time,loc)
        % integrate trans over bins
        runs.csflux.west.slopewater.envelope(:,kk) = ...
            nanmin(avg1(binmat,2) .* fillnan( ...
                runs.csflux.west.slopewater.trans{kk} > 0, 0), [], 2);
        runs.csflux.east.slopewater.envelope(:,kk) = ...
            nanmin(avg1(binmat,2) .* fillnan( ...
                runs.csflux.east.slopewater.trans{kk} > 0, 0), [], 2);

        % integrated transport = fn(y-origin bins, loc)
        % integrate trans over time
        runs.csflux.west.slopewater.itrans{kk} = squeeze(nansum( ...
            bsxfun(@times, runs.csflux.west.slopewater.trans{kk}, ...
                   dt'), 1));
        runs.csflux.east.slopewater.itrans{kk} = squeeze(nansum( ...
            bsxfun(@times, runs.csflux.east.slopewater.trans{kk}, dt'), 1));
        toc;

        for ll=1:kk
            disp(['Parititioning slope water: ' num2str(ll) '/' num2str(kk)]);
            slopemask = (csdye < loc(ll));
            % transports - integrate in z - f(x,t)
            disp('Integrating in z...');
            runs.csflux.slopext(:,:,kk,ll) = squeeze(trapz(zvec, ...
                                                        slopemask .* csvel,2));

            if do_energy
                runs.csflux.ikefluxxt(:,:,kk,ll) = squeeze(trapz(zvec, ...
                                                              keflux,2));
                runs.csflux.ipefluxxt(:,:,kk,ll) = squeeze(trapz(zvec, ...
                                                              peflux,2));

                runs.csflux.eddy.ikefluxxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                                  eddymask ...
                                                                  .* keflux,2));
                runs.csflux.eddy.ipefluxxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                                  eddymask .* ...
                                                                  peflux,2));
            end

            if ll == kk
                temp = bsxfun(@times, squeeze(slopemask .* csvel), ...
                              permute(westmask, [1 3 2]));
                runs.csflux.west.slopewater.vmax(:,kk) = ...
                    squeeze(nanmax(nanmax(abs(temp), [], 1), [], 2));
            end

            % calculate transport as fn of vertical depth - west of
            % eddy only - f(z,t)
            disp('Binning vertically...');
            runs.csflux.west.slopezt(:,:,kk,ll) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(slopemask .* csvel), ...
                                                     permute(westmask, [1 3 2])), dx),1));
            runs.csflux.east.slopezt(:,:,kk,ll) = ...
                squeeze(nansum(bsxfun(@times, bsxfun(@times, squeeze(slopemask .* csvel), ...
                                                     permute(eastmask, [1 3 2])), dx),1));

            % transports as a function of z only
            runs.csflux.west.slopewater.vertitrans(:,kk,ll) = ...
                nansum(bsxfun(@times, runs.csflux.west.slopezt(:,:,kk,ll), dt), 2);
            runs.csflux.east.slopewater.vertitrans(:,kk,ll) = ...
                nansum(bsxfun(@times, runs.csflux.east.slopezt(:,:,kk,ll), dt), 2);

            % don't need east-west paritition for edge of northern sponge
            disp('Calculating fluxes...');
            if runs.csflux.ix(kk) ~= sy2
                % west of center - f(x,t)
                runs.csflux.west.slope(t0:tinf,kk,ll) = squeeze(nansum( ...
                    bsxfun(@times, runs.csflux.slopext(:,:,kk,ll) .* westmask, ...
                           dx),1))';

                % save average flux and itrans (west of center)
                [runs.csflux.west.itrans.slope(:,kk,ll), ...
                 runs.csflux.west.avgflux.slope(kk,ll)] = ...
                    runs.integrate_flux(time, runs.csflux.west.slope(:,kk,ll));

                % east of center
                runs.csflux.east.slope(t0:tinf,kk,ll) = squeeze(nansum( ...
                    bsxfun(@times, runs.csflux.slopext(:,:,kk,ll) .* eastmask, ...
                           dx),1))';

                % save average flux and itrans (east of center)
                [runs.csflux.east.itrans.slope(:,kk,ll), ...
                 runs.csflux.east.avgflux.slope(kk,ll)] = ...
                    runs.integrate_flux(time, runs.csflux.east.slope(:,kk,ll));
            end
        end

        % process pv
        if dopv
            % both at interior RHO points
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
            runs.csflux.west.pv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, pvcsflux .* westmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';
            runs.csflux.east.pv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, pvcsflux .* eastmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.west.rv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, rvcsflux .* westmask, ...
                       1./runs.rgrid.pm(1,2:end-1)'),1))';
            runs.csflux.east.rv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, rvcsflux .* eastmask, ...
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
    runs.csflux.time = time;

    hash = githash([mfilename('fullpath') '.m']);
    runs.csflux.hash = hash;

    csflux = runs.csflux;
    asflux = runs.asflux;

    disp('Saving data...');
    save([runs.dir '/fluxes.mat'], 'csflux', 'asflux');
    toc(ticstart);
end