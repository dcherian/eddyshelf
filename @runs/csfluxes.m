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
    do_energy = 1;

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
    Las = max(runs.rgrid.x_rho(:));
    %Lcs = runs.bathy.xsb;
    h = runs.bathy.h(2:end-1,2:end-1);

    % find sponge edges
    sz = size(runs.sponge);
    sy2 = find(runs.sponge(sz(1)/2, :) == 1, 1, 'first') - 1;

    % sort out isobaths across which to calculate transports
    loc = runs.bathy.xsb; %linspace(runs.bathy.xsb, runs.bathy.xsl, 4);
    if runs.params.bathy.axis == 'x'
        csvelid = 'u';
        asvelid = 'v';
        bathyax = 1;
        indices = vecfind(runs.eddy.xr(:,1),loc);
    else
        csvelid = 'v';
        asvelid = 'u';
        bathyax = 2;
        indices = vecfind(runs.eddy.yr(1,:),loc);
        %runs.rgrid.y_rho(vecfind(runs.bathy.h(1,:),[250 1000]),1)']);
    end

    % append sponge edge to existing values
    indices = [indices sy2];
    loc = [loc runs.eddy.yr(1, sy2)];

    % save locations
    runs.csflux.x = loc;
    % save indices for locations - w.r.t INTERIOR RHO POINTS
    runs.csflux.ix = indices;
    % save isobath values
    runs.csflux.h = ceil(runs.bathy.h(1,runs.csflux.ix));

    runs.csflux.comment = ['shelf / slope / eddy= which water mass am I ', ...
                        ' targeting? |\n (x,ix,h) = (location, index, depth) ' ...
                        'at which I''m calculating transport |\n '];

    % how much of the time vector should I read?
    if isfield(runs.eddy, 'tend')
        tinf = runs.eddy.tend;
    else
        tinf = Inf;
    end

    % interpolate center locations
    if strcmpi(ftype, 'his')
        time = dc_roms_read_data(runs.dir, 'ocean_time', [1 tinf], {}, [], ...
                                 [], 'his');
        if length(time) ~= length(runs.eddy.t)
            t0 = find_approx(time, runs.time(tstart), 1);
            if bathyax == 2
                cxi = interp1(runs.eddy.t(tstart:end)*86400, runs.eddy.vor.ee(tstart:end), ...
                              time(t0:end));
            else
                cxi = interp1(runs.eddy.t(tstart:end)*86400, runs.eddy.vor.se(tstart:end), ...
                              time(t0:end));
            end
        else
            t0 = tstart;
            if bathyax == 2
                cxi = runs.eddy.vor.ee(tstart:end);
            else
                cxi = runs.eddy.vor.se(tstart:end);
            end
        end
    else
        if strcmpi(ftype, 'avg')
            t0 = tstart;
            time = runs.time;
            if bathyax == 2
                cxi = runs.eddy.vor.ee(t0:end);
            else
                cxi = runs.eddy.vor.se(t0:end);
            end
        end
    end

    if isinf(tinf)
        tinf = length(time);
    end

    % size for initialization
    szfull = size(runs.bathy.h)
    szflux = [tinf length(loc)];
    szfluxxt = [szfull(1)-2 tinf length(loc)];

    % initialize - mass flux variables
    runs.csflux.west.shelf = nan(szflux);
    runs.csflux.west.slope = nan(szflux);
    runs.csflux.west.eddy = nan(szflux);

    runs.csflux.west.itrans.shelf = nan(szflux);
    runs.csflux.west.itrans.slope = nan(szflux);
    runs.csflux.west.itrans.eddy = nan(szflux);

    runs.csflux.east.shelf = nan(szflux);
    runs.csflux.east.slope = nan(szflux);
    runs.csflux.east.eddy = nan(szflux);

    runs.csflux.east.itrans.shelf = nan(szflux);
    runs.csflux.east.itrans.slope = nan(szflux);
    runs.csflux.east.itrans.eddy = nan(szflux);

    rr = runs.rrshelf;
    maxrr = ceil(runs.bathy.xsb/rr);
    runs.csflux.west.shelfwater.bins = (1:maxrr) * rr;
    runs.csflux.west.shelfwater.trans = nan([szflux maxrr]);
    binmat = repmat(runs.csflux.west.shelfwater.bins, [tinf 1]);

    runs.csflux.west.shelfwater.vertitrans = nan([runs.rgrid.N sz(2)]);

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

    % east and west (w.r.t eddy center) masks
    % cxi here = eastern edge because export occurs west of the eastern edge
    % dimensions - (along-shore) x time
    if bathyax == 2
        westmask = bsxfun(@times, ...
                          bsxfun(@lt, runs.eddy.xr(:,1), cxi(t0:tinf)), ...
                          spongemask);
    else
        westmask = bsxfun(@times, ...
                          bsxfun(@gt, runs.eddy.yr(1,:)', cxi(t0:tinf)), ...
                          spongemask);
    end
    eastmask = bsxfun(@times, 1 - westmask, ...
                      spongemask);

    % for integrated transport diagnostics
    dt = [time(2)-time(1) diff(time)];

    % loop over all isobaths
    for kk=1:length(indices)
        disp(['Doing isobath ' num2str(kk) '/', ...
              num2str(length(indices))]);

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

        % define water masses
        shelfmask = (csdye < runs.bathy.xsb);
        slopemask = (csdye >= runs.bathy.xsb) & ...
            (csdye <= runs.bathy.xsl);
        eddymask = eddye > runs.eddy_thresh;

        if do_energy
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

        % transports - integrate in z - f(x,t)
        if bathyax == 2
            zvec = runs.rgrid.z_r(:, runs.csflux.ix(kk)+1, 1);
        else
            zvec = runs.rgrid.z_r(:, 1, runs.csflux.ix(kk)+1);
        end
        runs.csflux.shelfxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                    shelfmask .* csvel,2));
        runs.csflux.slopext(:,:,kk) = squeeze(trapz(zvec, ...
                                                    slopemask .* csvel,2));
        runs.csflux.eddyxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                   eddymask .* csvel,2));

        if do_energy
            runs.csflux.ikefluxxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                        keflux,2));
            runs.csflux.ipefluxxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                          peflux,2));

            runs.csflux.eddy.ikefluxxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                          eddymask ...
                                                          .* keflux,2));
            runs.csflux.eddy.ipefluxxt(:,:,kk) = squeeze(trapz(zvec, ...
                                                        eddymask .* ...
                                                          peflux,2));
        end

        %            for ntt = 60:size(runs.csflux.shelfxt,2)
        %    vec = runs.csflux.shelfxt(:,ntt,kk);
        %    clf;
        %    plot(vec);
        %    linex(find(vec < 0, 1, 'first'));
        %    pause(1);
        %end

        % calculate transport as fn of vertical depth - west of
        % eddy only - f(z,t)
        if bathyax == 2
            dx = 1./runs.rgrid.pm(1,2:end-1)';
        else
            dx = 1./runs.rgrid.pn(2:end-1,1);
        end
        runs.csflux.shelfzt = squeeze(nansum(bsxfun(@times, ...
                                                    bsxfun(@times, squeeze(shelfmask .* csvel), ...
                                                          permute(westmask, [1 3 2])), dx),1));
        runs.csflux.west.shelfwater.vertitrans(:,kk) = nansum(bsxfun(@times, ...
                                                          runs.csflux.shelfzt, dt), 2);
        runs.csflux.west.shelfwater.vertbins(:,kk) = ...
            runs.rgrid.z_r(:, runs.csflux.ix(kk)+1, ...
                           1);

        % water mass analysis of fluxes for shelfbreak only
        tic;
        if kk == 1
            for mmm = 1:maxrr-1
                % calculate transport for each shelf water mass
                % bin - extra shelfmask multiplication is just a
                % check. adds 2 seconds.
                binmask = (csdye < runs.csflux.west.shelfwater.bins(mmm + 1)) ...
                          & (csdye >= runs.csflux.west.shelfwater.bins(mmm));
                bintrans = squeeze(csvel .* binmask);
                runs.csflux.west.shelfwater.trans(t0:tinf, kk, mmm) ...
                    = squeeze(trapz(runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1), ...
                                    nansum(bsxfun(@times, ...
                                                  bsxfun(@times, ...
                                                         bintrans,permute(westmask, [1 3 2])), ...
                                                  dx),1),2));
            end

            % save envelope for across-shelfbreak only
            runs.csflux.west.shelfwater.envelope = nanmin(binmat .* ...
                                                          fillnan(squeeze( ...
                                                              runs.csflux.west.shelfwater.trans(:,1,:)) ...
                                                              > ...
                                                              0, 0), [], 2);

            runs.csflux.west.shelfwater.itrans = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.west.shelfwater.trans(:,1,:), ...
                       dt'), 1));
        end
        toc;

        % don't need east-west paritition for edge of northern sponge
        if runs.csflux.ix(kk) ~= sy2
            % west of center - f(x,t)
            runs.csflux.west.shelf(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.shelfxt(:,:,kk) .* westmask, ...
                       dx),1))';

            runs.csflux.west.slope(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.slopext(:,:,kk) .* westmask, ...
                       dx),1))';

            runs.csflux.west.eddy(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.eddyxt(:,:,kk) .* westmask, ...
                       dx),1))';

            % save average flux and itrans (west of center)
            [runs.csflux.west.itrans.shelf(:,kk), ...
             runs.csflux.west.avgflux.shelf(kk)] = ...
                runs.integrate_flux(time, runs.csflux.west.shelf(:,kk));

            [runs.csflux.west.itrans.slope(:,kk), ...
             runs.csflux.west.avgflux.slope(kk)] = ...
                runs.integrate_flux(time, runs.csflux.west.slope(:,kk));

            [runs.csflux.west.itrans.eddy(:,kk), ...
             runs.csflux.west.avgflux.eddy(kk)] = ...
                runs.integrate_flux(time, runs.csflux.west.eddy(:,kk));

            % east of center
            runs.csflux.east.shelf(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.shelfxt(:,:,kk) .* eastmask, ...
                       dx),1))';

            runs.csflux.east.slope(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.slopext(:,:,kk) .* eastmask, ...
                       dx),1))';

            runs.csflux.east.eddy(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, runs.csflux.eddyxt(:,:,kk) .* eastmask, ...
                       dx),1))';

            % save average flux and itrans (east of center)
            [runs.csflux.east.itrans.shelf(:,kk), ...
             runs.csflux.east.avgflux.shelf(kk)] = ...
                runs.integrate_flux(time, runs.csflux.east.shelf(:,kk));

            [runs.csflux.east.itrans.slope(:,kk), ...
             runs.csflux.east.avgflux.slope(kk)] = ...
                runs.integrate_flux(time, runs.csflux.east.slope(:,kk));

            [runs.csflux.east.itrans.eddy(:,kk), ...
             runs.csflux.east.avgflux.eddy(kk)] = ...
                runs.integrate_flux(time, runs.csflux.east.eddy(:,kk));
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
    toc(ticstart);

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
    %runs.asflux.time = time;

    runs.csflux.westmask = westmask;
    runs.csflux.eastmask = eastmask;

    hash = githash([mfilename('fullpath') '.m']);
    runs.csflux.hash = hash;

    csflux = runs.csflux;
    asflux = runs.asflux;

    save([runs.dir '/fluxes.mat'], 'csflux', 'asflux');
end