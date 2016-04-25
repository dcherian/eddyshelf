% calculates and saves average velocity seen by shelf-slope water at a given isobath
function [] = avgStreamerVelSection(runs, iso)

    ticstart = tic;
    debug = 0;

    disp('===================================')
    disp([runs.name '.avgStreamerVelSection']);
    disp('===================================')

    if ~exist('iso', 'var') | isempty(iso)
        iso = [1:niso];
    end

    xivec = -200:runs.rgrid.dx/1000:200;
    xl = length(xivec);
    N = runs.rgrid.N;
    niso = length(runs.csflux.x);

    szmeanvel = [xl N niso];
    vmean = nan(szmeanvel);
    offvmean = nan(szmeanvel);
    onvmean = nan(szmeanvel);

    for isobath = iso
        source = isobath;
        [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:, isobath, source));
        tindices = [start stop];
        ix = runs.csflux.ix(isobath); % WORKAROUND FOR CSFLUXES BUG
        volr = {runs.bathy.axis ix ix};
        volv = {runs.bathy.axis ix-1 ix};
        zvec = runs.csflux.vertbins(:,isobath);

        if runs.bathy.axis == 'x'
            xvec = runs.rgrid.y_rho(2:end-1,1)/1000;
            mx = runs.eddy.my(start:stop)/1000;
            bathyax = 1;
        else
            xvec = runs.rgrid.x_rho(1,2:end-1)/1000;
            mx = runs.eddy.mx(start:stop)/1000;
            bathyax = 2;
        end

        disp([runs.name ' | isobath = ' num2str(isobath) ' | Reading data...'])
        v = squeeze(avg1(dc_roms_read_data(runs, runs.csvelname, tindices, volv), bathyax));
        csdye = dc_roms_read_data(runs, runs.csdname, tindices, volr);

        % interpolate onto common grid centered on eddy center
        % this is not working well!!! - unsure why
        intstart = tic;
        disp([runs.name ' | isobath = ' num2str(isobath) ' | Interpolating data...'])
        clear vi csdyei
        for tt=1:size(v,3)
            for zz=1:size(v,2)
                vi(:,zz,tt) = interp1(xvec-mx(tt), v(2:end-1,zz,tt), xivec);
                csdyei(:,zz,tt) = interp1(xvec-mx(tt), csdye(2:end-1,zz,tt), xivec);
            end
        end
        toc(intstart);

        % offshore transport mask
        if isobath ~= 1
            if runs.bathy.axis == 'y'
                offmask = xivec' < 0;
            else
                offmask = xivec' > 0;
            end
        else
            offmask = 1; %xivec' < runs.eddy.rhovor.dia(1)/2;
        end

        offvmean(:,:,isobath) = mean(bsxfun(@times, vi .* (csdyei < runs.csflux.x(source)), ...
                                            offmask), 3);
        onvmean(:,:,isobath) = mean(bsxfun(@times, vi .* (csdyei < runs.csflux.x(source)), ...
                                            ~offmask), 3);
        vmean(:,:,isobath) = mean(vi .* (csdyei < runs.csflux.x(source)), 3);
    end

    runs.streamer = [];
    runs.streamer.xivec = xivec';
    runs.streamer.zvec = runs.csflux.vertbins;

    % (x,z) fields
    runs.streamer.vmean = vmean;
    runs.streamer.off.vmean = offvmean;
    runs.streamer.on.vmean = onvmean;

    % integrated profiles
    runs.streamer.off.zprof = squeeze(trapz(xivec, repnan(offvmean,0), 1)); % offshore mean profile
    runs.streamer.on.zprof = squeeze(trapz(xivec, repnan(onvmean,0), 1)); % onshore mean profile
    runs.streamer.off.xprof = squeeze(trapz(zvec, offvmean, 2));% offshore vel : x-profile
    runs.streamer.on.xprof = squeeze(trapz(zvec, onvmean, 2));% onshore vel : x-profile
    runs.streamer.xprof = squeeze(trapz(zvec, vmean, 2)); % full vel : x-profile

    % diagnostics
    [xmax,xind] = nanmax(runs.streamer.off.xprof, [], 1);
    runs.streamer.off.ixpeak = xind;
    runs.streamer.off.xpeak = xivec(xind);
    [zmax,zind] = nanmax(runs.streamer.off.zprof, [], 1);
    runs.streamer.off.izpeak = zind;
    runs.streamer.off.zpeak = zvec(zind);

    [xmax,xind] = nanmax(runs.streamer.on.xprof, [], 1);
    runs.streamer.on.ixpeak = xind;
    runs.streamer.on.xpeak = xivec(xind);
    [zmax,zind] = nanmax(runs.streamer.on.zprof, [], 1);
    runs.streamer.on.izpeak = zind;
    runs.streamer.on.zpeak = zvec(zind);

    for iso=1:niso
        runs.streamer.off.xwidth(iso) = ...
            findProfilePeakWidth(runs.streamer.off.xprof(:,iso), xivec);
        runs.streamer.on.xwidth(iso) = ...
            findProfilePeakWidth(runs.streamer.on.xprof(:,iso), xivec);

        zvec = runs.csflux.vertbins(:,iso);
        [~,ilo,~] = findProfilePeakWidth(runs.streamer.off.zprof(:,iso), zvec);
        if ~isempty(ilo) & ~isnan(ilo)
            runs.streamer.off.zwidth(iso) = zvec(ilo);
        else
            runs.streamer.off.zwidth(iso) = NaN;
        end
        %[~,ilo,~] = findProfilePeakWidth(runs.streamer.on.zprof(:,iso), zvec);
        %runs.streamer.on.zwidth(iso) = zvec(ilo);
    end

    runs.streamer.hash = githash([mfilename('fullpath') '.m']);

    streamer = runs.streamer;
    save([runs.dir '/avgstreamer.mat'], 'streamer');

    disp('Done.');
    toc(ticstart);
end