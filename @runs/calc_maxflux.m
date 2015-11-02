% return magnitude and time index of max flux
function [maxflux, maxloc] = calc_maxflux(runs, fluxin, debug)

    if ~exist('debug', 'var'), debug = 0; end

    if numel(fluxin) == 1
        fluxin = runs.csflux.off.slope(:,fluxin,fluxin);
    end

    iflux = cumtrapz(runs.csflux.time, fluxin);
    [start,stop] = runs.flux_tindices(fluxin);

    % This smoothing is used only to make peak detection work
    % better. i.e., this makes maxloc a better estimate neglecting
    % some smaller peaks that might occur earlier.
    % The maxflux value is from the actual *unsmoothed*  time
    % series at time index "maxloc"
    nsmooth = 10;

    % flux vector for applicable time
    % convert everything to double since I'm
    % dealing with large numbers here (time in
    % seconds) and integrated flux (m^3)
    fluxvec = double(smooth(fluxin(start:stop), nsmooth));
    tvec = double(runs.csflux.time(start:stop));

    mpp = 0.1; % min. peak prominence
    mpw = 6; % peaker wider than mpw points
    [~,locs] = findpeaks(fluxvec, ...
                         'MinPeakProminence', mpp*max(fluxvec(:)), ...
                         'MinPeakWidth', mpw);

    if debug %| isempty(locs)
        if isempty(locs), warning('No peaks found!'); end

        hfig = figure;
        findpeaks(fluxvec, ...
                  'MinPeakProminence', mpp*max(fluxvec(:)), ...
                  'MinPeakWidth', mpw ...
                  ); %'annotate', 'extents');
        %        plot(fluxin);
        title(runs.name);
    end

    if isempty(locs)
        [maxflux, maxloc] = max(fluxin);
        warning(['Couldn''t find maxflux for ' runs.name '. Using max!!!']);
        return;
    end

    % low Rh runs seem to lose annulus fluid earlier.
    % for really deep shelfbreak runs, the second peak
    % is artificially small because of the splitting.
    if ~isempty(findstr(runs.name, 'ew-0')) | ...
            ~isempty(findstr(runs.name, '1_wider')) | ...
            ~isempty(findstr(runs.name, '3_wider'))
        maxloc = locs(1);
    end
    if ~isempty(strfind(runs.name, 'ew-40')) ...
            | ~isempty(strfind(runs.name, 'ew-83'))
        maxloc = locs(1);
    end

    % For some reason there is a secondary peak.
    % This might be fixed with a peak width criterion
    % but that might screw up something else.
    if strcmpi(runs.name, 'ew-35')
        maxloc = locs(3);
    end
    if strcmpi(runs.name, 'ew-2365-75km')
        maxloc = locs(2);
    end
    % with bottom friction, choose first peak
    if ~isempty(strfind(runs.name, 'ew-5'))
        maxloc = locs(1);
    end

    if length(locs) == 1
        maxloc = locs(1);
    end

    % select the second peak by default
    % This should, in general, correspond to annulus fluid
    % being shed
    if ~exist('maxloc', 'var'), maxloc = locs(2); end

    % correct time shift
    maxloc = maxloc + start - 1;
    maxflux = fluxin(maxloc,1);

    if strcmpi(runs.name, 'ew-4341')
        [maxflux, maxloc] = max(fluxin);
    end

    %[maxflux, maxloc] = max(fluxin);

    if exist('hfig', 'var')
        linex(maxloc - start + 1);
    end

    runs.csflux.maxflux = maxflux;
    runs.csflux.maxloc = maxloc;
end