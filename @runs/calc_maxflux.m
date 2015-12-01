% return magnitude and time index of max flux
function [maxflux, maxloc, err] = calc_maxflux(runs, fluxin, isobath, debug)

    if ~exist('isobath', 'var'), isobath = 0; end
    if ~exist('debug', 'var'), debug = 0; end

    if numel(fluxin) == 1
        fluxin = runs.csflux.off.slope(:,fluxin,fluxin);
    end

    iflux = cumtrapz(runs.csflux.time, fluxin);
    [start,stop] = runs.flux_tindices(fluxin);
    err = 0;

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
    % if ~isempty(findstr(runs.name, 'ew-0')) | ...
    %         ~isempty(findstr(runs.name, 'ew-1')) | ...
    %         ~isempty(findstr(runs.name, '1_wider')) | ...
    %         ~isempty(findstr(runs.name, '3_wider'))
    %     maxloc = locs(1);
    % end
    %if ~isempty(strfind(runs.name, 'ew-40')) ...
    %        | ~isempty(strfind(runs.name, 'ew-83'))
    %    maxloc = locs(1);
    %end

    % For some reason there is a secondary peak.
    % This might be fixed with a peak width criterion
    % but that might screw up something else.
    % if strcmpi(runs.name, 'ew-35') | strcmpi(runs.name, 'ew-36')
    %     maxloc = locs(1);
    % end
    % if strcmpi(runs.name, 'ew-37')
    %     maxloc = locs(3);
    % end

    % if strcmpi(runs.name, 'ew-2365-75km')
    %     maxloc = locs(1);
    % end

    % % with bottom friction, choose first peak
    % % sloping shelf too
    % if ~isempty(strfind(runs.name, 'ew-8'))
    %     maxloc = locs(1);
    % end

    % if ~isempty(strfind(runs.name, 'ew-5'))
    %     [maxflux, maxloc] = max(fluxin);
    % end

    %if length(locs) == 1
    %    maxloc = locs(1);
    %end

    % by default, plain maximum
    [~, maxloc] = max(fluxvec);

    % choose second peak everywhere for ew-4343 and ew-2041
    if strcmpi(runs.name, 'ew-2041') | strcmpi(runs.name, 'ew-4343')
        if length(locs) > 1
            maxloc = locs(2);
        else
            maxloc = locs(1);
        end
    end

    %if strcmpi(runs.name, 'ew-2365-75km')
    %    maxloc = locs(1);
    %end

    if strcmpi(runs.name, 'ew-04')
        maxloc = locs(1);
    end
    if strcmpi(runs.name, 'ew-06')
        maxloc = locs(1);
    end
    if strcmpi(runs.name, 'ew-2340') | strcmpi(runs.name, 'ew-2345')
        maxloc = locs(3);
    end

    maxloc = maxloc + start - 1;
    maxflux = fluxin(maxloc,1);

    % maxflux = median(fluxvec(locs));
    if length(locs) > 2
        err = std(fluxvec(locs))/sqrt(length(locs));
    else
        err = 0;
    end

    if exist('hfig', 'var')
        linex(maxloc - start + 1);
    end

    runs.csflux.maxflux = maxflux;
    runs.csflux.maxloc = maxloc;
end