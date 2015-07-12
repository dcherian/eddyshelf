% return magnitude and time index of max flux
function [maxflux, maxloc] = calc_maxflux(runs)

    debug = 0;
    
    ind = runs.csflux.tscaleind;
    nsmooth = 10;

    % flux vector for applicable time
    % convert everything to double since I'm
    % dealing with large numbers here (time in
    % seconds) and integrated flux (m^3)
    fluxvec = double(smooth( ...
        runs.csflux.west.shelf(ind:end,1), nsmooth));
    ifluxvec = double(smooth( ...
        runs.csflux.west.itrans.shelf(ind:end,1), nsmooth));
    tvec = double(runs.csflux.time(ind:end));

    [~,locs] = findpeaks(fluxvec);
    
    maxloc = locs(1);
    maxflux = fluxvec(maxloc);

    % correct time shift
    maxloc = maxloc + ind - 1;

    if debug
        figure;
        %findpeaks(fluxvec, 'annotate', 'extents');
        plot(runs.csflux.west.shelf(:,1));
        linex((maxloc));
        title(runs.name);
    end

    runs.csflux.maxflux = maxflux;
    runs.csflux.maxloc = maxloc;
end