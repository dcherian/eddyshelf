function [avgflux, err] = calc_avgflux(runs, fluxvec, debug, ...
                                       whenstart, whenstop)

    if ~exist('debug', 'var'), debug = 0; end
    if ~exist('whenstart', 'var'), whenstart = []; end
    if ~exist('whenstop', 'var'), whenstop = []; end
    use_wunsch = 0;

    [start,stop] = runs.flux_tindices(fluxvec, whenstart, whenstop);
    nsmooth = 1;

    % flux vector for applicable time
    % convert everything to double since I'm
    % dealing with large numbers here (time in
    % seconds) and integrated flux (m^3)
    fluxvec = double(smooth(fluxvec(start:stop), nsmooth));
    tvec = double(runs.csflux.time(start:stop));
    tvec = (tvec - tvec(1))';

    [tvec, uind] = unique(tvec);
    fluxvec = fluxvec(uind);

    if mean(diff(tvec)) ~= diff(tvec(1:2))
        % interpolate to constant grid
        warning('flux time series Δt is not constant');
        tvecnew = tvec(1):diff(tvec(1:2)):tvec(end);
        fluxvecnew = interp1(tvec, fluxvec, tvecnew);

        tvec = tvecnew;
        fluxvec = fluxvecnew;
    end

    avgflux = mean(fluxvec);

    [dof,IT] = calcdof(fluxvec);
    err = abs(conft(0.05, dof-1) * std(fluxvec) / sqrt(dof));

    % check error bounds
    if debug
        tvec = tvec / 86400;
        figure;
        subplot(2,1,1);
        plot(tvec, fluxvec);
        liney([avgflux-err avgflux avgflux+err]);
        title(runs.name);

        subplot(2,1,2);
        plot(tvec,cumtrapz(tvec*86400, fluxvec));
        hold on
        plot([0 tvec(end)], [0 (avgflux - err)*86400*tvec(end)]);
        plot([0 tvec(end)], [0 avgflux*86400*tvec(end)]);
        plot([0 tvec(end)], [0 (avgflux + err)*86400*tvec(end)]);
    end
end