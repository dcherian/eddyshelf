function [avgflux, err] = calc_avgflux(runs, fluxvec, debug)

    if ~exist('debug', 'var'), debug = 0; end
    use_wunsch = 0;

    [start,stop] = runs.flux_tindices(fluxvec);
    nsmooth = 1;

    % flux vector for applicable time
    % convert everything to double since I'm
    % dealing with large numbers here (time in
    % seconds) and integrated flux (m^3)
    fluxvec = double(smooth(fluxvec(start:stop), nsmooth));
    tvec = double(runs.csflux.time(start:stop));
    tvec = (tvec - tvec(1))';

    avgflux = mean(fluxvec);
    dof = calcdof(fluxvec);
    err = abs(conft(0.05, dof-1) * std(fluxvec) / sqrt(dof));

    % check error bounds
    if debug
        tvec = tvec / 86400;
        figure;
        subplot(2,1,1);
        plot(tvec, fluxvec);
        liney([avgflux-err avgflux avgflux+err]);

        subplot(2,1,2);
        plot(tvec,ifluxvec);
        hold on
        plot([0 tvec(end)], [0 (avgflux - err)*86400*tvec(end)]);
        plot([0 tvec(end)], [0 avgflux*86400*tvec(end)]);
        plot([0 tvec(end)], [0 (avgflux + err)*86400*tvec(end)]);
    end
end