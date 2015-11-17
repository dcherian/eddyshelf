function [avgflux, err] = calc_avgflux(runs, fluxvec, debug)

    if ~exist('debug', 'var'), debug = 0; end
    use_wunsch = 0;

    [start,stop] = runs.flux_tindices(fluxvec);
    nsmooth = 1;

    ifluxvec = cumtrapz(runs.csflux.time, fluxvec);

    % flux vector for applicable time
    % convert everything to double since I'm
    % dealing with large numbers here (time in
    % seconds) and integrated flux (m^3)
    fluxvec = double(smooth(fluxvec(start:stop), nsmooth));
    ifluxvec = double(smooth(ifluxvec(start:stop), nsmooth));
    tvec = double(runs.csflux.time(start:stop));

    % change origin
    ifluxvec = (ifluxvec - ifluxvec(1));
    tvec = (tvec - tvec(1))';

    % E = [ones(size(tvec))' tvec'];

    % if use_wunsch
    %     %%%%%%%%%%% See Wunsch(1996) pg. 116
    %     x = E\ifluxvec;
    %     intercept = x(1);
    %     avgflux = x(2);
    % else
    %     %%%%%%%%%%% use MATLAB regress
    %     [b, bint, r, rint, stats] = regress(ifluxvec, E);
    %     intercept = b(1);
    %     avgflux = b(2);
    %     err = abs(bint(2) - b(2));
    % end

    % true = ifluxvec;
    % est = intercept + avgflux .* (tvec-tvec(1))';
    % res = true-est;

    % if use_wunsch
    %     % (E' * E) ^-1
    %     %ETEI = inv(E'*E);
    %     % from http://blogs.mathworks.com/loren/2007/05/16/purpose-of-inv/
    %     [Q,R] = qr(E,0);
    %     S = inv(R);
    %     ETEI = S*S';
    %     % assuming noise vector (res) is white
    %     P = ETEI * E' * var(res) * E * ETEI;
    %     err = sqrt(diag(P));
    %     err = err(2); % standard error
    % end

    % runs.csflux.avgflux = avgflux;
    % runs.csflux.err = err;

    % % plot fit
    % if debug
    %     figure; hold all;
    %     plot(tvec/86400, true, '*');
    %     plot(tvec/86400, est); plot(tvec/86400, res); liney(0);
    %     title(runs.name);
    % end

    %[c,lags] = xcorr(fluxvec - mean(fluxvec(:)), 'coef');
    %plot(lags, c); linex(0); liney(0);

    %%%%%%%%%%% mean of instantaneous flux
    % find number of peaks
    mpd = 6;
    % crests
    [~,pl] = findpeaks(fluxvec, 'MinPeakDistance', mpd);
    % troughs
    [~,nl] = findpeaks(-1*fluxvec, 'MinPeakDistance', mpd); ...

    % make sure peak to trough distance is not
    % smaller than mpd
    indices = sort([pl; nl]);
    mask = [0; diff(indices) < mpd];
    filtered = indices(~isnan(fillnan(indices .* ~mask,0)));
    dof = length(filtered) + 1; % (crude) degrees of freedom;

    % check dof calculation
    %figure; plot(fluxvec); linex(filtered); title(num2str(dof));pause;

    % error bounds
    if dof == 2
        warning('Only 2 degrees of freedom. skipping');
    end

    err = abs(conft(0.05, dof-1) * std(fluxvec) / sqrt(dof));
    avgflux = mean(fluxvec);

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
