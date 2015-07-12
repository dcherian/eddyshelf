function [avgflux, err] = calc_avgflux(runs)

    ind = runs.csflux.tscaleind;
    nsmooth = 1;

    % flux vector for applicable time
    % convert everything to double since I'm
    % dealing with large numbers here (time in
    % seconds) and integrated flux (m^3)
    fluxvec = double(smooth( ...
        runs.csflux.west.shelf(ind:end,1), nsmooth));
    ifluxvec = double(smooth( ...
        runs.csflux.west.itrans.shelf(ind:end,1), nsmooth));
    tvec = double(runs.csflux.time(ind:end));

    % change origin
    ifluxvec = (ifluxvec - ifluxvec(1));
    tvec = (tvec - tvec(1));

    E = [ones(size(tvec))' tvec'];

    %%%%%%%%%%% See Wunsch(1996) pg. 116
    % P matrix
    x = E\ifluxvec;
    intercept = x(1);
    avgflux = x(2);
    true = ifluxvec; est = intercept + avgflux .* ...
           (tvec-tvec(1))';
    res = true-est;
    % (E' * E) ^-1
    %ETEI = inv(E'*E);
    % from http://blogs.mathworks.com/loren/2007/05/16/purpose-of-inv/
    [Q,R] = qr(E,0);
    S = inv(R);
    ETEI = S*S';
    % assuming noise vector (res) is white
    P = ETEI * E' * var(res) * E * ETEI;
    err = sqrt(diag(P));
    err = err(2); % standard error

    %%%%%%%%%%% use MATLAB regress
    [b, bint, r, rint, stats] = ...
        regress(ifluxvec, E);

    avgflux = b(2);
    err = abs(bint(2) - b(2));

    runs.csflux.avgflux = avgflux;
    runs.csflux.err = err;

    % plot fit
    %figure; hold all;
    %plot(tvec/86400, true, '*');
    %plot(tvec/86400, est); plot(tvec/86400, res); liney(0);

    %[c,lags] = xcorr(fluxvec - mean(fluxvec(:)), 'coef');
    %plot(lags, c); linex(0); liney(0);

    %%%%%%%%%%% mean of instantaneous flux
    % find number of peaks
    %mpd = 6;
    % crests
    %[~,pl] = findpeaks(fluxvec, 'MinPeakDistance', mpd);
    % troughs
    %[~,nl] = findpeaks(-1*fluxvec, 'MinPeakDistance', mpd); ...

    % make sure peak to trough distance is not
    % smaller than mpd
    %indices = sort([pl; nl]);
    %mask = [0; diff(indices) < mpd];
    %filtered = indices(~isnan(fillnan(indices .* ~mask,0)));
    %dof = length(filtered) + 1; % (crude) degrees of freedom;

    % check dof calculation
    %figure; plot(fluxvec); linex(filtered); title(num2str(dof));pause;

    %flx = mean(max(fluxvec/1000));
    %flx = runs.csflux.west.avgflux.shelf(1)/1000;
    % standard deviation
    %sdev = sqrt(1./(length(fluxvec)-1) .* sum((fluxvec - flx*1000).^2))/1000;
    % error bounds
    %errmean = abs(conft(0.05, dof-1) * sdev / sqrt(dof));
% $$$ % $$$
% $$$                     % check error bounds with itrans
% $$$                     hfig2 = figure;
% $$$                     set(gcf, 'renderer', 'opengl');
% $$$                     subplot(2,1,1);
% $$$                     plot(runs.csflux.time/runs.tscale, ...
% $$$                          runs.csflux.west.shelf(:,1)/1000);
% $$$                     liney([flx-err flx flx+err]);
% $$$                     subplot(2,1,2);
% $$$                     plot(runs.csflux.time/runs.tscale, ...
% $$$                          runs.csflux.west.itrans.shelf(:,1));
% $$$                     Ln = createLine(1, ...
% $$$                                    runs.csflux.west.itrans.shelf(runs.tscaleind,1), ...
% $$$                                    1, (flx-err)*1000*runs.tscale);
% $$$                     L = createLine(1, ...
% $$$                                    runs.csflux.west.itrans.shelf(runs.tscaleind,1), ...
% $$$                                    1, flx*1000*runs.tscale);
% $$$                     Lp = createLine(1, ...
% $$$                                    runs.csflux.west.itrans.shelf(runs.tscaleind,1), ...
% $$$                                    1, (flx+err)*1000*runs.tscale);
% $$$                     hold on; drawLine(L);drawLine(Ln,'Color','g'); ...
% $$$                         drawLine(Lp,'Color','r');
% $$$                     limy = ylim; ylim([0 limy(2)]);
% $$$                     %pause;
% $$$                     try
% $$$                         close(hfig2);
% $$$                     catch ME; end
% $$$ % $$$
end
