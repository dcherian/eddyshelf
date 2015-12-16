% estimate start and end of flux time series that I will use to
% calculate diagnostics: max and avg.

function [start, stop] = flux_tindices(runs, flux)

    itrans = cumtrapz(flux);

    ind = 1;

    % start where integrated transport is 5% of max. and stop where
    % integrated transport is 95% of max.
    %start = ind + find(abs(flux(ind:end)) > 0.25 * max(abs(flux(ind:end))), ...
    %                   1, 'first');
    % Sep 09, 2015: change from 0.02 because of 3341-2-big
    start = find_approx(abs(itrans), 0.05 * max(abs(itrans)), 1);
    stop = find_approx(abs(itrans), 0.90 * max(abs(itrans)), 1);
end