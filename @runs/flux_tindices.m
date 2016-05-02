% estimate start and end of flux time series that I will use to
% calculate diagnostics: max and avg.
% [start, stop] = flux_tindices(runs, flux, whenstart, whenstop)
function [start, stop] = flux_tindices(runs, flux, whenstart, whenstop)

    if ~exist('whenstart', 'var'), whenstart = 0.05; end
    if ~exist('whenstop', 'var'), whenstop = 0.90; end

    itrans = cumtrapz(flux);

    ind = 1;

    % start where integrated transport is 5% of max. and stop where
    % integrated transport is 95% of max.
    %start = ind + find(abs(flux(ind:end)) > 0.25 * max(abs(flux(ind:end))), ...
    %                   1, 'first');
    % Sep 09, 2015: change from 0.02 because of 3341-2-big
    start = find_approx(abs(itrans), whenstart * max(abs(itrans)), 1);
    stop = find_approx(abs(itrans), whenstop * max(abs(itrans)), 1);
end