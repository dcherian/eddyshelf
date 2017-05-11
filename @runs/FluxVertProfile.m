function [profile] = FluxVertProfile(runs, isobath)

    if ~exist('isobath', 'var'), isobath = 1; end

    [start, stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1));
    offflux = runs.csflux.off.slopezt(:,start:stop,1,1);
    zivec = runs.csflux.vertbins(:,isobath);
    profile = trapz(runs.csflux.time(start:stop)*86400, ...
                    offflux, 2);
    profile = profile./max(abs(profile));
