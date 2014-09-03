
% calculate integrated transport and avg flux given a flux
% vector and time vector - used by fluxes().
function [itrans, avgflux] = integrate_flux(runs, time, flux)

    if any(isnan(flux))
        error('FLUX IS NAN!');
    end

    ind = runs.eddy.tscaleind;
    itrans = cumtrapz(time, flux);

    % start where flux is 5% of max. and stop where
    % integrated transport is 95% of max.
    start = ind + find(abs(flux(ind:end)) > 0.05 * max(abs(flux(ind:end))), ...
                       1, 'first');
    stop = find_approx(abs(itrans), 0.95 * max(abs(itrans)), 1);

    % calc average flux
    avgflux = (itrans(stop)-itrans(start)) ./ ...
              (time(stop)-time(start));
end