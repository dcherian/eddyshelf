% calculate integrated transport and avg flux given a flux
% vector and time vector - used by fluxes().
function [itrans, avgflux] = integrate_flux(runs, time, flux)

    if any(isnan(flux))
        warning('FLUX IS NAN!');
        itrans = zeros(size(flux));
        flux = zeros(size(flux));
        return;
    end

    itrans = cumtrapz(time, flux);

    [start,stop] = runs.flux_tindices(flux);

    if isempty(start)
        start = 1;
        warning('Fluxes are zero?');
        avgflux = 0;
        size(itrans)
        return;
    end

    % calc average flux
    avgflux = (itrans(stop)-itrans(start)) ./ ...
              (time(stop)-time(start));
end