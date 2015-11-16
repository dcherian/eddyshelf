% integrate flux down to 'depth' at 'isobath' for 'source'
function [fluxvec] = recalculateFlux(runs, depth, isobath, source, debug)

    if ~exist('depth'), error('Specify depth!'); end
    if ~exist('isobath'), error('Specify isobath!'); end
    if ~exist('source', 'var'), source = isobath; end
    if ~exist('debug', 'var'), debug = 0; end

    vertbins = runs.csflux.vertbins(:,isobath);
    slopezt = runs.csflux.off.slopezt(:,:, isobath, source);

    zind = find_approx(vertbins, -abs(depth), 1);

    if ~strcmpi(runs.name, 'ew-2041')
        nsmth = 7;
    else
        nsmth = 4;
    end
    fluxvec = smooth(trapz(vertbins(zind:end), slopezt(zind:end,:), 1)', nsmth);

    if debug
        figure;
        plot(fluxvec);
        hold on
        plot(runs.csflux.off.slope(:,isobath,source));
        legend(num2str(depth), 'Total');
    end
end