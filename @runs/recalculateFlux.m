% integrate flux down to 'depth' at 'isobath' for 'source'
function [fluxvec] = recalculateFlux(runs, depth, isobath, source, debug)

    if ~exist('depth'), error('Specify depth!'); end
    if ~exist('isobath'), error('Specify isobath!'); end
    if ~exist('source', 'var'), source = isobath; end
    if ~exist('debug', 'var'), debug = 0; end

    vertbins = runs.csflux.vertbins(:,isobath);
    slopezt = runs.csflux.west.slopezt(:,:, isobath, source);

    zind = find_approx(vertbins, -abs(depth), 1);

    fluxvec = trapz(vertbins(zind:end), slopezt(zind:end,:), 1);

    if debug
        figure;
        plot(fluxvec);
        hold on
        plot(runs.csflux.west.slope(:,isobath,source));
        legend(num2str(depth), 'Total');
    end
end