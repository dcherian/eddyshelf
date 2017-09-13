% makes streamer mask for numerical parameterization
function [v, mask, rho, xvec, zvec] = ...
        makeStreamerSection(runs, isobath, maxloc, V0, L0, Lz0)

    debug = 0;

    phys = runs.params.phys;

    if ~exist('V0', 'var') || isempty(V0)
        if ~exist('maxloc', 'var') || isempty(maxloc)
            [~,maxloc] = runs.calc_maxflux(isobath);
        end

        %V0 = runs.eddy.rhovor.Vke(maxloc);
        vel = smooth(hypot(runs.eddy.fitx.V0, runs.eddy.fity.V0), 20) / 2.3 / sqrt(2);
        V0 = vel(maxloc);
    end

    if ~exist('L0', 'var') || isempty(L0)
        if ~exist('maxloc', 'var') || isempty(maxloc)
            [~,maxloc] = runs.calc_maxflux(isobath);
        end

        L0 = median(runs.eddy.rhovor.dia(1:maxloc))/2;
    end

    if ~exist('Lz0', 'var') || isempty(Lz0)
        if ~exist('maxloc', 'var') || isempty(maxloc)
            [~,maxloc] = runs.calc_maxflux(isobath);
        end

        Lz0 = runs.eddy.Lgauss(maxloc);
    end

    zvec = runs.csflux.vertbins(:, isobath);
    xvec = runs.rgrid.x_rho(1,2:end-1) - mean(runs.rgrid.x_rho(1,:));

    % normalized grid matrices to create mask
    [xmat, zmat] = ndgrid(xvec/L0, zvec/Lz0);

    R = runs.csflux.R;
    yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
                                      % y0oL =  R/L0 * (1 - yoR); % y0/L - used in derivation
                                      %xfrac = sqrt(1 - y0oL^2);
                                      %y0oL = (runs.eddy.my(maxloc) - runs.csflux.x(isobath))/L0;

    % eddy fields
    a = 2;
    v = -sqrt(2*exp(1)) * runs.sgntamp * V0 * xmat.^(a-1) .* exp(-xmat.^a) .* (1-erf(-zmat));

    xline = 0;

    eddymask = ((xmat.^a + zmat.^a) > 1.0^a);
    mask = (xmat < xline) & (eddymask);

    if debug
        figure;
        pcolorcen(xmat, zmat, v);
        center_colorbar;
        hold on
        try
            contour(xmat, zmat, kinkmask, 'k');
            contour(xmat, zmat, eddymask, 'r');
            linex(xline);
        catch ME
        end
        contour(xmat, zmat, mask, 'b');
    end
end
