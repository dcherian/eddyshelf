% makes streamer mask for numerical parameterization
function [v, mask] = makeStreamerSection(runs, isobath, maxloc, V0, L0, Lz0)

    debug = 0;

    if ~exist('maxloc', 'var') || isempty(maxloc)
        [~,maxloc] = runs.calc_maxflux(isobath);
    end

    if ~exist('V0', 'var') || isempty(V0)
        V0 = runs.eddy.rhovor.Vke(maxloc);
    end

    if ~exist('L0', 'var') || isempty(L0)
        L0 = median(runs.eddy.rhovor.dia(1:maxloc))/2;
    end

    if ~exist('Lz0', 'var') || isempty(Lz0)
        Lz0 = runs.eddy.Lgauss(maxloc);
    end

    zvec = runs.csflux.vertbins(:, isobath);
    xvec = (runs.rgrid.x_rho(1,2:end-1) - runs.eddy.mx(maxloc));

    % normalized grid matrices to create mask
    [xmat, zmat] = ndgrid(xvec/L0, zvec/Lz0);

    R = runs.csflux.R;
    yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
    y0oL =  R/L0 * (1 - yoR); % y0/L - used in derivation
    xfrac = sqrt(1 - y0oL^2);

    %y0oL = (runs.eddy.my(maxloc) - runs.csflux.x(isobath))/L0;
    %if isobath < 3, y0oL = 0; end
    v = -2.3 * V0 * xmat .* exp(-xmat.^2 - y0oL.^2) .* (1-erf(-zmat));

    [width, zpeak] = runs.predict_zpeak(isobath, 'use');
    width = abs(width/Lz0); zpeak = abs(zpeak/Lz0);

    %kzrad = width/2; % kink radius - z
    %kxrad = kzrad; % kink radius - x
    %x0 = -xfrac-kxrad; -xfrac-kxrad;
    %z0 = -1 * width/3;
    xline = 0; -xfrac;

    a = 2;
    if ~isreal(xfrac)
        % complex xfrac -- cannot be trusted
        % make the kink (semi-circle) intersect the eddy contour
        xfrac = sqrt(1 - (width)^2);
        x0 = -xfrac;
        xline = 0;
    end

    %if abs(runs.csflux.x(isobath) - runs.bathy.xsb) < 2000
        % if close to shelfbreak use barotropic mask
        %    mask = xmat < 0;
        %v = repmat(v(:,1), [1 size(v,2)]);
        %else
        eddymask = ((xmat.^a + zmat.^a) > 1.0^a) .* (zmat < -width);
        %kinkmask = (((xmat-x0)/kxrad).^2 + ((zmat-z0)/kzrad).^2) <= 1;
        kinkmask = ((xmat.^a + zmat.^a) > 0.7^a) .* (zmat >= -width);
        mask = (xmat < xline) & (eddymask | kinkmask);
        %end

    if debug
        figure;
        hold on
        contour(xmat, zmat, kinkmask, 'k');
        contour(xmat, zmat, eddymask, 'r');
        linex(xline);
        contour(xmat, zmat, mask, 'b');
        keyboard;
    end
end
