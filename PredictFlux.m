function [] = PredictFlux()
% Predict flux based on provided parameters

    % physical properties
    phys.f0 = 2*(2*pi/86400)*sind(38); % Coriolis parameter (1/s)

    % eddy properties
    eddy.V0 = 2; % maximum velocity in eddy (m/s)
    eddy.L0 = 70e3; % horizontal scale (m)
    eddy.Lz0 = 700; % vertical scale (m)
    eddy.sign = -1; % Anticyclone = -1, cyclone = +1

    % flux properties
    flux.IntegrationDepth = 100; % depth to which to integrate (m)
    flux.IsobathLocation = 0e3; % from shelfbreak (m)
    flux.IsobathDepth = 100; % (m)

    plot_mask = 0; % plot mask?

    [v,mask,xvec, zvec] = makeEddyStreamerSection(phys, eddy, flux, plot_mask);

    zind = find_approx(zvec, -1 * abs(flux.IntegrationDepth));
    Flux = trapz(xvec, trapz(zvec(zind:end), v(:,zind:end) .* mask(:,zind:end), 2), 1);

    yoR = [0 0.17 0.33 0.5 0.67 0.83 1 1.17];
    %RegressionSlopes = [0.02 0.12 0.16 0.2 0.24 0.26 0.28 0.3];
    %Uncertainty = [0.0 0.01 0.02 0.03 0.03 0.03 0.03 0.03];
    RegressionSlopes = [0.1394 0.1161 0.1628 0.2026 0.2409 0.2600 0.2830 0.3054];
    Uncertainty = [0.0225 0.0213 0.0298 0.0354 0.0369 0.0338 0.0338 0.0359];

    Slope = interp1(yoR, RegressionSlopes, flux.IsobathLocation./eddy.L0);
    Error = interp1(yoR, Uncertainty, flux.IsobathLocation./eddy.L0);

    fprintf('Flux is %0.2f ± %0.2f Sv\n', Flux * Slope/1e6, Flux * Error/1e6);
end

function [v, mask, xvec, zvec] = makeEddyStreamerSection(phys, eddy, flux, plot_mask)

    zvec = linspace(-1*abs(flux.IsobathDepth), 0, 80);
    xvec = linspace(-5*eddy.L0,0,100);

    % normalized grid matrices to create mask
    [xmat, zmat] = ndgrid(xvec/eddy.L0, zvec/eddy.Lz0);

    % eddy velocity field at latitude (cross-isobath location) of eddy center
    v = 2.3 * eddy.sign * eddy.V0 * xmat .* exp(-xmat.^2) .* (1-erf(-zmat));

    eddymask = ((xmat.^2 + zmat.^2) > 1.0^2);

    mask = (xmat < 0) & (eddymask);

    if plot_mask
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

function g = find_approx(m,v,n)
%     g = find_approx(m,v,[n])
% find index (g) of matrix (m) that is most nearly equal to value (v); this
% is similar to g=find(m==v), except that the nearest approximate equality
% is found if no exact equality exists.
% the third argument (default n=1) tells how many values to find; n=3 means
% the nearest 3 indices in order of descending nearness.

    g=find(m==v);
    if isempty(g)
        [nul g]=min(abs(m-v));
        if isnan(nul), g=nan; end
    end

    if nargin>2
        g=zeros(n,1)*nan;
        for nn=1:n
            [nul g(nn)]=min(abs(m-v));
            m(g(nn))=nan;
        end
    end
end