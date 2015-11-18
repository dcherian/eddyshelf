% Predict depth of peak and width of peak in streamer structure
% flag == use : apply parameterization
% flag == detect : constant = 1;
% By defualt, flag = use.
function [width, zpeak] = predict_zpeak(runs, iso, flag, maxloc)

    if ~exist('flag', 'var'), flag = 'use'; end
    if ~exist('maxloc', 'var'), maxloc = []; end

    eddy = runs.params.eddy;
    phys = runs.params.phys;

    if ~isfield('eddy', 'rhovor'), width = NaN; zpeak = NaN; end
    if isempty(maxloc)
        [~,maxloc] = runs.calc_maxflux(iso);
    end
    Ly = runs.eddy.rhovor.dia(1)/2; % using maxloc does not change much
    %Lz = smooth(Lz, 30);
    %Ro = smooth(Ro, 1);
    %Ly = runs.eddy.rhovor.dia(1)/2;
    R = runs.csflux.R;
    yoR = runs.csflux.ndloc(iso);
    y0oL = R/Ly * (1 - yoR); % y0/L - used in derivation
    xfrac = 0.7;

    % calculate density profiles
    [rhoshelf, zshelf] = runs.getDensityProfile(runs.bathy.isb);
    [rhoslope, zslope] = runs.getDensityProfile(runs.csflux.ix(iso));

    rhobot = rhoshelf(1);
    eddy.temp = eddy.tamp .* exp(-xfrac^2 - y0oL^2) .* ...
        exp(-(zslope/eddy.depth).^2);
    eddy.rho = phys.R0*(- phys.TCOEF * eddy.temp);
    zind = find_approx(eddy.rho + rhoslope, rhobot, 1);

    if strcmpi(flag, 'use')
        zpw = load('./params/param_zpeakwidth.mat');
        index = zpw.isobath == iso;
        slope = zpw.slope(index);
        intercept = zpw.intercept(index);
    end

    if ~strcmpi(flag, 'use') || isempty(slope)
        slope = 1;
        intercept = 0;
    end

    width = slope * zslope(zind) + intercept;
    zpeak = width/2;
end