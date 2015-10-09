% shelf water density profile at shelfbreak
function [rho, zr] = getShelfDensityProfile(runs, index)

    if runs.bathy.axis == 'y'
        zrmat = runs.rgrid.z_r(:,:,1);
    else
        zrmat = squeeze(runs.rgrid.z_r(:,1,:));
    end

    zr = zrmat(:,index);
    
    flags = runs.params.flags;
    phys = runs.params.phys;
    
    % Create background state (assumes uniform horizontal grid)
    % assign initial stratification
    temp = phys.T0.*ones(size(zr));

    % this accounts for when I'm extracting parameters from runs that
    % didn't have flags.conststrat
    if ~isfield(flags, 'conststrat')
        flags.conststrat = 1;
    end

    % N2 here is phys.N2 = strat.N2
    if flags.conststrat
        % constant stratification
        Tz = phys.N2/phys.g/phys.TCOEF * ones(size(zr) - [1 0]);
    else
        error('doesn''t work yet');
        if strat.z0 > 0 , strat.z0 = strat.z0 * -1; end
        % non-constant stratification.
        zmat = zwmat(:,:,2:end-1);
        N2mat = strat.N2max .* ...
                ((exp(-(zmat - strat.z0)./strat.Lp) .* (zmat >  strat.z0)) + ...
                 (exp( (zmat - strat.z0)./strat.Lm) .* (zmat <= strat.z0)));
        % clamp min N² to 1e-6
        N2mat(N2mat < 1e-6) = 1e-6;
        Tz = N2mat./phys.g./phys.TCOEF;

        % depth-averaged N² is calculated later

        % estimate slope Burger number using max. N² at the _bottom_
        N2bot = N2mat(:,:,1);
        bathy.S_sl = sqrt(max(N2bot(:))) .* bathy.sl_slope ./ phys.f0;
    end

    for k=length(zr)-1:-1:1
        temp(k) = temp(k+1) - Tz(k).*(zr(k+1)-zr(k));
    end

    rho = phys.R0*(1 - phys.TCOEF * (temp-phys.T0));
end