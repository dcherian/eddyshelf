% function [itfit, tfit, T, BadFitFlag] = FitCenterVelocity(runs, debug)
function [itfit, tfit, T, BadFitFlag] = FitCenterVelocity(runs, debug)

    if ~exist('debug', 'var'), debug = 0; end

    sgn = runs.sgntamp * sign(runs.params.phys.f0);
    % eddy center (SSH max) velocity
    cvy = sgn* runs.smoothCenterVelocity([], 'max');

    [~,imin] = min(cvy);
    imin = imin;
    tvec = runs.eddy.t*86400;

    [v0, T, t0, v1, exitflag] = gauss_fit(tvec(imin:end) - tvec(imin), ...
                                          (cvy(imin:end) - cvy(end))./cvy(imin), ...
                                          debug);
    tfit = T + tvec(imin) + t0;
    itfit = find_approx(tvec, tfit, 1);
    BadFitFlag = ~exitflag;

    if debug
        hax = gca;
        subplot(2,1,1,hax);
        linex(T); title(runs.name);
        subplot(2,1,2);
        plot(tvec, cvy);
        linex(tfit);
    end
end
