% function [itfit, tfit, T, BadFitFlag] = FitCenterVelocity(runs, debug)
% returns vectors itfit = [estimate, lower bound, upper bound] for all variables
function [itfit, tfit, T, BadFitFlag] = FitCenterVelocity(runs, debug)

    if ~exist('debug', 'var'), debug = 0; end

    sgn = runs.sgntamp * sign(runs.params.phys.f0);
    % eddy center (SSH max) velocity
    cvy = sgn* runs.smoothCenterVelocity([], 'max');

    [~,imin] = min(cvy);
    imin = imin;
    tvec = runs.eddy.t*86400;

    [v0, T, t0, v1, exitflag, conf] = gauss_fit(tvec(imin:end) - tvec(imin), ...
                                                (cvy(imin:end) - cvy(end))./cvy(imin), ...
                                                debug);
    Tconf = conf(:,2);

    % ASSUMES THAT ALL VELOCITY POINTS ARE INDEPENDENT ESTIMATES
    T = [T Tconf(1) Tconf(2)];
    tfit = T + tvec(imin) + t0;
    itfit = vecfind(tvec, tfit);
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
