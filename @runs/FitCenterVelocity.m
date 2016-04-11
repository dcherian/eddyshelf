% function [itfit, tfit, T, BadFitFlag] = FitCenterVelocity(runs, debug)
% returns vectors itfit = [estimate, lower bound, upper bound] for all variables
function [itfit, tfit, T, BadFitFlag] = FitCenterVelocity(runs, debug)

    if ~exist('debug', 'var'), debug = 0; end

    sgn = runs.sgntamp * sign(runs.params.phys.f0);
    % eddy center (SSH max) velocity
    cvy = sgn* runs.smoothCenterVelocity(5, 'max');

    [~,imin] = min(cvy);
    imin = imin;
    tvec = runs.eddy.t*86400;

    [tvec, uind] = unique(tvec);
    cvy = cvy(uind);

    % OK default choice.
    % For some special cases (below), I set this to zero to get better results.
    % I think if the constant is ≳ 20% of max amplitude,
    % then the fit is hard to get right
    cvy0 = cvy(end);
    if abs(cvy0./cvy(imin)) > 0.1, cvy0 = 0; end

    if strcmpi(runs.name, 'ew-64461-4') | strcmpi(runs.name, 'ew-64461-7');
        % this one appears to split, so cut off the end
        cvy = cvy(1:220);
        tvec = tvec(1:220);
        cvy0 = 0;
    end

    if strcmpi(runs.name, 'ew-64461-5')
        % gets to sponge?
        cvy = cvy(1:250);
        tvec = tvec(1:250);
        cvy0 = 0;
    end

    if strcmpi(runs.name, 'ew-6341-2')
        % not long enouhg. large oscillations at the end
        % fit to whole record to constrain
        imin = 11;
    end
    %if strcmpi(runs.name, 'ew-6341-4')
        %cvy = smooth(cvy, 10);
        %end

    %if strcmpi(runs.name, 'ew-6362-1')
        % needs extra smoothing
        % cvy = smooth(cvy, 5);
    %end

    %if strcmpi(runs.name, 'ew-56341-2'), cvy0 = 0; end

    [v0, T, t0, v1, exitflag, conf] = gauss_fit(tvec(imin-10:end) - tvec(imin), ...
                                                (cvy(imin-10:end) - cvy0)./cvy(imin), ...
                                                debug);
    Tconf = conf(:,2);

    % ASSUMES THAT ALL VELOCITY POINTS ARE INDEPENDENT ESTIMATES
    T = [T Tconf(1) Tconf(2)];
    tfit = T + tvec(imin) + t0;
    itfit = vecfind(runs.eddy.t*86400, tfit);
    BadFitFlag = ~exitflag;

    if any(isnan(itfit)), BadFitFlag = 1; end

    if debug
        hax = gca;
        subplot(2,1,1,hax);
        linex(T); title(runs.name);
        subplot(2,1,2);
        plot(tvec, cvy);
        linex(tfit);
    end
end
