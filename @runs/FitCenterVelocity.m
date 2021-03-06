% function [itfit, tfit, T, BadFitFlag, FitTimeSeries] = FitCenterVelocity(runs, debug)
% returns vectors itfit = [estimate, lower bound, upper bound] for all variables
function [itfit, tfit, T, BadFitFlag, FitTimeSeries] = FitCenterVelocity(runs, debug)

    if ~exist('debug', 'var'), debug = 0; end

    sgn = runs.sgntamp * sign(runs.params.phys.f0);
    % eddy center (SSH max) velocity
    cvy = sgn* runs.smoothCenterVelocity(5, 'max');

    tvec = runs.eddy.t*86400;
    [tvec, uind] = unique(tvec);
    cvy = cvy(uind);

    [~,imin] = min(cvy);

    if strcmpi(runs.name, 'ew-8384')
        % weird stuff at the end
        [~,imin] = min(cvy(1:end-20));
    end
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

    %if strcmpi(runs.name, 'ew-6341-4')
        %cvy = smooth(cvy, 10);
        %end

    % if strcmpi(runs.name, 'ew-64461-6') ...
    %         | strcmpi(runs.name, 'ew-6341-2') ...
    %         | strcmpi(runs.name, 'ew-6341-4')
    %     imin = imin - 0;
    % end

    % if strcmpi(runs.name, 'ew-6441')
    %     imin = 51;
    % end

    % I want to use a few values before the peak as a constraint.
    % Sometimes this won't work if you don't average the velocity enough.
    dimin = 10; ceil(20/diff(tvec(1:2))*86400);
    if imin-dimin < 1, imin = dimin + 1; end

    %if strcmpi(runs.name, 'ew-56341-2'), cvy0 = 0; end
    [v0, T, t0, v1, exitflag, conf] = gauss_fit(tvec(imin-dimin:end) - tvec(imin), ...
                                                (cvy(imin-dimin:end) - cvy0)./cvy(imin), ...
                                                debug);
    Tconf = conf(:,2);
    t0conf = conf(:,3);

    % ASSUMES THAT ALL VELOCITY POINTS ARE INDEPENDENT ESTIMATES
    T = [T Tconf(1) Tconf(2)];
    tfit = T(1) + tvec(imin) + t0(1) + [0 -1 1]* sqrt( (t0conf(2)-t0(1))^2 + (Tconf(2)-T(1))^2);
    itfit = vecfind(runs.eddy.t*86400, tfit);
    BadFitFlag = ~exitflag;

    FitTimeSeries(:,2) = (v0 * exp(-((tvec(imin-dimin:end)-t0-tvec(imin))./T(1)).^2)  ...
                          + v1)*cvy(imin) + cvy0;
    FitTimeSeries(:,1) = tvec(imin-dimin:end);

    if any(isnan(itfit)), BadFitFlag = 1; end

    if debug
        hax = gca;
        subplot(2,1,1,hax);
        linex(T); title(runs.name);
        subplot(2,1,2);
        plot(tvec, cvy);
        linex([tfit t0]); liney(0);
        plot(FitTimeSeries(:,1), FitTimeSeries(:,2));

        fprintf('T = '); disp(T/86400);
        fprintf('t0 = '); disp([t0 t0conf']/86400);
    end
end
