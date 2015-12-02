% locate point of "resistance"
% defined as start of half the maximum southward speed.
% returns center locations and time index
% nsmooth is number of turnover periods to smooth over
% factor * v_min is what I'm looking for.
function [xx,yy,tind,cvy] = locate_resistance(runs, nsmooth, factor)

    debug_plot = 0;

    if ~exist('factor', 'var'), factor = 1/4; end
    if ~exist('nsmooth', 'var'), nsmooth = []; end

    % support cyclones too
    sgn = sign(runs.params.eddy.tamp) * ...
          sign(runs.params.phys.f0);
    vel = sgn * runs.smoothCenterVelocity(nsmooth);

    vel = vel(1:runs.eddy.tend-1);
    tvec = runs.eddy.t(1:runs.eddy.tend-1);
    ndtime = tvec*86400./runs.eddy.turnover;

    % locate minimum, i.e., max. southward/westward speed
    % but limit search to first 'n' turnover time scales
    it = find_approx(ndtime, 10);
    ind0 = length(ndtime); %find(vel(it:end) > 0, 1, 'first');
    if strcmpi(runs.name, 'ew-8352')
        ind0 = ind0 - 30;
    end
    [mn, imin] = min(vel(5:ind0));
    %[mn, imin] = min(vel(:));

    vv = vel(imin:end) - mn * factor;
    tind = find(vv > 0, 1, 'first');
    tind = tind + imin - 1;

    if tind >= (length(runs.eddy.my)-10)
        disp(['WARNING: ' runs.name ' - Location within ten timesteps ' ...
              'of end | factor = ' num2str(factor)]);
    end

    xx = runs.eddy.mx(tind);
    yy = runs.eddy.my(tind);
    cvy = vel(tind);

    runs.res.xx = xx;
    runs.res.yy = yy;
    runs.res.tind = tind;
    runs.res.cvy = cvy;

    if debug_plot
        figure;
        insertAnnotation('runArray.locate_resistance');
        subplot(221)
        plot(ndtime, runs.eddy.my(1:runs.eddy.tend-1));
        title(runs.name);
        linex(ndtime(tind));

        subplot(223)
        plot(ndtime, vel);
        linex(ndtime(imin)); liney(mn); liney(mn*factor);
        linex(ndtime(tind));

        subplot(2,2,[2 4])
        plot(runs.eddy.mx/1000, runs.eddy.my/1000);
        hold on; axis image
        plot(xx/1000, yy/1000, 'x', 'MarkerSize', 14);
    end
end