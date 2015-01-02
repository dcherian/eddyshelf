% locate point of "resistance"
% defined as start of half the maximum southward speed.
function [xx,yy,tind] = locate_resistance(runs)

    debug_plot = 0;

    % number of points to smooth over.
    npts = 15;

    % smooth velocity a lot!
    if bathy.axis == 'y'
        vel = smooth(runs.eddy.mvy, npts);
    else
        vel = smooth(runs.eddy.mvx, npts);
    end

    % locate minimum, i.e., max. southward/westward speed
    [mn, imin] = min(vel);

    tind = find_approx(vel(imin:end), mn/2, 1);
    tind = tind + imin;

    xx = runs.eddy.mx(tind);
    yy = runs.eddy.my(tind);

    if debug_plot
        figure;
        subplot(211)
        plot(runs.eddy.my/1000);
        subplot(212)
        plot(vel);
        linex(imin); liney(mn); liney(mn/2);

        linex(tind);
        subplot(211); linex(tind);
    end
end