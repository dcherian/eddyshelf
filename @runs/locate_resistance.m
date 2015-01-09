% locate point of "resistance"
% defined as start of half the maximum southward speed.
% returns center locations and time index
function [xx,yy,tind] = locate_resistance(runs)

    debug_plot = 0;

    % number of points to smooth over.
    npts = 15;

    % smooth velocity a lot!
    if runs.bathy.axis == 'y'
        vel = smooth(runs.eddy.mvy, npts);
    else
        vel = smooth(runs.eddy.mvx, npts);
    end

    % locate minimum, i.e., max. southward/westward speed
    [mn, imin] = min(vel);

    vv = vel(imin:end) - mn * 1/2;
    tind = find(vv > 0, 1, 'first');
    tind = tind + imin;

    xx = runs.eddy.mx(tind);
    yy = runs.eddy.my(tind);

    if debug_plot
        figure;
        subplot(221)
        plot(runs.eddy.my/1000);
        title(runs.name);
        linex(tind);

        subplot(223)
        plot(vel);
        linex(imin); liney(mn); liney(mn/2);
        linex(tind);

        subplot(2,2,[2 4])
        plot(runs.eddy.mx/1000, runs.eddy.my/1000);
        hold on; axis image
        plot(xx/1000, yy/1000, 'x', 'MarkerSize', 14);
    end

end