% Diagnose peak location and peak width from vertical profile of
% streamer.
% Best is to look for depth where slope of profile drops to 1/2 its
% max value.
function [zmax, ivmax, zwidth, idiff] = streamer_peak(runs, isobath, debug)

    if ~exist('debug', 'var'), debug = 0; end

    nsmooth = 3;

    if isobath > length(runs.csflux.x), return; end

    [~,maxloc] = runs.calc_maxflux(isobath);
    vec = smooth(runs.csflux.off.slopewater.vertitrans(:,isobath,isobath), ...
                 nsmooth);
    %vec = smooth(runs.csflux.off.slopezt(:,maxloc,isobath, isobath), ...
    %             nsmooth);
    zvec = runs.csflux.vertbins(:,isobath);

    [vmax,ivmax] = max(abs(vec));
    zmax = zvec(ivmax);

    % for z-peak, idifflo is what I'm interested in since idiffhi is the surface grid cell
    [~,idiff,~] = findProfilePeakWidth(vec, zvec, debug);
    zwidth = zvec(idiff);

    if debug
        liney(runs.bathy.hsb*-1, 'hsb', 'k');
        title(runs.name);
    end
end
