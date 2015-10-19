% Diagnose peak location and peak width from vertical profile of
% streamer.
% Best is to look for depth where slope of profile drops to 1/2 its
% max value.
function [zmax, imax, zwidth, idiff] = streamer_peak(runs, isobath, debug)

    if ~exist('debug', 'var'), debug = 0; end
    nsmooth = 3;

    kk = 1;

    zmax = NaN; imax = NaN; zwidth = NaN;
    if isobath > length(runs.csflux.x), return; end

    [~,maxloc] = runs.calc_maxflux(isobath);
    vec = smooth(runs.csflux.west.slopewater.vertitrans(:,isobath,isobath), ...
                 nsmooth);
    %vec = smooth(runs.csflux.west.slopezt(:,maxloc,isobath, isobath), ...
    %             nsmooth);
    zvec = runs.csflux.vertbins(:,isobath);
    dvec = smooth(diff(vec,1,1)./diff(zvec), nsmooth);
    dzvec = avg1(zvec);

    % peak
    [vmax,ivmax] = max(abs(vec));
    zmax(kk) = zvec(ivmax);

    % peakwidth based on surface value
    ibot = find_approx(vec(1:end-2), vec(end), 1);
    zbot(kk) = zvec(ibot);

    % peakwidth based on derivative drop to half it's max
    % value — look for first depth after maximum
    [dvmax, idmax] = max(dvec);
    idiff = find(dvec(1:idmax) - dvmax/2 < 0, 1, 'last');

    if isempty(idiff) & ~isempty(findstr(runs.name, 'ew-8'))
        % max slope is near bottom, so search up to location of actual maximum.
        %idiff = find(dvec(idmax+10:ivmax) - dvmax/2 < 0, 1, 'first');
        % just look for (peak value)/2
        idiff = find(vec - vmax/2 < 0, 1, 'last');
    end

    if ~isempty(idiff), zdiff(kk) = dzvec(idiff); end

    if debug
        figure;
        plot(vec, zvec);
        liney(zmax, 'max');
        liney(zbot(kk), 'bot');
        try
            liney(zdiff(kk), 'der');
        catch ME; end
        liney(runs.bathy.hsb*-1, 'hsb', 'k');
        linex(0);
        plot(diff(vec,1,1)./diff(zvec)*10, avg1(zvec));
        title(runs.name); drawnow;
    end

    if isempty(idiff)
        warning('streamer_peak: Derivative criterion not satisfied.');
        return;
    end

    zwidth = zdiff;
end
