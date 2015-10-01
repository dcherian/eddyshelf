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
    dvec = diff(vec,1,1)./diff(zvec);
    dzvec = avg1(zvec);

    % peak
    [~,imax] = max(abs(vec));
    zmax(kk) = zvec(imax);

    % peakwidth based on surface value
    ibot = find_approx(vec(1:end-2), vec(end), 1);
    zbot(kk) = zvec(ibot);

    % peakwidth based on derivative drop to half it's max
    % value — look for first depth after maximum
    [dvmax, imax] = max(dvec);
    idiff = find(dvec(1:imax) - dvmax/2 < 0, 1, 'last');
    if isempty(idiff), return; end

    zdiff(kk) = dzvec(idiff);

    if debug
        figure;
        plot(vec, zvec);
        liney(zbot(kk), 'bot');
        liney(zdiff(kk), 'der');
        plot(diff(vec,1,1)./diff(zvec)*10, avg1(zvec));
        title(runs.name); drawnow;
    end
    kk = kk + 1;

    zwidth = zdiff;
end
