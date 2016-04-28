% time series of eddy diagnostics
function [handles] = EddyDiags(runs, hfig)

    eddy = runs.eddy;
    tvec = runs.eddy.t;

    nsmooth = 10;

    mvx = smooth(eddy.mvx, nsmooth);
    mvy = smooth(eddy.mvy, nsmooth);

    [itfit,~,~,~,FitTimeSeries] = runs.FitCenterVelocity; itfit = itfit(1);
    [itsl, itse, tsl, tse] = runs.getEddyCenterTimeScales;

    figure(hfig); maximize;
    hax = packfig(2,1);
    axes(hax(1));
    [axyy, hplt(1), hplt(2)] = ...
        plotyy(tvec, eddy.vor.dia/2000, ...
               tvec, eddy.Lgauss);
    axes(axyy(1)); hold on;
    handles.hsl(1) = plot(tvec(itsl), eddy.vor.dia(itsl)/2000, ...
                          'ko', 'MarkerSize', 12);
    handles.hfit(1) = plot(tvec(itfit), eddy.vor.dia(itfit)/2000, ...
                           'kx', 'MarkerSize', 12);
    handles.hlabel(1) = text(0.05,0.1,'a)','Units','normalized');
    %linex(tvec([itsl itfit]));
    axes(axyy(2)); hold on
    handles.hsl(2) = plot(tvec(itsl), eddy.Lgauss(itsl), 'ko', 'MarkerSize', 12);
    handles.hfit(2) = plot(tvec(itfit), eddy.Lgauss(itfit), 'kx', 'MarkerSize', 12);

    axyy(1).YLabel.String = 'Horizontal Scale, L_0 (km)';
    axyy(2).YLabel.String = 'Vertical Scale, L_z (m)';
    axyy(1).XTickLabel = {};
    axes(axyy(1)); beautify;
    axes(axyy(2)); beautify;

    axes(hax(2));
    hold on;
    hplt(3) = plot(tvec, mvx);
    hplt(4) = plot(tvec, mvy);
    handles.hsl(3) = plot(tvec([itsl itsl]), [mvx(itsl) mvy(itsl)], ...
                          'ko', 'MarkerSize', 12);
    handles.hfit(3) = plot(tvec([itfit itfit]), [mvx(itfit) mvy(itfit)], ...
                           'kx', 'MarkerSize', 12);
    hplt(5) = plot(FitTimeSeries(:,1)/86400, FitTimeSeries(:,2), 'k--');
    %linex(tvec([itsl itfit]));
    liney(0);
    hleg = legend([hplt(4) hplt(3) hplt(5)], {'V_{cen}^y', 'V_{cen}^x', 'Fit'}, ...
                  'Location', 'NorthEast');
    hleg.Position(1) = 0.86;
    hleg.Position(2) = 0.23;
    ylabel({'Eddy center'; 'translation velocity';  '(km/day)'});
    xlabel('Time (days)');
    handles.hlabel(2) = text(0.05,0.1,'b)','Units','normalized');
    beautify;

    handles.hax = hax;
    handles.axyy = axyy;
    handles.hplt = hplt;
    handles.hleg = hleg;
end