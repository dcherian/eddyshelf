% time series of eddy diagnostics
function [handles] = EddyDiags(runs, hfig)

    eddy = runs.eddy;
    tvec = runs.eddy.t;

    nsmooth = 10;

    itfit = runs.FitCenterVelocity;
    [itsl, itse, tsl, tse] = runs.getEddyCenterTimeScales;

    figure(hfig); maximize;
    hax = packfig(2,1);
    axes(hax(1));
    [axyy, hplt(1), hplt(2)] = ...
        plotyy(tvec, eddy.vor.dia/2000, ...
               tvec, eddy.Lgauss);
    linex(tvec([itsl itfit]));

    axyy(1).YLabel.String = 'Horizontal Scale (km)';
    axyy(2).YLabel.String = 'Vertical Scale, L^z (km)';
    axyy(1).XTickLabel = {};
    axes(axyy(1)); beautify;
    axes(axyy(2)); beautify;

    axes(hax(2));
    hold on;
    hplt(3) = plot(tvec, smooth(eddy.mvx, nsmooth));
    hplt(4) = plot(tvec, smooth(eddy.mvy, nsmooth));
    linex(tvec([itsl itfit]));
    liney(0);
    hleg = legend([hplt(4) hplt(3)], {'V_{cen}^y', 'V_{cen}^x'}, 'Location', 'NorthEast');
    hleg.Position(1) = 0.84;
    hleg.Position(2) = 0.25;
    ylabel({'Eddy center'; 'translation velocity';  '(m/s)'});
    xlabel('Time (days)');
    beautify;

    handles.hax = hax;
    handles.axyy = axyy;
    handles.hplt = hplt;
    handles.hleg = hleg;
end