function [hplt] = plotBottomRhoPV(runs, tindex)

    hplt = runs.overlay_section('pv', 'rho', tindex, {}, 's', 1);

    insertAnnotation([runs.name '.plotBottomRhoPV']);
    clim = caxis;
    caxis([0 clim(end)]);
    warning off;
    colormap(cbrewer('seq', 'Reds', 40));
    warning on;

    xlim([-1 1]*1.3*runs.eddy.rhovor.dia(tindex)/2000 ...
         + runs.eddy.mx(tindex)/1000);
    ylim([-1.5 -0.25]*runs.eddy.rhovor.dia(tindex)/2000 ...
         + runs.eddy.my(tindex)/1000);

    hplt.h_rho.LineWidth = 2;
end