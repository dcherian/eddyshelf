function [hplt] = plotBottomRhoPV(runs, tindices)

    for tindex = flip(tindices)
        ind = find(tindex == tindices);
        if length(tindices) > 1
            hax = subplot(2,ceil(length(tindices)/2),ind);
        else
            hax = gca;
        end
        hplt = runs.overlay_section('pv', 'rho', tindex, {}, 's', 1, [], hax);

        insertAnnotation([runs.name '.plotBottomRhoPV']);
        if tindex == tindices(end)
            clim = caxis;
        end
        caxis([0 clim(end)]);
        warning off;
        colormap(cbrewer('seq', 'Reds', 40));
        warning on;

        xlim([-1 1]*1.3*runs.eddy.rhovor.dia(tindex)/2000 ...
             + runs.eddy.mx(tindex)/1000);
        ylim([runs.bathy.xsb/1000-10 runs.eddy.my(tindex)/1000]);

        hplt.h_rho.LineWidth = 2;

        beautify;
    end
end