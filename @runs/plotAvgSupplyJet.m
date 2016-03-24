function [] = plotAvgSupplyJet(runs)

    supply = runs.supply;

    figure;
    insertAnnotation([runs.name '.plotAvgSupplyJet']);
    subplot(211);
    contourf(supply.ymat/1000 - runs.bathy.xsb/1000, supply.zmat, supply.asvmean, ...
            40, 'EdgeColor', 'none');
    center_colorbar;
    linex([supply.xmin supply.xmin-supply.xscale]/1000 - runs.bathy.xsb/1000)
    title(['Averaged along-shelf supply jet (m/s) | ' runs.name]);
    xlabel('Y - Y_{sb} (km)');
    ylabel('Z (m)');
    set(gcf, 'Position', [455 387 1013 710]);
    pbaspect([1.618 1 1]);
    beautify([22 24 28]);

    subplot(212)
    contourf(supply.ymat/1000 - runs.bathy.xsb/1000, supply.zmat, supply.csdmean, ...
            40, 'EdgeColor', 'none');
    caxis([0 runs.bathy.xsl]);
    title(['Dye | ' runs.name]);
    xlabel('Y - Y_{sb} (km)');
    ylabel('Z (m)');
    set(gcf, 'Position', [455 387 1013 710]);
    pbaspect([1.618 1 1]);
    beautify([22 24 28]);

end