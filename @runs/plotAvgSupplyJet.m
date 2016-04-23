function [] = plotAvgSupplyJet(runs)

    supply = runs.supply;

    xloc = -[supply.xscale supply.IntersectScale]/1000;

    figure;
    hax = packfig(3,1);
    insertAnnotation([runs.name '.plotAvgSupplyJet']);

    axes(hax(1));
    plot(supply.ymat(:,1)/1000 - runs.bathy.xsb/1000, supply.vint);
    hold on;
    plot(supply.ymat(:,1)/1000 - runs.bathy.xsb/1000, supply.shelf.vint);
    ylabel({'Depth-integrated'; 'velocity'});
    legend('All water', 'Shelf water', 'Location', 'NorthWest');
    title(['Averaged along-shelf supply jet (m/s) | ' runs.name]);
    beautify([22 24 28]);
    linex(xloc);
    liney(0);
    hax(1).XTickLabel = [];

    axes(hax(2));
    contourf(supply.ymat/1000 - runs.bathy.xsb/1000, supply.zmat, supply.vmean, ...
            40, 'EdgeColor', 'none');
    hcb(1) = center_colorbar;
    linex(xloc);
    ylabel('Z (m)');
    %pbaspect([1.618 1 1]);
    beautify([22 24 28]);
    hax(2).XTickLabel = [];

    axes(hax(3));
    contourf(supply.ymat/1000 - runs.bathy.xsb/1000, supply.zmat, ...
             supply.csdmean/1000 - runs.bathy.xsb/1000, 40, 'EdgeColor', 'none');
    hold on;
    contour(supply.ymat/1000 - runs.bathy.xsb/1000, supply.zmat, ...
            supply.csdmean/1000 - runs.bathy.xsb/1000, [1 1] * 0, 'k', ...
            'LineWidth', 2);
    hcb(2) = center_colorbar;
    text(0.1, 0.1, 'Averaged cross-shelf dye - Y_{sb} (km)', 'Units', 'normalized');
    xlabel('Y - Y_{sb} (km)');
    ylabel('Z (m)');
    linex(xloc);
    %pbaspect([1.618 1 1]);
    beautify([22 24 28]);

    linkaxes(hax, 'x');
    xlim([min(supply.ymat(:) - runs.bathy.xsb)/1000 ...
          max(supply.ymat(:) - runs.bathy.xsb)/1000]);

    set(gcf, 'Position', [455 387 1013 710]);
    hcb(1).Position(1) = 0.92;
    hcb(2).Position(1) = 0.92;
end