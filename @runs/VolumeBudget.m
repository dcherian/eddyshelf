function [] = VolumeBudget(runs)

    runs.read_velbar;
    runs.read_zeta;

    use runs.spng

    ubar = avg1(runs.ubar(:,2:end-1,:), 1);
    vbar = avg1(runs.vbar(2:end-1,:,:), 2);

    xrho = runs.rgrid.x_rho(2:end-1,2:end-1);
    yrho = runs.rgrid.y_rho(2:end-1,2:end-1);

    if runs.bathy.axis == 'y'
        if runs.sgntamp < 1
            sy1 = runs.bathy.isb;
            sy2 = size(ubar, 2);
            isb = 1;
        else
            sy1 = 1;
            sy2 = runs.bathy.isb;
            isb = 2;
        end
        sx1 = sx1+2; sx2 = sx2-2;

        asvel = bsxfun(@times, runs.ubar([sx1 sx2], sy1:sy2, :), ...
                       runs.bathy.h([sx1 sx2], sy1:sy2));
        csvel = bsxfun(@times, runs.vbar([sx1:sx2], [sy1 sy2], :), ...
                       runs.bathy.h(sx1:sx2, [sy1 sy2]));

        asvec = xrho(1,sx1:sx2)';
        csvec = yrho(sy1:sy2,1)';

        csax = 2; asax = 1;
    else
        sx1 = 1;
        sx2 = runs.bathy.isb;
        isb = 2;

        sy1 = sy1 + 2; sy2 = sy2 - 2;

        csvel = bsxfun(@times, runs.ubar([sx1 sx2], sy1:sy2, :), ...
                       runs.bathy.h([sx1 sx2], sy1:sy2));
        asvel = bsxfun(@times, runs.vbar([sx1:sx2], [sy1 sy2], :), ...
                       runs.bathy.h(sx1:sx2, [sy1 sy2]));

        asvec = yrho(sy1:sy2,1)';
        csvec = xrho(1,sx1:sx2);

        csax = 1; asax = 2;
    end
    zt = bsxfun(@rdivide, diff(runs.zeta(sx1:sx2,sy1:sy2,:), 1, 3), ...
                    permute(diff(runs.time), [3 2 1]));

    intAS = squeeze(trapz(asvec, csvel, asax));
    intCS = squeeze(trapz(csvec, asvel, csax));
    intZT = squeeze(trapz(asvec, trapz(csvec, zt, csax), asax));

    figure; maximize;
    insertAnnotation([runs.name '.VolumeBudget']);
    hold on;
    plot(runs.csflux.time/86400, runs.csflux.off.slope(:,1,1), ...
         'Color', [1 1 1]*0.75);
    plot(runs.time/86400, intAS(isb,:)');
    plot(runs.time/86400, -1*intCS(1,:)');
    plot(runs.time/86400, intCS(2,:)');
    plot(avg1(runs.time/86400), ...
         intZT - avg1(intAS(isb,:) + intCS(2,:) - intCS(1,:))', 'k');
    plot(runs.time/86400, -runs.csflux.off.slope(:,1,1) + intAS(isb,:)', ...
         'Color', [1 1 1]*0.75, 'LineStyle', '--');
    xlabel('Time (days)');
    ylabel('Total volume flux (m^3/s)');
    xlim([0 max(runs.time/86400)]);
    htxt(1) = text(0.5,0.1, 'Import', 'Units', 'Normalized');
    htxt(2) = text(0.5, 0.9, 'Export', 'Units', 'Normalized');
    linkprop(htxt,{'FontSize', 'Color'});
    htxt(1).FontSize = 20;
    htxt(1).Color = [1 1 1]*0.55;
    hax = gca;
    hax.XAxisLocation = 'origin';
    hax.XLabel.Position(1) = max(runs.time/86400);
    legend('Offshore flux of shelf water', 'At shelfbreak', ...
           'At sponge (low)', 'At sponge (high)', 'budget residual', ...
           'Non-shelf water flow at shelfbreak', 'Location', 'NorthWest');
    title([runs.name ' | Volume budget']);
    beautify;

end