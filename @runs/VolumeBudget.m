function [handles] = VolumeBudget(runs)

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

    eddyeColor = runs.eddyeColormap;

    figure; maximize;
    insertAnnotation([runs.name '.VolumeBudget']);
    hold on;
    handles.totalsb = plot(runs.time/86400, intAS(isb,:)', 'Color', 'k', ...
                           'DisplayName', 'Cross-shelf: total');
    handles.shelfonly = plot(runs.csflux.time/86400, runs.csflux.off.slope(:,1,1), ...
                             'Color', runs.shelfSlopeColor('light'), 'LineStyle', '--', ...
                             'DisplayName', 'Cross-shelf: shelf water');
    handles.nonshelf = plot(runs.time/86400, -runs.csflux.off.slope(:,1,1) + intAS(isb,:)', ...
                            'Color',eddyeColor(end-3,:,:), 'LineStyle', '--', ...
                            'DisplayName', 'Cross-shelf: eddy/slope water');
    handles.hisponge = plot(runs.time/86400, intCS(2,:)', ...
                            'DisplayName', 'Along-shelf: high sponge');
    handles.lowsponge = plot(runs.time/86400, -1*intCS(1,:)', ...
                             'DisplayName', 'Along-shelf: low sponge');
    handles.residual = plot(avg1(runs.time/86400), ...
                            intZT - avg1(intAS(isb,:) + intCS(2,:) - intCS(1,:))', ...
                            'Color', [1 1 1]*0.65, ...
                            'DisplayName', 'Budget residual');
    xlabel('Time (days)');
    ylabel('Total volume flux (m^3/s)');
    xlim([0 max(runs.time/86400)]);
    hax = gca;
    htxt(1) = text(0.2,0.2, 'domain gains water', 'Units', 'Normalized');
    htxt(2) = text(0.2, 0.8, 'domain loses water', 'Units', 'Normalized');
    htxt(1).Units = 'data'; htxt(1).Position(2) = hax.YTick(find(hax.YTick < 0, 1, 'last'))*1.5;
    htxt(2).Units = 'data'; htxt(2).Position(2) = hax.YTick(find(hax.YTick > 0, 1, 'first'))*1.5;
    linkprop(htxt,{'FontSize', 'Color'});
    htxt(1).FontSize = 20;
    htxt(1).Color = [1 1 1]*0.55;

    hax.XAxisLocation = 'origin';
    hleg = legend('show', 'Location', 'NorthWest');
    title([runs.name ' | Volume budget']);
    beautify([22 24 28]);
    hax.XLabel.Position(1) = max(runs.time/86400);

    handles.htxt = htxt;
    handles.hax = hax;
    handles.hleg = hleg;
end