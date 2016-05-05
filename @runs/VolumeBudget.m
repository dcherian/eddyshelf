function [handles] = VolumeBudget(runs)

    runs.read_velbar;
    runs.read_zeta;

    use runs.spng

    ubar = avg1(runs.ubar(:,2:end-1,:), 1);
    vbar = avg1(runs.vbar(2:end-1,:,:), 2);

    xrho = runs.rgrid.x_rho(2:end-1,2:end-1);
    yrho = runs.rgrid.y_rho(2:end-1,2:end-1);

    tind = 1:size(runs.csflux.off.slope, 1);
    tvec = runs.csflux.time/86400;

    if runs.bathy.axis == 'y'
        if runs.sgntamp < 1
            sy1 = runs.bathy.isb;
            sy2 = size(ubar, 2);
            indsb = 1;
        else
            sy1 = 1;
            sy2 = runs.bathy.isb;
            indsb = 2;
        end
        sx1 = sx1+2; sx2 = sx2-2;

        asvelz(1,:,:,:) = dc_roms_read_data(runs, 'u', [tind(1) tind(end)], ...
                                            {'x' sx1 sx1; 'y' 1 runs.bathy.isb});
        asvelz(2,:,:,:) = dc_roms_read_data(runs, 'u', [tind(1) tind(end)], ...
                                            {'x' sx2 sx2; 'y' 1 runs.bathy.isb});
        csdyez(1,:,:,:) = dc_roms_read_data(runs, runs.csdname, [tind(1) tind(end)], ...
                                            {'x' sx1 sx1; 'y' 1 runs.bathy.isb});
        csdyez(2,:,:,:) = dc_roms_read_data(runs, runs.csdname, [tind(1) tind(end)], ...
                                            {'x' sx2 sx2; 'y' 1 runs.bathy.isb});

        dz = permute(diff(runs.rgrid.z_w(:,1:runs.bathy.isb,sx2),1,1), [3 2 1]);
        shasvel = squeeze(sum(bsxfun(@times, asvelz .* (csdyez <= runs.bathy.xsb), dz), 3));

        asvel = bsxfun(@times, runs.ubar([sx1 sx2], sy1:sy2, tind), ...
                       runs.bathy.h([sx1 sx2], sy1:sy2));
        csvel = bsxfun(@times, runs.vbar([sx1:sx2], [sy1 sy2], tind), ...
                       runs.bathy.h(sx1:sx2, [sy1 sy2]));

        asvec = xrho(1,sx1:sx2)';
        csvec = yrho(sy1:sy2,1)';

        csax = 2; asax = 1;
    else
        sx1 = 1;
        sx2 = runs.bathy.isb;
        indsb = 2;

        sy1 = sy1 + 2; sy2 = sy2 - 2;

        csvel = bsxfun(@times, runs.ubar([sx1 sx2], sy1:sy2, tind), ...
                       runs.bathy.h([sx1 sx2], sy1:sy2));
        asvel = bsxfun(@times, runs.vbar([sx1:sx2], [sy1 sy2], tind), ...
                       runs.bathy.h(sx1:sx2, [sy1 sy2]));

        asvec = yrho(sy1:sy2,1)';
        csvec = xrho(1,sx1:sx2);

        csax = 1; asax = 2;
    end
    % zt = bsxfun(@rdivide, diff(runs.zeta(sx1:sx2,sy1:sy2,tind), 1, 3), ...
    %                 permute(diff(tvec*86400), [3 2 1]));

    VolAlongIsobathSh = squeeze(trapz(csvec, shasvel, csax));
    VolCrossIsobath = squeeze(trapz(asvec, csvel, asax)); % integrated along-shore
    VolAlongIsobath = squeeze(trapz(csvec, asvel, csax)); % integrated cross-shore
    %intZT = squeeze(trapz(asvec, trapz(csvec, zt, csax), asax));

    eddyeColor = runs.eddyeColormap;

    figure; maximize;
    insertAnnotation([runs.name '.VolumeBudget']);
    hold on;
    handles.totalsb = plot(tvec, VolCrossIsobath(indsb,:)', 'Color', 'k', ...
                           'DisplayName', 'Cross-shelf: total');
    handles.shelfonly = plot(runs.csflux.time/86400, runs.csflux.off.slope(:,1,1), ...
                             'Color', runs.shelfSlopeColor('light'), 'LineStyle', '-', ...
                             'DisplayName', 'Cross-shelf: shelf water');
    handles.nonshelf = plot(tvec, ...
                            -runs.csflux.off.slope(:,1,1) + VolCrossIsobath(indsb,:)', ...
                            'Color',eddyeColor(end-3,:,:), 'LineStyle', '-', ...
                            'DisplayName', 'Cross-shelf: eddy/slope water');
    handles.hisponge = plot(tvec, VolAlongIsobath(2,:)', 'k--', ...
                            'DisplayName', 'Along-shelf: high sponge: total');
    handles.hispongesh = plot(tvec, VolAlongIsobathSh(2,:)', ...
                              'Color', runs.shelfSlopeColor('light'), 'LineStyle', '--', ...
                              'DisplayName', 'Along-shelf: high sponge: shelf water');
    handles.lowsponge = plot(tvec, -1*VolAlongIsobath(1,:)', ...
                             'DisplayName', 'Along-shelf: low sponge');
    handles.residual = plot(tvec, ...
                            VolCrossIsobath(indsb,:) + ...
                            VolAlongIsobath(2,:) - VolAlongIsobath(1,:), ...
                            'Color', [1 1 1]*0.65, 'DisplayName', 'Budget residual');
    xlabel('Time (days)');
    ylabel('Total volume flux (m^3/s)');
    xlim([0 max(tvec)]);
    ylim([-1 1] * max(abs(ylim)));
    hax = gca;
    htxt(1) = text(0.7, 0.05, 'domain gains water', 'Units', 'Normalized');
    htxt(2) = text(0.7, 0.95, 'domain loses water', 'Units', 'Normalized');
    % htxt(1).Units = 'data';
    % htxt(1).Position(2) = hax.YTick(find(hax.YTick < 0, 1, 'last'))*1.5;
    % htxt(2).Units = 'data';
    % htxt(2).Position(2) = hax.YTick(find(hax.YTick > 0, 1, 'first'))*1.5;
    linkprop(htxt,{'FontSize', 'Color'});
    htxt(1).FontSize = 20;
    htxt(1).Color = [1 1 1]*0.55;

    hax.XAxisLocation = 'origin';
    hleg = legend('show', 'Location', 'SouthWest');
    title([runs.name ' | Volume budget']);
    beautify([22 24 28]);
    hax.XLabel.Position(1) = max(tvec);

    handles.htxt = htxt;
    handles.hax = hax;
    handles.hleg = hleg;

    % VolShelfWaterLost = runs.csflux.off.slope(:,1,1) + VolAlongIsobathSh(2,:)';
    % VolEddyWaterGained = runs.csflux.off.slope(:,1,1) - VolCrossIsobath(indsb,:)';

    % figure;
    % plot(tvec, VolShelfWaterLost-VolEddyWaterGained);
end