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
            indsb = 1;
        else
            sy1 = 1;
            sy2 = runs.bathy.isb;
            indsb = 2;
        end
        sx1 = sx1+10; sx2 = sx2-10;

        asvelz(1,:,:,:) = dc_roms_read_data(runs, 'u', [], {'x' sx1 sx1; 'y' 1 runs.bathy.isb});
        asvelz(2,:,:,:) = dc_roms_read_data(runs, 'u', [], {'x' sx2 sx2; 'y' 1 runs.bathy.isb});
        csdyez(1,:,:,:) = dc_roms_read_data(runs, runs.csdname, [], {'x' sx1 sx1; 'y' 1 ...
                            runs.bathy.isb});
        csdyez(2,:,:,:) = dc_roms_read_data(runs, runs.csdname, [], {'x' sx2 sx2; 'y' 1 ...
                            runs.bathy.isb});

        dz = permute(diff(runs.rgrid.z_w(:,1:runs.bathy.isb,sx2),1,1), [3 2 1]);
        shasvel = squeeze(sum(bsxfun(@times, asvelz .* (csdyez <= runs.bathy.xsb), dz), 3));

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
        indsb = 2;

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

    intCSsh = squeeze(trapz(csvec, shasvel, csax));
    intAS = squeeze(trapz(asvec, csvel, asax)); % integrated along-shore
    intCS = squeeze(trapz(csvec, asvel, csax)); % integrated cross-shore
    intZT = squeeze(trapz(asvec, trapz(csvec, zt, csax), asax));

    eddyeColor = runs.eddyeColormap;

    figure; maximize;
    insertAnnotation([runs.name '.VolumeBudget']);
    hold on;
    handles.totalsb = plot(runs.time/86400, intAS(indsb,:)', 'Color', 'k', ...
                           'DisplayName', 'Cross-shelf: total');
    handles.shelfonly = plot(runs.csflux.time/86400, runs.csflux.off.slope(:,1,1), ...
                             'Color', runs.shelfSlopeColor('light'), 'LineStyle', '--', ...
                             'DisplayName', 'Cross-shelf: shelf water');
    handles.nonshelf = plot(runs.time/86400, -runs.csflux.off.slope(:,1,1) + intAS(indsb,:)', ...
                            'Color',eddyeColor(end-3,:,:), 'LineStyle', '--', ...
                            'DisplayName', 'Cross-shelf: eddy/slope water');
    handles.hisponge = plot(runs.time/86400, intCS(2,:)', ...
                            'DisplayName', 'Along-shelf: high sponge');
    handles.hispongesh = plot(runs.time/86400, intCSsh(2,:)', ...
                              'Color', runs.shelfSlopeColor('light'), 'LineStyle', ':', ...
                              'DisplayName', 'Along-shelf: high sponge: shelf water');
    handles.lowsponge = plot(runs.time/86400, -1*intCS(1,:)', ...
                             'DisplayName', 'Along-shelf: low sponge');
    handles.residual = plot(avg1(runs.time/86400), ...
                            intZT - avg1(intAS(indsb,:) + intCS(2,:) - intCS(1,:))', ...
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