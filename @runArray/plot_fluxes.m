function [] = plot_fluxes(runArray)
    hfig1 = figure;
    subplot(2,1,1); hold all
    subplot(2,1,2); hold all

    hfig2 = figure;
    hold all
    % subplot(2,2,1); hold all
    % subplot(2,2,2); hold all
    % subplot(2,2,3); hold all
    % subplot(2,2,4); hold all

    hfig3 = figure;
    hold all;

    hfig4 = figure;
    subplot(2,1,1); hold all
    subplot(2,1,2); hold all

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

        run = runArray.array(ii);
        names{ff} = getname(runArray, ii);

        if run.params.flags.flat_bottom
            continue;
        end

        tind = find_approx(run.ndtime, 1.5, 1);
        Ue = run.eddy.V(tind);
        He = run.bathy.hsb;
        Le = run.eddy.vor.dia(tind)/2;

        fluxscl = 1e6;Ue * Le * He;
        transscl = 1;fluxscl * (2*Le/Ue);

        ndtime = run.ndtime; %run.csflux.time/run.tscale
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHELF WATER
        figure(hfig1)
        subplot(2,1,1)
        hgplt1(ff) = plot(ndtime(1:end-2), ...
                     smooth((run.csflux.west.shelf(1:end-2, ...
                                                   1))/fluxscl, 3));

        subplot(2,1,2)
        plot(ndtime, ...
             run.csflux.west.itrans.shelf(:,1)/transscl);

        % total transport
        ttrans = max(abs(run.csflux.west.itrans.shelf(:,1)));
        figure(hfig2)
        try
            % find location of center at t=1.5 (arbitrary
            % choice)
            %                    xnd = (run.csflux.west.shelfwater.bins)./run.rrshelf;
            %                    subplot(2,2,1)
            %hgplt = plot(xnd, run.csflux.west.shelfwater.itrans);

            %subplot(2,2,2)
            %plot(xnd, run.csflux.west.shelfwater.itrans ...
            %             ./ttrans);

            %subplot(2,2,3)
            %plot(run.csflux.west.shelfwater.vertitrans, ...
            %     run.csflux.west.shelfwater.vertbins);

            %subplot(2,2,4)
            profile = ...
                run.csflux.west.shelfwater.vertitrans./ ...
                ttrans;
            zvec = run.csflux.west.shelfwater.vertbins ./ ...
                   run.bathy.hsb;
            bc = baroclinicity(zvec, profile);
            hgplt2(ff) = plot(profile, zvec);
            names2{ff} =  [names{ff} ' | bc = ' num2str(bc,'%.3f')];
        catch ME
            disp(ME)
        end

        figure(hfig3)
        % change from envelope to depth
        env = run.csflux.west.shelfwater.envelope;
        env(isnan(env)) = max(env);
        ind = vecfind(run.rgrid.y_rho(:,1), env);
        metric = run.bathy.h(1,ind)./run.bathy.hsb .* ...
                 (1+run.rgrid.f(run.bathy.isb,1)./run.rgrid.f(ind,1))';

        hgplt3(ff) = plot(ndtime, ...
                          (run.bathy.xsb - env)./run.rrshelf);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDDY WATER
        figure(hfig4)
        subplot(2,1,1)
        hgplt4(ff) = plot(ndtime(1:end-2), ...
                     smooth((run.csflux.east.eddy(1:end-2, ...
                                                  1))/fluxscl, 3));

        subplot(2,1,2)
        plot(ndtime, ...
             run.csflux.east.itrans.eddy(:,1)/transscl);
    end

    figure(hfig1)
    subplot(2,1,1)
    ylabel('Flux of shelf water (Sv)');
    liney(0, [], [1 1 1]*0.75);
    subplot(2,1,2)
    ylabel('Total volume transported');
    xlabel('Non-dimensional time');
    legend(hgplt1, names);

    figure(hfig2)
    %subplot(2,2,1)
    %ylabel('Total volume transported (m^3)');
    %xlabel('cs location (km)');

    %subplot(2,2,2)
    %limy = ylim;
    %ldefbt = sqrt(9.81 * run.bathy.hsb)/run.params.phys.f0;
    %xnd = xnd .* run.rrshelf ./ ldefbt;
    %plot(xnd(2:end-1), 1./xnd(2:end-1).^2, 'color', [1 1 ...
    %                    1]*1);
    %ylim(limy)
    %ylabel('Normalized volume transported');
    %xlabel('cs location / shelfbreak rossby radius');

    %subplot(2,2,3)
    %xlabel('Total volume transported (m^3)');
    %ylabel('Vertical bin (m)');

    %subplot(2,2,4)
    ylim([-1 0]);
    limx = xlim;
    xlim([0 limx(2)]);
    xlabel('Normalized volume transported');
    ylabel('Vertical bin / Shelfbreak depth');
    %legend(hgplt2, names2, 'Location', 'SouthWest');

    figure(hfig3)
    xlabel('Non-dimensional time');
    ylabel('isobath most on-shore source of water / h_{sb}');
    legend(hgplt3, names{ff});

    figure(hfig4)
    subplot(2,1,1);
    ylabel('Flux of eddy water (Sv)');
    legend(hgplt4, names{ff});
    subplot(2,1,2);
    ylabel('Total volume transported (m^3)');
    xlabel('Non-dimensional time');
end
