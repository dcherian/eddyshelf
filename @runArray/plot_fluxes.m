function [] = plot_fluxes(runArray, isobath, source, factor, figs)

    if ~exist('isobath', 'var') | isempty(isobath), isobath = 2; end
    if ~exist('source', 'var') | isempty(source), source = isobath; end
    if ~exist('factor', 'var') | isempty(factor), factor = 2; end
    if ~exist('figs', 'var'), figs = zeros([10 1]); figs(1) = 1; end

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    corder_backup = runArray.sorted_colors;

    useAvgStreamer = 1;

    figs(length(figs)+1:10) = 0;

    if figs(1)
        hfig1 = figure; % shelf water flux time series
        ax1 = packfig(2,1, 'rows');
        axes(ax1(1)); hold all; axes(ax1(2)); hold all;
    end

    if figs(2)
        hfig2 = figure; hold all % integrated vertical structure / baroclinicity
    end

    if figs(3)
        hfig3 = figure; hold all; % envelope
    end

    if figs(4)
        figure; subplot(2,1,1); hold all; % eddy water flux
        subplot(2,1,2); hold all;
    end

    if figs(5)
        hfig5 = figure; hold on;% along shelf structure
    end

    if figs(6)
        hfig6 = figure; hold on; % streamer velocity
    end

    if figs(7)
        hfig7 = figure; hold on; % cross-shore bins
    end

    if figs(8)
        hfig8 = figure; hold on; % streamer max. vel.
    end

    if figs(9)
        hfig9 = figure; hold all;% instantaneous streamer vertical structure
    end

    if figs(10)
        hfig10 = figure; subplot(211); hold all; % horizontal velocity profile at surface
        subplot(212); hold all;
    end

    nsmooth = 1;

    nn = 1;
    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

        run = runArray.array(ii);

        if isempty(run.csflux) | (isobath > length(run.csflux.x))
            disp(['Skipping ' run.name]);
            continue;
        end

        names{nn} = getname(runArray, ii); nn = nn + 1;

        locstr = num2str(run.csflux.ndloc(isobath), 2);

        tind = find_approx(run.ndtime, 1.5, 1);
        Lz = run.eddy.Lgauss(1);
        hsb = run.bathy.hsb;

        He = hsb;
        Le = run.eddy.vor.dia(1)/2;

        fluxscl = 1000; %run.eddyfluxscale;
        transscl = 1; He * Le^2;

        %fluxvec = run.csflux.off.slope(:,isobath, source);
        fluxvec = run.recalculateFlux(-factor*hsb,isobath);
        [maxf, maxi] = run.calc_maxflux(fluxvec);
        [start,stop] = run.flux_tindices(fluxvec);
        ifluxvec = run.csflux.off.itrans.slope(:,isobath, source);
        ttrans = max(abs(ifluxvec));
        ttransv = abs(trapz(run.csflux.vertbins(:,isobath), ...
                        run.csflux.off.slopewater.vertitrans(:,isobath,source)));

        % ndtime = run.csflux.time/run.eddy.turnover;
        ndtime = run.csflux.time/86400;

        R = run.csflux.R;
        [~,~,restind] = run.locate_resistance;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHELF WATER
        if figs(1)
            figure(hfig1)
            axes(ax1(1));
            hgplt1(ff) = plot(ndtime, smooth(fluxvec/fluxscl, nsmooth));
            plot(ndtime(maxi), maxf/fluxscl, 'x', 'Color', hgplt1(ff).Color);
            axes(ax1(2));
            plot(ndtime, ifluxvec/transscl, 'Color', hgplt1(ff).Color);
            plot(ndtime([start stop]), ifluxvec([start stop])/transscl, ...
                 'o', 'Color', hgplt1(ff).Color, 'MarkerSize', 8);
        else
            colors = runArray.sorted_colors;%get(gca, 'ColorOrder');
            hgplt1(ff).Color = colors(ff,:);
        end

        %%%%%% BAROCLINICITY
        if figs(2)
            figure(hfig2)
            if useAvgStreamer & ~isnan(run.streamer.off.xwidth(isobath))
                profile = run.streamer.off.zprof(:,isobath);
                vertbins = run.streamer.zvec(:,isobath);
            else
                profile = ...
                    run.csflux.off.slopewater.vertitrans(:,isobath,source);
                vertbins = run.csflux.vertbins(:,isobath);
            end
            zvec = vertbins ./ hsb; %max(abs(vertbins));
            zind = find_approx(vertbins, -1*hsb);

            %[~,~,zwidth,~] = run.streamer_peak(isobath);
            %if ~isnan(zwidth)
            %    zind = [zind find_approx(vertbins, zwidth)];
            %end

            bc = baroclinicity(zvec, profile);
            profile = profile ./ max(profile);
            hgplt2(ff) = plot(profile, zvec, 'Color',hgplt1(ff).Color);
            %plot(profile(zind), zvec(zind), 'x', 'Color', hgplt1(ff).Color);
            names2{ff} =  [names{ff} ' | bc = ' num2str(bc, '%.3f')];
        end

        %%%%%%%%%%%%% ENVELOPE
        if figs(3)
            figure(hfig3)
            % change from envelope to depth
            env = run.csflux.off.slopewater.envelope(:,isobath,1);
            env(isnan(env)) = max(env);
            ind = vecfind(run.rgrid.y_rho(:,1), env);
            %metric = run.bathy.h(1,ind)./run.bathy.hsb .* ...
            %         (1+run.rgrid.f(run.bathy.isb,1)./run.rgrid.f(ind,1))';

            hgplt3(ff) = plot(ndtime, (run.csflux.x(isobath) - env)/1000, ...
                              'Color', hgplt1(ff).Color);
            plot(0, run.bathy.Lbetash/1000, 'x', ...
                 'MarkerSize', 20, 'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDDY WATER
        if figs(4)
            figure(hfig4)
            subplot(2,1,1)
            hgplt4(ff) = plot(ndtime, ...
                              smooth((run.csflux.on.eddy(:, isobath))/fluxscl, ...
                                     nsmooth), 'Color', hgplt1(ff).Color);

            subplot(2,1,2)
            plot(ndtime, ...
                 run.csflux.on.itrans.eddy(:,isobath)/transscl, ...
                 'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%% ALONG SHELF STRUCTURE
        if figs(5)
            figure(hfig5)
            vmean = run.streamer.off.vmean(:,:,isobath);
            xivec = run.streamer.xivec * 1000;
            zvec = run.streamer.zvec(:,isobath);
            L = run.eddy.rhovor.dia(1)/2;
            ivmean = squeeze(trapz(zvec, vmean, 2));

            hgplt5(ff) = plot(xivec./L, ivmean, 'Color', ...
                              hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%% STREAMER VELOCITY
        if figs(6)
            figure(hfig6);

            yoR = run.csflux.ndloc; % y/R - used in csflux

            clear vnd y0oL
            for kk=1:length(yoR)
                [start, stop] = ...
                    run.flux_tindices(run.csflux.off.slopewater.vmax(:,kk));
                [vsmax,maxind] = ...
                    max(run.csflux.off.slopewater.vmax(start:stop,kk), [],1);
                % swirl velocity at ζ = 0 contour
                vscale = sqrt(2/exp(1))*run.eddy.V(start+maxind-1);
                vnd(kk) = vsmax/vscale;

                L = run.eddy.vor.dia(1)/2;
                y0oL(kk) =  R/L * (yoR(kk)-1); % y0/L - used in derivation
            end
            color =  hgplt1(ff).Color;[1 1 1]*0.75;
            hgplt6(ff) = plot(y0oL, vnd, 'Color', color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%% CROSS-SHORE BINS
        if figs(7)
            figure(hfig7)
            xsb = run.bathy.xsb;

            bins = avg1(run.csflux.off.bins{isobath}, 2);
            itrans = run.csflux.off.slopewater.itrans{isobath};
            hgplt7(ff) = plot((bins-xsb)./R, itrans./ttrans);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%% STREAMER MAX VEL
        if figs(8)
            figure(hfig8)
            [start,~] = run.flux_tindices(run.csflux.off.slope(:,isobath,isobath));
            hgplt8(ff) = plot(ndtime, ...
                              run.csflux.off.slopewater.vmax(:,isobath) ...
                              / run.eddy.V(1), 'Color', hgplt1(ff).Color);
            plot(ndtime(maxi), run.csflux.off.slopewater.vmax(maxi, isobath) ...
                 / run.eddy.V(1), 'x', 'Color', hgplt1(ff).Color);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%% INSTANTANEOUS VERTICAL STRUCTURE
        if figs(9)
            figure(hfig9)

            [~,tind] = run.calc_maxflux(run.csflux.off.slope(:,isobath,source));
            zvec = run.csflux.vertbins(:,isobath)./run.bathy.hsb;
            actual = run.csflux.off.slopezt(:,tind,isobath,source);
            actual = actual./max(actual);
            ideal = run.streamer_ideal_profile(isobath);
            hgplt9(ff) = plot(actual, zvec);
            plot(max(actual)*ideal, zvec, 'k--');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% velocity at shelfbreak
        if figs(10)
            tindex = 'max flux';
            figure(hfig10)
            ax10(1) = subplot(211);
            hgplt10(ff) = run.plot_profilex([run.csvelname 'bar'], tindex, run.bathy.axis, 'sb', gca);
            ax10(2) = subplot(212);
            hh = run.plot_profilex([run.asvelname 'bar'], tindex, run.bathy.axis, 'sb', gca);
            hh.Color = hgplt10(ff).Color;
        end
    end

    if figs(1)
        pba = [1.9 1 1];

        figure(hfig1)
        insertAnnotation('runArray.plot_fluxes');
        axes(ax1(1));
        title(['Slope water | ND isobath = ', locstr]);
        ylabel('Flux (mSv)');
        ylim([0 max(ylim)]);
        liney(0, [], [1 1 1]*0.75);
        legend(hgplt1, names);
        pbaspect(pba);
        beautify;

        axes(ax1(2));
        ylabel('Volume transported (m^3)');
        xlabel('Time (days)');
        ylim([0 max(ylim)]);
        pbaspect(pba);
        beautify;

        % positioning
        ax1(1).Position(2) = ax1(1).Position(2) + 0.01;
        ax1(1).Position(4) = ax1(1).Position(4) - 0.01;
        ax1(2).Position(2) = ax1(2).Position(2) - 0.01;
        ax1(2).Position(4) = ax1(2).Position(4) - 0.01;
    end

    if figs(2)
        figure(hfig2)
        set(gca, 'XAxisLocation', 'Top');
        insertAnnotation('runArray.plot_fluxes');
        liney(-1);%ylim([-1 0]);
        limx = xlim;
        xlim([0 limx(2)]);
        xlabel('Normalized area transported');
        ylabel('Z / H_{sb}');
        text(0.1, 0.1, ['y/R = ' num2str(run.csflux.ndloc(isobath))], ...
             'Units', 'normalized');
        legend(hgplt2, names2, 'Location', 'NorthWest');
        beautify;
    end

    if figs(3)
        figure(hfig3)
        insertAnnotation('runArray.plot_fluxes');
        xlabel('Non-dimensional time');
        ylabel('Distance from shelfbreak (km)');
        legend(hgplt3, names);
        beautify;
    end

    if figs(4)
        figure(hfig4)
        insertAnnotation('runArray.plot_fluxes');
        subplot(2,1,1);
        title(['Eddy water | ND isobath = ', locstr]);
        ylabel('Flux / Eddy flux above sb');
        liney(0);
        legend(hgplt4, names);
        beautify;

        subplot(2,1,2);
        ylabel('Vol_{transport}/Vol_{eddy}');
        xlabel('Non-dimensional time');
        liney(0);
        beautify;
    end

    if figs(5)
        figure(hfig5)
        insertAnnotation('runArray.plot_fluxes');
        xlabel('(X - X_e)/L_x');
        ylabel('Normalized flux');
        title('\int v(x,z,t) dz dt');
        linex([-1 0 1]); liney(0);
        xlim([-1 1]*6);
        legend(hgplt5, names);
        beautify;
    end

    if figs(6)
        figure(hfig6)
        insertAnnotation('runArray.plot_fluxes');
        ylabel('Max shelf water velocity / Eddy velocity (vor=0)');
        xlabel('y_0/L');
        liney(0);
        % legend(hgplt6, names);

        yL = linspace(min(xlim),max(xlim), 40);
        amp = 0.8;
        da = 0.15;
        hp = plot(yL, bsxfun(@times, [amp amp-da amp+da]', exp(-yL.^2)), ...
                  'Color', 'k', 'LineStyle', '-');
        legend(hp, [num2str(amp) '\pm' num2str(da) ' e^{-(y/' ...
                            'l)^2}']);
        ylim([0 1]);
        beautify;
    end

    if figs(7)
        figure(hfig7)
        insertAnnotation('runArray.plot_fluxes');
        linex(0, 'shelfbreak');
        linex(str2double(locstr), locstr);
        legend(hgplt7, names, 'Location', 'NorthWest');
        ylabel('Transport');
        xlabel('Source (y/R)');
        beautify;
    end

    if figs(8)
        figure(hfig8)
        legend(hgplt8, names, 'Location', 'NorthWest');
        xlabel('Time (days)');
        ylabel('Max streamer velocity / eddy velocity');
        beautify;
    end

    if figs(9)
        figure(hfig9)
        legend(hgplt9, names, 'Location', 'SouthEast');
        xlabel('Instantaneous Flux (m^3/s)');
        ylabel('Depth / H_{sb}');
        title('At instant of max flux');
        liney(-1);
        xlim([0 max(xlim)]);
        beautify;
    end

    if figs(10)
        figure(hfig10)
        subplot(211);
        legend(hgplt10, names, 'Location', 'NorthEast');
        ylabel('Cross-shelf velocity (m/s)');
        linex(0); liney(0);
        beautify;

        subplot(212);
        ylabel('Along-shelf velocity (m/s)');
        linex(0); liney(0);
        beautify;

        linkaxes(ax10, 'xy');
    end

    runArray.reset_colors(corder_backup);
end