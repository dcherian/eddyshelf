function [] = plot_fluxes(runArray, isobath, source)

    if ~exist('isobath', 'var'), isobath = 2; end
    if ~exist('source', 'var'), source = 2; end

    hfig1 = []; %figure; ax1(1) = subplot(2,1,1); hold all; % shelf water flux time series
                % ax1(2) = subplot(2,1,2); hold all;

    hfig2 = []; %figure; hold all % vertical structure / baroclinicity

    hfig3 = []; %figure; hold all; % envelope

    hfig4 = []; %figure; subplot(2,1,1); hold all; % eddy water flux
                % subplot(2,1,2); hold all

    hfig5 = []; %figure; hold on;% along shelf structure

    hfig6 = figure; hold on; % streamer velocity

    hfig7 = []; %figure; hold on; % cross-shore bins

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    corder_backup = runArray.sorted_colors;

    nsmooth = 1;

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

        run = runArray.array(ii);
        names{ff} = getname(runArray, ii);

        if run.params.flags.flat_bottom
            continue;
        end

        if isempty(run.csflux) || ...
            isobath > length(run.csflux.x)
            disp(['Skipping ' run.name]);
            continue;
        end

        locstr = num2str(run.csflux.ndloc(isobath), 2);

        tind = find_approx(run.ndtime, 1.5, 1);
        Lz = run.eddy.Lgauss(1);
        hsb = run.bathy.hsb;

        He = hsb;
        Le = run.eddy.vor.dia(1)/2;

        fluxscl = 1; %run.eddyfluxscale;
        transscl = 1; He * Le^2;

        fluxvec = run.csflux.west.slope(:,isobath, source);
        ifluxvec = run.csflux.west.itrans.slope(:,isobath, source);
        ttrans = max(abs(ifluxvec));

        ndtime = run.csflux.time/run.eddy.turnover;

        R = run.csflux.R;
        [~,~,restind] = run.locate_resistance;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHELF WATER
        if ~isempty(hfig1)
            figure(hfig1)
            axes(ax1(1));
            [maxf, maxi] = run.calc_maxflux(fluxvec);
            hgplt1(ff) = plot(ndtime, ...
                              smooth(fluxvec/fluxscl, nsmooth));
            plot(ndtime(maxi), maxf/fluxscl, 'x', ...
                 'Color', hgplt1(ff).Color);
            axes(ax1(2));
            plot(ndtime, ifluxvec/transscl, 'Color', hgplt1(ff).Color);
        else
            colors = get(gca, 'ColorOrder');
            hgplt1(ff).Color = colors(ff,:);
        end

        %%%%%% BAROCLINICITY
        if ~isempty(hfig2)
            figure(hfig2)
            profile = ...
                run.csflux.west.slopewater.vertitrans(:,isobath,source)./ ttrans;
            vertbins = run.csflux.vertbins(:,isobath);
            zvec = vertbins./ max(abs(vertbins));
            zind = find_approx(vertbins, -1*run.bathy.hsb);
            bc = baroclinicity(zvec, profile);
            hgplt2(ff) = plot(profile, zvec, 'Color', hgplt1(ff).Color);
            plot(profile(zind), zvec(zind), 'x', ...
                 'Color', hgplt1(ff).Color);
            names2{ff} =  [names{ff} ' | bc = ' num2str(bc,'%.3f')];
        end

        %%%%%%%%%%%%% ENVELOPE
        if ~isempty(hfig3)
            figure(hfig3)
            % change from envelope to depth
            env = run.csflux.west.slopewater.envelope(:,isobath,1);
            env(isnan(env)) = max(env);
            ind = vecfind(run.rgrid.y_rho(:,1), env);
            %metric = run.bathy.h(1,ind)./run.bathy.hsb .* ...
            %         (1+run.rgrid.f(run.bathy.isb,1)./run.rgrid.f(ind,1))';

            hgplt3(ff) = plot(ndtime, (run.csflux.x(isobath) - env), ...
                              'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDDY WATER
        if ~isempty(hfig4)
            figure(hfig4)
            subplot(2,1,1)
            hgplt4(ff) = plot(ndtime, ...
                              smooth((run.csflux.east.eddy(:, isobath))/fluxscl, ...
                                     nsmooth), 'Color', hgplt1(ff).Color);

            subplot(2,1,2)
            plot(ndtime, ...
                 run.csflux.east.itrans.eddy(:,isobath)/transscl, ...
                 'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%% ALONG SHELF STRUCTURE
        if ~isempty(hfig5)
            if ~isfield(run.csflux, 'shelfx')
                run.streamerstruct;
            end

            figure(hfig5)
            R = run.csflux.R;
            L = run.eddy.vor.dia(1)/2;
            Lz0 = Lz(1);
            xi = run.csflux.slopex.xi;

            yoR = run.csflux.ndloc(isobath); % y/R - used in csflux
            y0oL = R/L * (1 - yoR); % y0/L - used in derivation
            xfrac = sqrt(1 - y0oL^2);

            xloc = find_approx(xi./L, -xfrac, 1);

            xi = run.csflux.slopex.xi;
            slopex = run.csflux.slopex.flux(:,isobath, source);
            hgplt5(ff) = plot(xi./L, slopex./ttrans, 'Color', ...
                              hgplt1(ff).Color);
            plot(xi(xloc)./L, slopex(xloc)./ttrans, 'x', ...
                 'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%% STREAMER VELOCITY
        if ~isempty(hfig6)
            figure(hfig6);

            yoR = run.csflux.ndloc; % y/R - used in csflux

            clear vnd y0oL
            for kk=1:length(yoR)
                [start, stop] = ...
                    run.flux_tindices(run.csflux.west.slopewater.vmax(:,kk));
                [vsmax,maxind] = ...
                    max(run.csflux.west.slopewater.vmax(start:stop,kk), [],1);
                % swirl velocity at ζ = 0 contour
                vscale = sqrt(2/exp(1))*run.eddy.V(start+maxind-1);
                vnd(kk) = vsmax/vscale;

                L = run.eddy.vor.dia(1)/2;
                y0oL(kk) =  R/L * (yoR(kk)-1); % y0/L - used in derivation
            end
            color = [1 1 1]*0.75;hgplt1(ff).Color;
            hgplt6(ff) = plot(y0oL, vnd, 'Color', color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%% CROSS-SHORE BINS
        if ~isempty(hfig7)
            figure(hfig7)
            xsb = run.bathy.xsb;

            bins = avg1(run.csflux.west.bins{isobath}, 2);
            itrans = run.csflux.west.slopewater.itrans{isobath};
            hgplt7(ff) = plot((bins-xsb)./R, itrans./ttrans);
        end
    end

    pause(1); drawnow;
    if ~isempty(hfig1)
        figure(hfig1)
        insertAnnotation('runArray.plot_fluxes');
        subplot(2,1,1)
        title(['Slope water | ND isobath = ', locstr]);
        ylabel('Flux / Eddy Flux above sb');
        liney(0, [], [1 1 1]*0.75);
        beautify;

        subplot(2,1,2)
        ylabel('Vol_{transport}/ Vol_{eddy}');
        xlabel('Non-dimensional time');
        legend(hgplt1, names);
        liney(0);
        beautify;
    end

    if ~isempty(hfig2)
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
        insertAnnotation('runArray.plot_fluxes');
        ylim([-1 0]);
        limx = xlim;
        xlim([0 limx(2)]);
        xlabel('Normalized volume transported');
        ylabel('Vertical bin / Water depth');
        title(['ND isobath = ', locstr]);
        legend(hgplt2, names2, 'Location', 'SouthEast');
        beautify;
    end

    if ~isempty(hfig3)
        figure(hfig3)
        insertAnnotation('runArray.plot_fluxes');
        xlabel('Non-dimensional time');
        ylabel('Distance from shelfbreak');
        legend(hgplt3, names);
        beautify;
    end

    if ~isempty(hfig4)
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

    pause(0.5);
    if ~isempty(hfig5)
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

    if ~isempty(hfig6)
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

    pause(0.5);
    if ~isempty(hfig7)
        figure(hfig7)
        insertAnnotation('runArray.plot_fluxes');
        linex(0, 'shelfbreak');
        linex(str2double(locstr), locstr);
        legend(hgplt7, names, 'Location', 'NorthWest');
        ylabel('Transport');
        xlabel('Source (y/R)');
        beautify;
    end

    runArray.reset_colors(corder_backup);
end