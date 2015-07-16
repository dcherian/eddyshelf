function [] = plot_fluxes(runArray, index)

    if ~exist('index', 'var'), index = 2; end
    hfig1 = figure; ax1(1) = subplot(2,1,1); hold all; % shelf water flux time series
            ax1(2) = subplot(2,1,2); hold all;

    hfig2 = figure; hold all % vertical structure / baroclinicity
    % subplot(2,2,1); hold all
    % subplot(2,2,2); hold all
    % subplot(2,2,3); hold all
    % subplot(2,2,4); hold all

    hfig3 = figure; hold all; % envelope

    hfig4 = figure; subplot(2,1,1); hold all; % eddy water flux
    subplot(2,1,2); hold all

    hfig5 = figure; hold on;% along shelf structure

    hfig6 = figure; hold on; % streamer velocity

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

        if isempty(run.csflux)
            disp(['Skipping ' run.name]);
            continue;
        end

        tind = find_approx(run.ndtime, 1.5, 1);
        Lz = run.eddy.Lgauss(1);
        hsb = run.bathy.hsb;

        He = hsb;
        Le = run.eddy.vor.dia(1)/2;

        fluxscl = 1; %run.eddyfluxscale;
        transscl = 1; He * Le^2;

        fluxvec = run.csflux.west.slope(:,index);
        ifluxvec = run.csflux.west.itrans.slope(:,index);
        ttrans = max(abs(ifluxvec));

        ndtime = run.csflux.time/run.eddy.turnover;

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
                    run.csflux.west.slopewater.vertitransw(:,index)./ ttrans;
                vertbins = run.csflux.west.vertbins(:,index);
                zvec = vertbins./ max(abs(vertbins));
                bc = baroclinicity(zvec, profile);
                hgplt2(ff) = plot(profile, zvec, 'Color', hgplt1(ff).Color);
                names2{ff} =  [names{ff} ' | bc = ' num2str(bc,'%.3f')];
            catch ME
                keyboard;
                disp(ME)
            end
        end

        %%%%%%%%%%%%% ENVELOPE
        if ~isempty(hfig3)
            figure(hfig3)
            % change from envelope to depth
            env = run.csflux.west.slopewater.envelope(:,index);
            env(isnan(env)) = max(env);
            ind = vecfind(run.rgrid.y_rho(:,1), env);
            %metric = run.bathy.h(1,ind)./run.bathy.hsb .* ...
            %         (1+run.rgrid.f(run.bathy.isb,1)./run.rgrid.f(ind,1))';

            hgplt3(ff) = plot(ndtime, (run.csflux.x(index) - env), ...
                              'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDDY WATER
        if ~isempty(hfig4)
            figure(hfig4)
            subplot(2,1,1)
            hgplt4(ff) = plot(ndtime, ...
                              smooth((run.csflux.east.eddy(:, ...
                                                           index))/fluxscl, ...
                                     nsmooth), 'Color', hgplt1(ff).Color);

            subplot(2,1,2)
            plot(ndtime, ...
                 run.csflux.east.itrans.eddy(:,index)/transscl, ...
                 'Color', hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%% ALONG SHELF STRUCTURE
        if ~isempty(hfig5)
            figure(hfig5)
            Lx = run.eddy.vor.dia(1)/2;
            %if ~isfield(run.csflux, 'shelfx')
                run.streamerstruct(index);
                %end
            xi = run.csflux.slopex.xi;
            slopex = run.csflux.slopex.flux;
            hgplt5(ff) = plot(xi./Lx, slopex./ttrans, 'Color', ...
                              hgplt1(ff).Color);
        end

        %%%%%%%%%%%%%%%%%%%%%%%% STREAMER VELOCITY
        if ~isempty(hfig6)
            figure(hfig6)
            vxt = run.csflux.slopext(:,:,1) .* ...
                  run.csflux.westmask./run.bathy.hsb;

            hgplt6(ff) = plot(ndtime, max(vxt,[],1)./run.eddy.V(1), ...
                              'Color', hgplt1(ff).Color);
        end
    end

    pause(1); drawnow;
    locstr = num2str(run.csflux.ndloc(index), 2);
    if ~isempty(hfig1)
        figure(hfig1)
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
        ylim([-1 0]);
        limx = xlim;
        xlim([0 limx(2)]);
        xlabel('Normalized volume transported');
        ylabel('Vertical bin / Water depth');
        legend(hgplt2, names2, 'Location', 'SouthWest');
        beautify;
    end

    if ~isempty(hfig3)
        figure(hfig3)
        xlabel('Non-dimensional time');
        ylabel('Distance from shelfbreak');
        legend(hgplt3, names);
        beautify;
    end

    if ~isempty(hfig4)
        figure(hfig4)
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

    if ~isempty(hfig5)
        figure(hfig5)
        xlabel('(X - X_e)/L_x');
        ylabel('Normalized flux');
        title('\int v(x,z,t) dz dt');
        linex([-1 0 1]); liney(0);
        xlim([-2 2]);
        legend(hgplt5, names);
        beautify;
    end

    if ~isempty(hfig6)
        figure(hfig6)
        ylabel('Max shelf water velocity');
        liney(0);
        legend(hgplt6, names);
        beautify;
    end

    runArray.reset_colors(corder_backup);
end