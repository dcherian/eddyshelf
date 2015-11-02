function [] = plot_fluxes(runs, isobath, source)

    n = length(runs.csflux.x);
    ti = runs.csflux.time/86400;

    % source of water
    if ~exist('source', 'var'), source = 1; end
    % across which isobath
    if ~exist('isobath', 'var'), isobath = 1:n; end

    n = length(isobath);

    if any(isobath > length(runs.csflux.x))
        error(['Choose a lower isobath number (max = ' ...
               num2str(length(runs.csflux.x)) ')']);
    end


    [~,~,restind] = runs.locate_resistance;

    hfig1 = []; %figure; % streamer vertical structure
    hfig2 = []; %figure; % hovmoeller plots (x,t)
    hfig3 = []; %figure; % hovmoeller plots (z,t)
    hfig4 = []; %figure; % flux time-series
    hfig5 = []; %figure; % flux time-series at isobath
    hfig6 = figure; % flux cross-sections
    hfig7 = []; %figure; % streamer vmax
    hfig8 = []; %figure; % debug flux parameterization
    hfig9 = []; %figure; % vertical structure at maxloc

    hsb = runs.bathy.hsb;
    Lz = runs.eddy.Lgauss(1);
    Ro = runs.eddy.Ro(1);

    if ~isempty(hfig1)
        if length(source) > 1
            error('Specify 1 source isobath only!');
        end
        figure(hfig1); % streamer vertical structure
        insertAnnotation([runs.name '.plot_fluxes']);
        subplot(121)
        lightDarkLines(n);
        plot(runs.csflux.off.slopewater.vertitrans(:,isobath,source), ...
             runs.csflux.vertbins(:,isobath));
        liney(-1 * runs.bathy.hsb);
        legend(cellstr(num2str(runs.csflux.h(isobath)', '%.2f')), ...
               'Location', 'SouthEast');
        beautify;
        title([runs.name ' | Source = ' ...
               num2str(runs.csflux.h(source), '%d') ...
               'm |  index = ' num2str(source)]);
        axis square;

        subplot(122)
        lightDarkLines(length(runs.csflux.ix));
        hold on;
        hplt = nan([length(runs.csflux.ix) 1]);
        isochoose = 2:5;
        for iso=isochoose %length(runs.csflux.ix)
                          % this doesn't account for eddy transiting through isobath
                          %fluxvec = runs.csflux.off.slopewater.vertitrans(:, ...
                          %                                                 iso,iso);
                          % [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:, iso,iso));
            [~,~,start] = runs.locate_resistance;
            fluxvec = trapz(runs.csflux.time(start:end), ...
                            runs.csflux.off.slopezt(:,start:end,iso,iso), 2);

            zvec = runs.csflux.vertbins(:,iso)./hsb;
            trans = trapz(zvec, fluxvec);

            hplt(iso) = plot(fluxvec, zvec);
            [~,~,~,idiff] = runs.streamer_peak(iso);
            plot(fluxvec(idiff), zvec(idiff), 'kx');

            %ideal = runs.streamer_ideal_profile(iso);
            %plot(ideal*max(fluxvec), zvec, 'k--');
        end
        %ideal = exp(-((zvec + hsb/2)/(hsb)).^2) .* (zvec >= -hsb) ...
        %        + exp(-1/4) * (1 - erf(-zvec/Lz)) .* (zvec < -hsb);
        %plot(ideal, zvec, 'k--');
        %liney(-1 * runs.bathy.hsb);
        legend(cut_nan(hplt), cellstr(num2str(runs.csflux.h(isochoose)', '%.0f')), ...
               'Location', 'SouthEast');
        beautify;
        title([runs.name]);
        axis square;
    end

    if ~isempty(hfig2)
        runs.streamerstruct;

        figure(hfig2); % hovmoeller plot of fluxes (x,t)
        insertAnnotation([runs.name '.plot_fluxes']);
        for index=isobath
            slopexti = runs.csflux.slopex.slopexti(:,:,index,source);
            xi = runs.csflux.slopex.xi;

            ax(index) = subplot(2,ceil(n/2),index);
            pcolorcen(xi/1000, ti, slopexti');
            xlabel('X - X_{eddy} (km)'); ylabel('Time (days)');
            title(['Y = ' num2str(runs.csflux.x(index)/1000, '%.2f'), ...
                   ' km | H = ' num2str(runs.csflux.h(index), '%d') ...
                   ' m | ' runs.name]);
            hold on;
            plot((runs.eddy.vor.we-runs.eddy.mx)/1000, ti);
            plot((runs.eddy.vor.ee-runs.eddy.mx)/1000, ti);
            center_colorbar;
            linex(0);
        end
    end

    if ~isempty(hfig3)
        figure(hfig3); % hovmoeller plot of fluxes (z,t)
        insertAnnotation([runs.name '.plot_fluxes']);
        for index=isobath
            zvec = runs.csflux.vertbins(:,index);
            slopezt = runs.csflux.off.slopezt(:,:,index,source);

            subplot(2,ceil(n/2),index)
            pcolorcen(ti, zvec, slopezt);
            center_colorbar;
            liney(-hsb);
            title(['Y = ' num2str(runs.csflux.x(index)/1000, '%.2f'), ...
                   ' km | H = ' num2str(runs.csflux.h(index), '%d') ...
                   ' m | ' runs.name]);
        end
    end

    if ~isempty(hfig4)
        figure(hfig4)
        insertAnnotation([runs.name '.plot_fluxes']);
        lightDarkLines(n);
        plot(ti, runs.csflux.off.slope(:,isobath, source));
        legend(cellstr(num2str(runs.csflux.h', '%.0f')), ...
               'Location', 'NorthWest');
        title(runs.name);
        [~,~,tind] = runs.locate_resistance;
        linex(ti(tind));
        beautify;
    end

    if ~isempty(hfig5)
        figure(hfig5)
        insertAnnotation([runs.name '.plot_fluxes']);
        lightDarkLines(n);
        for ii=1:length(isobath)
            hplt(ii) = plot(ti, ...
                            runs.csflux.off.slope(:,isobath(ii), ...
                                                   isobath(ii)));
            [~, maxloc] = runs.calc_maxflux( ...
                runs.csflux.off.slope(:,isobath(ii), ...
                                       isobath(ii)));
            plot(ti(maxloc), runs.csflux.off.slope(maxloc,isobath(ii), ...
                                                    isobath(ii)), ...
                 'kx')
        end
        linex(ti(restind), 'resistance');
        legend(hplt, cellstr(num2str(runs.csflux.h', '%.0f')), ...
               'Location', 'NorthWest');
        title(runs.name);
        ylabel('Flux (m^3/s)');
        xlabel('Time');
        beautify;
    end

    if ~isempty(hfig6) || ~isempty(hfig8)
        if length(isobath) > 1
            close(hfig6);
            error('Specify a *single* isobath!');
        end

        ix = runs.csflux.ix(isobath);
        [~,tindex] = runs.calc_maxflux(...
            runs.csflux.off.slope(:,isobath,isobath));;

        % tindex = 55;
        % vnorm = runs.eddy.V(tindex);

        % copied from csfluxes.m
        if runs.csvelname == 'v'
            bathyax = 2;
            xvec = runs.rgrid.x_rho(1,2:end-1) - runs.eddy.mx(tindex);
            zvec = runs.rgrid.z_r(:, ix+1, 1);
        else
            bathyax = 1;
            xvec = (runs.rgrid.y_rho(2:end-1,1) - runs.eddy.my(tindex))';
            zvec = runs.rgrid.z_r(:, 1, ix+1);
        end
        csvel = squeeze(avg1( ...
            dc_roms_read_data(runs.dir, runs.csvelname, tindex, ...
                              {runs.bathy.axis ix-1 ix}, [], runs.rgrid, ...
                              'his', 'single'), bathyax));
        csvel = csvel(2:end-1,:,:,:);

        % process cross-shelf dye
        csdye = dc_roms_read_data(runs.dir, runs.csdname, ...
                                  tindex, {runs.bathy.axis ix+1 ix+1}, ...
                                  [], runs.rgrid, 'his', 'single');
        csdye = csdye(2:end-1,:,:);

        eddye = dc_roms_read_data(runs.dir, runs.eddname, ...
                                  tindex, {runs.bathy.axis ix+1 ix+1}, ...
                                  [], runs.rgrid, 'his', 'single');
        eddye = eddye(2:end-1,:,:);

        rho = dc_roms_read_data(runs.dir, 'rho', ...
                                tindex, {runs.bathy.axis ix+1 ix+1}, ...
                                [], runs.rgrid, 'his', 'single');
        rho = rho(2:end-1,:,:);

        rback = dc_roms_read_data(runs.dir, 'rho', ...
                                  1, {runs.bathy.axis ix+1 ix+1}, ...
                                  [], runs.rgrid, 'his', 'single');
        % rho = rho - rback(2:end-1,:,:);

        mask = fillnan(bsxfun(@times, csdye < runs.csflux.x(isobath), ...
                              runs.csflux.westmask(:,tindex,isobath)),0);

        % profile I am assuming
        [videal, idmask] = runs.makeStreamerSection(isobath);

        tind = tindex; % FOR PARAMETERISATION
                       % syms x z;
        a = 2; % 2 for gaussian
               % V0 = runs.eddy.V(tind) * 2.33;
        R = runs.csflux.R;
        L = median(runs.eddy.rhovor.dia(1:tind))/2;
        % Lz = runs.eddy.Lgauss(tindex);
        % H = runs.csflux.h(isobath);
        yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
        y0oL =  R/L * (1 - yoR); % y0/L - used in derivation
        xfrac = -sqrt(1 - y0oL^a);

        % xvec = xvec ./ L;

        % %xfrac2 = -sqrt(1 - y0oL^2 - (zvec./Lz).^2);
        % %xfrac2(~isreal(xfrac2)) = 0;
        % % eddy profile
        % eddmask = bsxfun(@le, abs(xvec), sqrt(1-y0oL^a-(zvec./Lz).^2));
        % % mask to integrate over
        % inmask = bsxfun(@and, xvec < 0, 1 - eddmask);

        %if runs.bathy.hsb/Lz < 0.5
        %    inmask = repmat(inmask(end,:), [runs.rgrid.N 1]);
        %end

        % videal = bsxfun(@times, ...
        %                 -V0/vnorm * (xvec).^(a-1) .* exp(-xvec.^a -y0oL^a), ...
        %                 (1 - erf(-zvec/Lz)));
        % % idmask = repmat(xvec < xfrac, [runs.rgrid.N 1]);

        % diagnosed flux
        flux = runs.csflux.off.slope(tindex,isobath,isobath)
        % should agree with above
        vtrue = trapz(xvec, trapz(zvec, repnan(csvel .* mask,0), 2), 1)

        % idealized velocity with real mask
        vitruemask = trapz(xvec, trapz(zvec, repnan(videal .* mask, 0), 2), 1)

        % idealized parameterization
        vest = trapz(xvec, trapz(zvec, repnan(videal .* idmask, 0), 2), 1)

        if ~isempty(hfig6)
            figure(hfig6)
            insertAnnotation([runs.name '.plot_fluxes']);
            maximize;

            ax(1) = subplot(221);
            [~,h1] = contourf(xvec/1000, zvec, csvel', 20);
            hold on
            contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
            runs.add_timelabel(tindex);
            xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
            title(['Cross-shelf velocity (m/s) | ' runs.name]);
            %linex(xfrac);
            liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
            liney(-1 * runs.bathy.hsb, 'h_{sb}');
            caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
            beautify;

            figure(hfig6)
            ax(2) = subplot(224);
            pcolorcen(xvec/1000, zvec, runs.sgntamp*(csdye'-runs.bathy.xsb)/1000); % .* mask
            hold on; shading interp
            contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
            hcb = colorbar;
            xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
            title(['Cross-shelf dye (km) | ' runs.name]);
            runs.add_timelabel(tindex);
            linkaxes(ax, 'xy');
            %linex(xfrac);
            liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
            liney(-1 * runs.bathy.hsb, 'h_{sb}');
            hcb.TickLabels{1} = 'Shelf Water';
            hcb.TickLabels{end-1} = 'Eddy Water';
            hcb.TickLabels{floor(length(hcb.Ticks)/2)-2} = 'Slope Water';
            beautify;

            figure(hfig6)
            ax(3) = subplot(223);
            contour(xvec/1000, zvec, rho', 30, 'k'); % .* mask
            clim = caxis;
            hold on; shading interp;
            contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
            caxis(clim);
            colorbar;
            xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
            title(['\rho (kg/m^3) | ' runs.name]);
            runs.add_timelabel(tindex);
            linkaxes(ax, 'xy');
            %linex(xfrac);
            liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
            liney(-1 * runs.bathy.hsb, 'h_{sb}');
            beautify;

            figure(hfig6)
            % ax(4) = subplot(224);
            % pcolorcen(xvec/1000, zvec, eddye'); % .* mask
            % hold on; shading interp;
            % contour(xvec/1000, zvec, repnan(mask,0), [1 1], 'k', 'LineWidth', 2);
            % colorbar;
            % xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
            % title(['Eddy dye (km) | ' runs.name]);
            % runs.add_timelabel(tindex);
            % linkaxes(ax, 'xy');
            % linex(xfrac);
            % liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
            % liney(-1 * runs.bathy.hsb, 'h_{sb}');
            % caxis([0 1]);

            figure(hfig6)
            subplot(222);
            plot(abs(trapz(xvec, repnan(csvel.*mask,0), 1)), zvec);
            beautify;
            ylabel('Z (m)');
            xlabel('|Transport| (m^2/s)');
            liney(-runs.bathy.hsb, 'shelfbreak depth');
            liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
            ylim(ax(1).YLim);
            beautify;
        end

        if ~isempty(hfig8)
            figure(hfig8);
            insertAnnotation([runs.name '.plot_fluxes']);
            maximize;

            ax(1) = subplot(2,3,[1 2]);
            [~,h1] = contourf(xvec/1000, zvec, csvel', 20);
            hold on
            contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
            runs.add_timelabel(tindex);
            xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
            title(['Cross-shelf velocity (m/s) | ' runs.name]);
            linex(xfrac*L/1000, 'xfrac');
            liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
                  {'vertical scale'; 'h_{sb}'});
            caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
            beautify;

            ax(2) = subplot(2,3,[4 5]);
            [~,h1] = contourf(xvec/1000, zvec, videal', 20);
            hold on
            contour(xvec/1000, zvec, repnan(idmask',0), [1 1], 'k', 'LineWidth', 2);
            runs.add_timelabel(tindex);
            xlabel('(X - X_{eddy})/L_{eddy}'); ylabel('Z (m)');
            title('Idealized');
            linex(xfrac*L/1000, 'xfrac');
            liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
                  {'vertical scale'; 'h_{sb}'});
            caxis([-1 1] * max(abs(csvel(:)))); center_colorbar;
            beautify;

            subplot(2,3,[3 6]);
            plot(trapz(xvec, repnan(csvel.*mask,0), 1), zvec);
            hold on;
            plot(trapz(xvec, videal.*idmask, 1), zvec);
            legend('Actual', 'Ideal', 'Location', 'SouthEast');
            beautify;

            linkaxes(ax, 'xy');
        end
    end

    if ~isempty(hfig7)
        figure(hfig7)
        insertAnnotation([runs.name '.plot_fluxes']);

        lightDarkLines(length(runs.csflux.ix));
        plot(runs.csflux.time/86400, ...
             runs.csflux.off.slopewater.vmax(:,isobath) ./ ...
             runs.eddy.V(1));
        liney(1);
        legend(cellstr(num2str(runs.csflux.ndloc(isobath)', '%.2f')), ...
               'Location', 'NorthWest');
        title(runs.name);
        ylabel('Max streamer velocity / eddy velocity scale');
        xlabel('Time (days)');
        beautify;
    end

    if ~isempty(hfig9)
        figure(hfig9)
        insertAnnotation([runs.name '.plot_fluxes']);

        subplot(121);
        lightDarkLines(length(isobath));
        for ii=isobath
            [~,maxloc] = runs.calc_maxflux(ii);
            fluxvec = runs.csflux.off.slopezt(:,maxloc,ii,ii);
            zvec = runs.csflux.vertbins(:,ii);
            ideal = runs.streamer_ideal_profile(ii);

            plot(fluxvec./max(fluxvec), zvec);
            plot(ideal, zvec, 'k--');
            if length(isobath) == 1
                fluxvec = runs.csflux.off.slopewater.vertitrans(:,ii,ii);
                plot(fluxvec./max(fluxvec), zvec);
                legend('Instantaneous', 'Ideal', 'Total', ...
                       'Location', 'SouthEast');
            end
        end
        title([runs.name ' | At t=maxflux, source=isobath']);
        ylabel('Z (m)');
        beautify;

        subplot(122);
        lightDarkLines(length(isobath));
        for ii=isobath
            [~,maxloc] = runs.calc_maxflux(ii);
            plot(runs.csflux.off.slopezt(:,maxloc,ii,source), ...
                 runs.csflux.vertbins(:,ii));
        end
        legend(num2str(runs.csflux.h(isobath)', '%d'), 'Location', 'SouthEast');
        title([runs.name ' | Source = ' num2str(runs.csflux.h(source), '%d') ...
               ' m']);
        ylabel('Z (m)');
        beautify;
    end
end
