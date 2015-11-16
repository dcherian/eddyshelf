function [] = plot_fluxes(runs, isobath, source)

    n = length(runs.csflux.x);
    ti = runs.csflux.time/86400;

    % source of water
    if ~exist('source', 'var'), source = 1; end
    % across which isobath
    if ~exist('isobath', 'var') | isempty(isobath), isobath = 1:n; end

    n = length(isobath);

    if any(isobath > length(runs.csflux.x))
        error(['Choose a lower isobath number (max = ' ...
               num2str(length(runs.csflux.x)) ')']);
    end

    [~,~,restind] = runs.locate_resistance;

    hfig1 = []; %figure; % streamer vertical structure
    hfig2 = []; %figure; % hovmoeller plots (x,t)
    hfig3 = []; %figure; % hovmoeller plots (z,t)
    hfig4 = []; %figure; % flux time-series for given source
    hfig5 = figure; % flux time-series at isobath
    hfig7 = []; %figure; % streamer vmax
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
