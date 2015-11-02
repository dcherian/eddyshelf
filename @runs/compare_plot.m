% comparison plots
function [] = compare_plot(runs,num)
    eddy = runs.eddy;
    % 86400 since eddy.t is already in days
    eddy.t = eddy.t./ (eddy.tscale/86400);
    ii = num;

    % line styles, markers & colors
    colors = cbrewer('qual', 'Dark2', 8);
    linestyle = {'-','--','-.','-'};
    markers = {'none','none','none','.'};

    aa = 6; bb = aa*2;
    tloc = [1:0.5:floor(max(eddy.t))];
    tind = vecfind(eddy.t, tloc);

    % plot eddy tracks
    % background velocity displacement
    if ~isfield(runs.eddy, 'bgvel')
        runs.eddy_bgflow();
    end
    displace = cumtrapz(runs.time, runs.eddy.bgvel);
    plotx = (eddy.mx - displace - eddy.mx(1))/1000;
    ploty = (eddy.my - eddy.my(1))/1000;
    figure(1)
    hold on
    %subplot(aa,2,bb);
    %subplot(aa,2,1); hold all
    %pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);            colorbar
    xlabel('X (km)'); ylabel('Y (km)');
    he = plot(plotx, ploty, 'Color',colors(ii,:),'LineWidth',2);
    hold all;
    addlegend(he, runs.name, 'NorthEast');
    try
        plot(plotx(tind), ploty(tind),'*', ...
             'MarkerSize',12,'Color',colors(ii,:),'LineWidth',2);
    catch ME
    end
    %if runs.bathy.axis == 'x'
    %    plot(eddy.we/1000,eddy.cy/1000,'Color',colors(ii,:),'LineStyle','--');
    %else
    %    plot(eddy.cx/1000,eddy.se/1000,'Color',colors(ii,:),'LineStyle','--');
    %end
    %axis image; axis tight

    % plot eddy properties
    figure(2)
    hold on
    subplot(aa,2,2); hold on
    limx = [0 max([xlim eddy.t])];
    plot(eddy.t,eddy.vor.amp./eddy.amp(1),'Color',colors(ii,:));
    ylabel('amp/amp(t=0) ');xlim(limx);

    if isfield(eddy, 'vol')
        subplot(aa,2,1); hold on
        plot(eddy.t, eddy.vol./eddy.vol(1),'Color', colors(ii,:));
        ylabel('Volume');xlim(limx);

        subplot(aa,2,3); hold on
        plot(eddy.t, eddy.PV./abs(eddy.PV(1)), 'Color', colors(ii,:));
        ylabel('PV/|PV0|');xlim(limx);

        subplot(aa,2,5); hold on
        plot(eddy.t, eddy.RV./abs(eddy.RV(1)), 'Color', colors(ii,:));
        ylabel('RV/|RV0|');xlim(limx);

        subplot(aa,2,7); hold on;
        plot(eddy.t, eddy.KE./eddy.KE(1), 'Color', colors(ii,:));
        ylabel('KE');xlim(limx);

        subplot(aa,2,9); hold on;
        plot(eddy.t, eddy.PE./eddy.PE(1), 'Color', colors(ii,:));
        ylabel('PE');xlim(limx);
    end

    subplot(aa,2,4); hold on
    plot(eddy.t, eddy.Ls/runs.rrdeep,'Color',colors(ii,:));
    ylabel('Ls/RRdeep');xlim(limx);

    subplot(aa,2,6); hold on
    plot(eddy.t,eddy.cvx,'Color',colors(ii,:));
    ylabel('cvx(km/day)');
    ylim([-5 5]);
    liney(0); xlim(limx);
    %plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
    %ylabel('x - center (km)');

    subplot(aa,2,8); hold on
    plot(eddy.t,eddy.cvy,'Color',colors(ii,:));
    ylabel('cvy (km/day)');xlim(limx);
    ylim([-5 5]);
    %plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
    %ylabel('y - center (km)');

    subplot(aa,2,10); hold on
    plot(eddy.t,eddy.Lgauss./max(eddy.Lgauss(1)),'Color',colors(ii,:));
    ylabel('H_{eddy}/H_{eddy0}');xlim(limx);
    %xlabel('time (days)');

    subplot(aa,2,12); hold on
    plot(eddy.t,eddy.prox/1000,'Color',colors(ii,:));
    xlabel('time (days)');
    ylabel('Proximity (km)');xlim(limx);

    subplot(aa,2,11); hold on
    hp = plot(eddy.t,eddy.hcen./max(runs.bathy.h(:)),'Color',colors(ii,:));
    addlegend(hp,runs.name,'SouthWest');
    %        plot(eddy.t,runs.params.phys.f0 / sqrt(runs.params.phys.N2) ...
    %             *  runs.eddy.dia,'Color',colors(ii,:),'LineStyle','--');
    %        legend('H_{center}','f/N*dia');
    xlabel('Time / Time at which center reaches slopebreak');
    ylabel('H_{center}(m)/H_{max}');
    xlim(limx);

    %% plot fluxes
    %{
    if isfield(runs.csflux,'west')
        ftime = runs.csflux.time/eddy.tscale;
        figure(4);
        subplot(4,1,1);
        hold on;
        plot(ftime, runs.csflux.off.shelf(:,1), 'Color', colors(ii,:));
        ylabel('Shelf water flux - sb');
        title('West');
        xlim(limx);

        subplot(4,1,2);
        hold on;
        plot(ftime, runs.csflux.off.slope(:,1), 'Color', colors(ii,:));
        ylabel('Slope water flux - sb');
        xlim(limx);

        try
            subplot(4,1,3);
            hold on;
            plot(ftime, runs.csflux.off.pv(:,1), 'Color', colors(ii,:));
            ylabel('PV flux');
            xlim(limx);
            catch ME
                end

                try
                    subplot(4,1,4);
                    hold on;
                    plot(ftime, runs.csflux.off.rv(:,1), 'Color', colors(ii,:));
                    ylabel('RV flux');
                    xlim(limx);
                    catch ME
                        end

                        figure(5);
                        subplot(4,1,1);
                        hold on;
                        plot(ftime, runs.csflux.east.shelf(:,1), 'Color', colors(ii,:));
                        ylabel('Shelf water flux - sb');
                        title('East');
                        xlim(limx);

                        subplot(4,1,2);
                        hold on;
                        plot(ftime, runs.csflux.east.slope(:,1), 'Color', colors(ii,:));
                        ylabel('Slope water flux - sb');
                        xlim(limx);

                        try
                            subplot(4,1,3);
                            hold on;
                            plot(ftime, runs.csflux.east.pv(:,1), 'Color', colors(ii,:));
                            ylabel('PV flux');
                            xlim(limx);

                            subplot(4,1,4);
                            hold on;
                            plot(ftime, runs.csflux.east.rv(:,1), 'Color', colors(ii,:));
                            ylabel('RV flux');
                            xlim(limx);
                            catch ME
                                end
                                end
                                %}
                                time = eddy.t;

                                %% plot water masses
                                if isfield(runs.water, 'off')
                                    % normalize volumes by initial eddy volume
                                    if isfield(runs.eddy, 'vol')
                                        evol0 = 1;runs.eddy.vol(runs.eddy.tscaleind);
                                    else
                                        evol0 = 1;
                                    end
                                    figure(3);
                                    set(gcf, 'Renderer', 'painters');
                                    % by regions
                                    % colors: off = r, sl = g, sh = b , edd = k, mix = m
                                    subplot(3,1,1)
                                    hold on;
                                    hw = plot(time, runs.water.sl.deep/evol0, 'Color', ...
                                              colors(ii,:), 'LineStyle', linestyle{1}, 'Marker', ...
                                              markers{1});
                                    plot(time, runs.water.sh.deep/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{2}, 'Marker', ...
                                         markers{2});
                                    ylabel('Deep region ');
                                    title('All volumes normalized by eddy volume at t=tscale');
                                    xlim(limx);
                                    addlegend(hw, runs.name, 'NorthWest');

                                    subplot(312)
                                    hold on;
                                    plot(time, runs.water.off.slope/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{1}, 'Marker', ...
                                         markers{1});
                                    plot(time, runs.water.sh.slope/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{3}, 'Marker', ...
                                         markers{3});
                                    plot(time, runs.water.edd.slope/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{4}, 'Marker', ...
                                         markers{4});
                                    %plot(time, runs.water.mix.slope/evol0, ['m' linestyle{num}]);
                                    legend('Offshore','Shelf','Eddy');
                                    ylabel('Slope region');
                                    xlim(limx);

                                    subplot(313)
                                    hold on;
                                    plot(time, runs.water.off.shelf/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{1}, 'Marker', ...
                                         markers{1});
                                    plot(time, runs.water.sl.shelf/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{2}, 'Marker', ...
                                         markers{2});
                                    plot(time, runs.water.edd.shelf/evol0, 'Color', ...
                                         colors(ii,:), 'LineStyle', linestyle{4}, 'Marker', ...
                                         markers{4});
                                    %plot(time, runs.water.mix.shelf/evol0,['m' linestyle{num}]);
                                    legend('Offshore','Slope','Eddy');
                                    ylabel('Shelf region');
                                    xlabel('Time (days)');
                                    xlim(limx);
                                end

                                % background flow velocity estimates
                                %{
                                figure(6)
                                subplot(211); hold on;
                                hbg = plot(time, runs.eddy.bgvel, 'Color', colors(ii,:));
                                ylabel('mean(vel. at x=eddy center)');
                                addlegend(hbg, runs.name, 'NorthWest');
                                subplot(212); hold on;
                                if runs.bathy.axis == 'y'
                                    plot(time, squeeze(mean(runs.ubar(3,:,:),2)), 'Color', ...
                                    colors(ii,:));
                                    else
                                        plot(time, squeeze(mean(runs.vbar(:,3,:),1)), 'Color', ...
                                        colors(ii,:));
                                        end
                                        xlabel('Time');
                                        ylabel('mean(inflow 2d vel)');
                                        %}

                                        % jet diagnostics
                                        if isfield('jet', runs)
                                            jtime = runs.time/86400;
                                            figure(7);
                                            subplot(3,2,1)
                                            plot(jtime, runs.jet.vscale, 'Color', colors(ii,:));
                                            ylabel('Max. velocity');
                                            subplot(3,2,2)
                                            plot(jtime, runs.jet.zscale, 'Color', colors(ii,:));
                                            hold on
                                            plot(jtime, -1*runs.jet.h, 'Color', colors(ii,:), 'LineStyle', ...
                                                 '-.');
                                            ylabel('z-loc of max vel');
                                            subplot(323)
                                            plot(jtime, runs.jet.bc, 'Color', colors(ii,:));
                                            liney(0);
                                            ylabel('Baroclinicity measure');
                                            subplot(324)
                                            plot(time, runs.jet.yscale, 'Color', colors(ii,:));
                                            ylabel('y-loc of max. vel');
                                        end

                                        figure(8);
                                        set(gcf, 'Renderer', 'painters');
                                        subplot(2,1,1)
                                        hold on
                                        plot(runs.csflux.time/eddy.tscale, ...
                                             runs.csflux.off.shelf(:,1)/1e6, 'Color', colors(ii,:));
                                        plot(runs.csflux.time/eddy.tscale, ...
                                             runs.csflux.east.slope(:,1)/1e6, 'Color', colors(ii,:), ...
                                             'LineStyle', linestyle{2});
                                        liney(0);
                                        ylabel('Transport (Sv)');
                                        legend('West - shelf water', 'East - slope water');

                                        subplot(2,1,2)
                                        hold on
                                        hline = plot(eddy.t,eddy.hcen, 'Color', colors(ii,:));
                                        xlabel('Time (days)');
                                        ylabel('center-isobath');
                                        %linex(tind);
                                        suplabel(runs.dir,'t');
                                        %        packrows(2,1);
                                        addlegend(hline, runs.name);

                                        % shelf water envelope
                                        if isfield(runs.csflux.off.shelfwater, 'envelope')
                                            normtrans = sum(runs.csflux.off.shelfwater.itrans);

                                            figure(9)
                                            subplot(2,1,1)
                                            hold on
                                            hline = plot(runs.csflux.time/86400, ...
                                                         runs.csflux.off.shelfwater.envelope/ runs.rrshelf, ...
                                                         'Color', colors(ii,:));
                                            addlegend(hline, runs.name);
                                            xlabel('Time');
                                            ylabel({'Location of water parcel farthest from shelfbreak' ...
                                                    'in terms of shelfbreak rossby radius'})

                                            subplot(2,1,2)
                                            hold on;
                                            plot(runs.csflux.off.shelfwater.bins/runs.rrshelf, ...
                                                 runs.csflux.off.shelfwater.itrans./normtrans, 'color', colors(ii,:));
                                            ylabel('Total volume transported');
                                            xlabel('Bin = location / RR_{shelf} ');
                                        end

                                        % shelf water vorticity budget
                                        if isfield(runs, 'vorbudget')
                                            figure(10)
                                            subplot(2,1,1)
                                            hold all
                                            hline = plot(runs.csflux.time/86400, runs.csflux.off.shelf, ...
                                                         'Color', colors(ii,:));
                                            addlegend(hline, runs.name);

                                            subplot(2,1,2)
                                            hold all
                                            hline = plot(runs.vorbudget.time/86400, ...
                                                         runs.vorbudget.shelf.str, 'Color', ...
                                                         colors(ii,:));
                                        end
                                    end
