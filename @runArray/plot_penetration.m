% diagnostic plots to show how much eddy penetrates slope
function [] = plot_penetration(runArray, ax)

    mark_timestamp = 1;

    corder_backup = runArray.sorted_colors;

    hfig1 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2);
                %hold all; insertAnnotation('runArray.plot_penetration');

    if ~exist('ax', 'var')
        hfig2 = figure;ax1 = gca; hold all; %subplot(2,2,[1 3]); hold all;
    else
        hfig2 = gcf;
        ax1 = ax; hold all;
    end
    figure(hfig2); insertAnnotation('runArray.plot_penetration');
    ax2 = []; %subplot(2,2,2); hold all;
    ax3 = []; %subplot(2,2,4); hold all;

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    len = length(runArray.filter);
    colors = cbrewer('qual','Dark2',4);

    for ff=1:len
        ii = runArray.filter(ff);

        color = colors(ff,:);
        run  = runArray.array(ii);
        name = getname(runArray, ii);
        %name = [getname(runArray, ii) ' | L_{sl} = ' ...
        %        num2str(run.bathy.L_slope/1000) ' km'];
        %name = [getname(runArray, ii) ' | S_\alpha = ' ...
        %        num2str(run.bathy.S_sl)];
        %name = [' H_{sb}/H_{eddy} = ' ...
        %        num2str(run.bathy.hsb./run.eddy.Lgauss(1))];

        ndtime = run.eddy.t * 86400 / run.eddy.turnover;
        if mark_timestamp
            dt = 75;
            tinds = vecfind(run.eddy.t, [dt:dt:max(run.eddy.t)]);
        else
            tinds = [];
        end

        % normalization for axes
        xnorm = run.eddy.vor.dia(1)/2;
        ynorm = run.eddy.vor.dia(1)/2; run.eddy.my(1) - run.bathy.xsb;

        if run.bathy.axis == 'y'
            mx = run.eddy.mx;
            my = run.eddy.my;

            edge1 = run.eddy.vor.ne;
            edge2 = run.eddy.vor.se;

            % reference location
            x0 = mx(1);
            y0 = run.bathy.xsb;
            legloc = 'NorthWest';
            x0str = 'X_0';
            y0str = 'Y_{sb}';
        else
            mx = run.eddy.mx;
            my = run.eddy.my;

            edge1 = run.eddy.vor.ee;
            edge2 = run.eddy.vor.we;

            % reference location
            y0 = my(1);
            x0 = run.bathy.xsb;
            legloc = 'NorthEast';
            x0str = 'X_{sb}';
            y0str = 'Y_0';
        end

        if ~isempty(hfig1)
            figure(hfig1);
            subplot(2,1,1);
            hgplt1(ff) = plot(ndtime, (edge1 - y0)./ynorm, 'Color', ...
                              color);

            subplot(2,1,2);
            plot(ndtime, (edge2 - y0)./ynorm, 'Color', color);
        end

        axes(ax1)
        x = (mx - x0)/xnorm;
        y = (my - y0)/ynorm;

        % plot track
        hgplt2(ff) = plot(x, y, 'Color', color, 'LineStyle', '-');
        names{ff} = name;

        % mark start of resistance
        %[~,~,tind] = run.locate_resistance;
        %plot(x(tind), y(tind), '.', ...
        %     'Color', color, 'Markersize', 22);
        if run.bathy.L_slope/run.eddy.vor.dia(1) > 1
            run.fit_traj(1.0);
            plot(x(run.traj.tind), y(run.traj.tind), 'x', ...
                 'Color', color, 'Markersize', 16);
        end

        % mark timestamps
        if mark_timestamp
            plot(x(tinds), y(tinds), '.', 'Color', color, 'MarkerSize', ...
                 22);
            %text(x(tinds), y(tinds), ...
            %     cellstr(num2str(ndtime(tinds)', 2)));
        end

        if ~isempty(ax2)
            axes(ax2)
            hgplt = plot(ndtime, run.eddy.Lgauss./run.eddy.hcen', ...
                         'Color', color);
        end
        if ~isempty(ax3)
            axes(ax3)
            hgplt = plot(ndtime, run.eddy.Lgauss, 'Color', color);
        end
    end

    fontSize = [18 18 20];
    if ~isempty(hfig1)
        figure(hfig1);
        subplot(2,1,1)
        ylabel('Northern edge - Y_{sb}');
        legend(hgplt1, names);
        beautify(fontSize);
        subplot(2,1,2)
        ylabel('Southern edge - Y_{sb}');
        liney(0);
        beautify(fontSize);
    end

    axes(ax1)
    axis image;
    if ynorm == run.rrdeep
        ylabel(['(Y - ' y0str ')/(deformation radius)']);
        dytick = 1;
    else
        if ynorm == run.eddy.vor.dia(1)/2
            ylabel(['(Y - ' y0str ')/(initial radius)']);
            dytick = 1;
        else
            ylabel(['Distance from shelfbreak / Initial distance from ' ...
                    'shelfbreak']);
            dytick = 0.1;
        end
    end
    if run.bathy.axis == 'y'
        xlim([-12 2]);
        ylim([0 max(ylim)]);
        liney(1);
        set(gca, 'YTick', [0:dytick:max(ylim)]);
        title('Southern Coast');
    else
        xlim([0 max(xlim)]);
        ylim([min(ylim) -1*min(ylim)/2]);
        linex(1);
        title('Western Coast');
    end

    if xnorm == run.rrdeep
        xlabel(['(X - ' x0str ')/(deformation radius)']);
    else
        if xnorm == run.eddy.vor.dia(1)/2
            xlabel(['(X - ' x0str ')/(initial radius)']);
        end
    end

    %axis image;
    %    ylim([-2 max(ylim)]);

    hlegend = legend(hgplt2, names, 'Location', legloc, 'Box', ...
                     'off');
    beautify(fontSize);

    if ~isempty(ax2)
        axes(ax2)
        ylabel('Eddy vertical scale / water depth');
        ylim([0 0.6]);
        xlabel('Non-dimensional time');
        beautify;
    end
    if ~isempty(ax3)
        axes(ax3);
        ylabel('Eddy vertical scale (m)');
        beautify;
    end

    runArray.reset_colors(corder_backup);
end
