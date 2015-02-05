% diagnostic plots to show how much eddy penetrates slope
function [] = plot_penetration(runArray)

    mark_timestamp = 0;

    corder_backup = runArray.sorted_colors;

    hfig1 = figure; subplot(2,1,1); hold all; subplot(2,1,2);
            hold all; insertAnnotation('runArray.plot_penetration');

    hfig2 = figure;
    insertAnnotation('runArray.plot_penetration');
    ax1 = gca; hold all; %subplot(2,2,[1 3]); hold all;
    ax2 = []; %subplot(2,2,2); hold all;
    ax3 = []; %subplot(2,2,4); hold all;

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

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
            tinds = vecfind(ndtime, [0.5:0.5:(max(ndtime))]);
        else
            tinds = [];
        end

        % normalization for axes
        xnorm = run.rrdeep; run.eddy.vor.dia(1)/2;
        ynorm = run.rrdeep; run.eddy.vor.dia(1)/2; run.eddy.my(1) - run.bathy.xsb;

        % reference location
        x0 = run.eddy.mx(1);
        %y0 = run.eddy.my(1);
        y0 = run.bathy.xsb;

        if ~isempty(hfig1)
            figure(hfig1);
            subplot(2,1,1);
            hgplt = plot(ndtime, (run.eddy.vor.ne - y0)./ynorm);
            addlegend(hgplt, name);

            subplot(2,1,2);
            plot(ndtime, (run.eddy.vor.se - y0)./ynorm);
        end

        figure(hfig2)
        axes(ax1)
        x = (run.eddy.mx - x0)/xnorm;
        y = (run.eddy.my - y0)/ynorm;

        % plot track
        hgplt = plot(x, y);
        addlegend(hgplt, name, 'NorthWest');

        % mark start of resistance
        if isfield(run.eddy, 'res')
            plot(x(run.eddy.res.tind), y(run.eddy.res.tind), '.', ...
                 'Color', get(hgplt, 'Color'), 'Markersize', 22);
        end
         
        % mark timestamps
        if mark_timestamp
            plot(x(tinds), y(tinds), 'kx');
            text(x(tinds), y(tinds), ...
                 cellstr(num2str(ndtime(tinds)', 2)));
        end

        if ~isempty(ax2)
            axes(ax2)
            hgplt = plot(ndtime, run.eddy.Lgauss./run.eddy.hcen');
        end
        if ~isempty(ax3)
            axes(ax3)
            hgplt = plot(ndtime, run.eddy.Lgauss);
        end
    end

    fontSize = [16 16 18];
    if ~isempty(hfig1)
        figure(hfig1);
        subplot(2,1,1)
        ylabel('Northern edge - Y_{sb}');
        beautify(fontSize);
        subplot(2,1,2)
        ylabel('Southern edge - Y_{sb}');
        liney(0);
        beautify(fontSize);
    end

    figure(hfig2)
    axes(ax1)
    if ynorm == run.rrdeep
        ylabel('(Y - Y_{sb})/(deformation radius)');
        liney(1);
        set(gca, 'YTick', [0:1:max(ylim)]);
    else
        if ynorm == run.eddy.vor.dia(1)/2
            ylabel('(Y - Y_{sb})/(initial radius)');
            liney(1);
            set(gca, 'YTick', [0:1:max(ylim)]);
        else
            ylabel(['Distance from shelfbreak / Initial distance from ' ...
                    'shelfbreak']);
            liney(1);
            set(gca, 'YTick', [0:0.1:max(ylim)]);
        end
    end

    if xnorm == run.rrdeep
        xlabel('(X-X_0)/(deformation radius)');
    else
        if xnorm == run.eddy.vor.dia(1)/2
            xlabel('(X-X_0)/(initial radius)');
        end
    end

    %axis image;
    %    ylim([-2 max(ylim)]);

    beautify(fontSize);

    hlegend = legend;
    set(hlegend, 'box', 'off');

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
