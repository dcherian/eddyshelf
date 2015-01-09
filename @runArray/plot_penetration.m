% diagnostic plots to show how much eddy penetrates slope
function [] = plot_penetration(runArray)

    mark_timestamp = 0;

    hfig1 = []; %figure;
    %subplot(2,1,1); hold all;
    %subplot(2,1,2); hold all;

    hfig2 = figure;
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

        ndtime = run.eddy.t * 86400 / run.tscale;
        if mark_timestamp
            tinds = vecfind(ndtime, [0.5:0.5:(max(ndtime))]);
        else
            tinds = [];
        end

        if ~isempty(hfig1)
            figure(hfig1);
            subplot(2,1,1);
            hgplt = plot(ndtime, run.eddy.vor.ne/1000 - run.bathy.xsb/1000);
            addlegend(hgplt, name);

            subplot(2,1,2);
            plot(ndtime, run.eddy.vor.se/1000 - run.bathy.xsb/ ...
                 1000);
        end

        figure(hfig2)
        axes(ax1)
        xnorm = run.rrdeep;
        ynorm = run.rrdeep; run.eddy.my(1) - run.bathy.xsb;

        x = (run.eddy.mx - run.eddy.mx(1))/xnorm;
        y = (run.eddy.my - run.bathy.xsb)/ynorm;

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
    else
        ylabel(['Distance from shelfbreak / Initial distance from ' ...
                'shelfbreak']);
    end
    if xnorm == run.rrdeep
        xlabel('(X-X_0)/(deformation radius)');
    end

    %axis image;
    ylim([0 max(ylim)]);
    set(gca, 'YTick', [0:1:max(ylim)]);
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
end
