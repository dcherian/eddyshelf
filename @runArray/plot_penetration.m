% diagnostic plots to show how much eddy penetrates slope
function [] = plot_penetration(runArray)
    hfig1 = figure;
    subplot(2,1,1); hold all;
    subplot(2,1,2); hold all;

    hfig2 = figure;
    ax1 = subplot(2,2,[1 3]); hold all;
    ax2 = subplot(2,2,2); hold all;
    ax3 = subplot(2,2,4); hold all;

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

        ndtime = run.eddy.t * 86400 / run.csflux.tscale;
        tinds = vecfind(ndtime, [0.5:0.5:(max(ndtime))]);

        figure(hfig1);
        subplot(2,1,1);
        hgplt = plot(ndtime, run.eddy.vor.ne/1000 - run.bathy.xsb/1000);
        addlegend(hgplt, name);

        subplot(2,1,2);
        plot(ndtime, run.eddy.vor.se/1000 - run.bathy.xsb/1000);

        figure(hfig2)
        axes(ax1)
        x = (run.eddy.mx - run.eddy.mx(1))/run.rrdeep;
        y = (run.eddy.my - run.bathy.xsb)/run.rrdeep;
        hgplt = plot(x, y);
        addlegend(hgplt, name, 'NorthWest');
        plot(x(tinds), y(tinds), 'kx');
        text(x(tinds), y(tinds), ...
             cellstr(num2str(ndtime(tinds)', 2)));

        axes(ax2)
        hgplt = plot(ndtime, run.eddy.Lgauss./run.eddy.hcen');

        axes(ax3)
        hgplt = plot(ndtime, run.eddy.Lgauss);
    end

    fontSize = [16 16 18];
    figure(hfig1);
    subplot(2,1,1)
    ylabel('Northern edge - Y_{sb}');
    beautify(fontSize);
    subplot(2,1,2)
    ylabel('Southern edge - Y_{sb}');
    liney(0);
    beautify(fontSize);

    figure(hfig2)
    axes(ax1)
    ylabel('(Y - Y_{sb})/(deformation radius)');
    axis image;
    ylim([0 max(ylim)]);
    set(gca, 'YTick', [0:1:max(ylim)]);
    liney(1);
    xlabel('(X-X_0)/(deformation radius)');
    beautify(fontSize);
    axes(ax2)
    ylabel('Eddy vertical scale / water depth');
    ylim([0 0.6]);
    xlabel('Non-dimensional time');
    beautify;
    axes(ax3);
    ylabel('Eddy vertical scale (m)');
    beautify;

end
