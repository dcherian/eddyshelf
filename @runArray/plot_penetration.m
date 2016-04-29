% diagnostic plots to show how much eddy penetrates slope
function [handles] = plot_penetration(runArray, ax, choices, clear_axis)

    if ~exist('ax', 'var'), ax = []; end
    if ~exist('clear_axis', 'var'), clear_axis = 0; end

    if ischar(ax)
        choices = ax;
        ax = [];
    end

    if ~exist('choices', 'var'), choices = []; end
    cmd_list = {'mark_timestamp', 'mark_slopebreak', 'mark_fittraj', ...
               'mark_resistance', 'mark_fitvel', 'all','dt'};

    [flag, ~] = parse_commands(cmd_list, choices);
    mark_timestamp = flag(1);
    mark_slopebreak = flag(2);
    mark_fittraj = flag(3);
    mark_resistance = flag(4);
    mark_fitvel = flag(5);

    if flag(6)
        mark_timestamp = 1;
        mark_slopebreak = 1;
        mark_fittraj = 0;
        mark_resistance = 0;
        mark_fitvel = 1;
    end

    % dt for mark_timestamp
    if flag(7) == 0
        dt = 100; % default
    else
        dt = flag(7);
    end

    % shelfbreak label options
    txtfactor = 0.15;
    fontsize = 18;

    corder_backup = runArray.sorted_colors;
    if clear_axis
        % if using multiple panels, then sorted_colors won't work because the axis is
        % probably already created.
        pos = ax.Position;
        ax.delete;
        warning('deleting provided axis.');
        ax = axes('position', pos);
    end

    hfig1 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2);
                %hold all; insertAnnotation('runArray.plot_penetration');
    if ~exist('ax', 'var') || isempty(ax)
        hfig2 = figure;ax1 = gca; hold all; %subplot(2,2,[1 3]); hold all;
    else
        hfig2 = gcf;
        ax1 = ax;
        axes(ax1);  hold on;
    end
    figure(hfig2); insertAnnotation('runArray.plot_penetration');
    ax2 = []; %subplot(2,2,2); hold all;
    ax3 = []; %subplot(2,2,4); hold all;

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    len = length(runArray.filter);
    colors = get(gca, 'ColorOrder');

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
            %tinds = cut_nan(vecfind(run.ndtime, 85));
            try
                [~,tinds] = run.calc_maxflux(2);
            catch ME
                tinds = vecfind(run.eddy.t, [dt:dt:max(run.eddy.t)]);
            end
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

            tsl = find_approx(my, run.bathy.xsl, 1);

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

            tsl = find_approx(mx, run.bathy.xsl, 1);

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

        htotal = [];
        axes(ax1)
        x = (mx - x0)/xnorm;
        y = (my - y0)/ynorm;

        % plot track
        handles.hgplt2(ff) = plot(x, y, 'Color', color, 'LineStyle', '-');
        names{ff} = name;

        htotal = handles.hgplt2(ff);
        % mark start of resistance
        if mark_resistance
            [~,~,tres] = run.locate_resistance;
            handles.hres(ff) = plot(x(tres), y(tres), 'x', ...
                                    'Color', color, 'Markersize', 22);
            htotal = [htotal handles.hres(ff)];
        end

        if mark_fitvel
            [tfit,~,~,BadFitFlag] = run.FitCenterVelocity;
            if ~BadFitFlag
                handles.hfit(ff) = ...
                    plot(x(tfit(1)), y(tfit(1)), 'x', 'Color', color, 'MarkerSize', 22);
                htotal = [htotal handles.hfit(ff)];
            end
            %keyboard;
            %h2 = scatter(x(tfit(2):tfit(3)), y(tfit(2):tfit(3)), 24, color, 'filled');
        end

        if run.bathy.L_slope/run.eddy.vor.dia(1) > 1 && ...
                mark_fittraj
            run.fit_traj();
            handles.htraj(ff) = plot(x(run.traj.tind), y(run.traj.tind), 'x', ...
                                     'Color', color, 'Markersize', 16);
            htotal = [htotal handles.htraj(ff)];
        end

        % mark timestamps
        if mark_timestamp & ~isempty(tinds)
            handles.htstampt(ff) = ...
                plot(x(tinds), y(tinds), '.', 'Color', color, 'MarkerSize', 26);
            htotal = [htotal handles.htstampt(ff)];
            %text(x(tinds), y(tinds), ...
            %     cellstr(num2str(ndtime(tinds)', 2)));
        end

        % mark slopebreak
        if mark_slopebreak
            handles.hslbreak(ff) = plot(x(tsl), y(tsl), 'o', 'Color', color, ...
                                        'MarkerSize', 10);
            htotal = [htotal handles.hslbreak(ff)];
        end

        handles.hprop(ff) = linkprop(htotal, 'Color');

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

    fontSize = [];
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

    % axes labels
    if ynorm == run.rrdeep
        ylabel(['(Y_{eddy} - ' y0str ')/(deformation radius)']);
        dytick = 1;
    else
        if ynorm == run.eddy.vor.dia(1)/2
            ylabel(['(Y_{eddy} - ' y0str ')/(initial radius)']);
            dytick = 1;
        else
            ylabel(['Distance from shelfbreak / Initial distance from ' ...
                    'shelfbreak']);
            dytick = 1;
        end
    end
    if xnorm == run.rrdeep
        xlabel(['(X_{eddy} - ' x0str ')/(deformation radius)']);
    else
        if xnorm == run.eddy.vor.dia(1)/2
            xlabel(['(X_{eddy} - ' x0str ')/(initial radius)']);
        end
    end

    if ynorm == 1
        axis square; return;
    end

    if run.bathy.axis == 'y'
        xlim([-12 2]);
        ylim([0 max(ylim)]);
        liney(1);
        set(gca, 'YTick', [0:dytick:max(ylim)]);
        title('Southern Coast');
        % mark shelfbreak
        limx = xlim;
        text(limx(1) + 0.15*diff(limx), 0.02, 'shelfbreak', 'Rotation', ...
             0, 'VerticalAlignment', 'Bottom');
    else
        xlim([0 max(xlim)]);
        ylim([min(ylim) -1*min(ylim)/2]);
        linex(1);
        title('Western Coast');
        % mark shelfbreak
        limy = ylim;
        text(0.06, limy(1) + 0.25*diff(limy), 'shelfbreak', 'Rotation', ...
             270, 'VerticalAlignment', 'bottom');
    end

    handles.hlegend = legend(handles.hgplt2, names, 'Location', legloc, 'Box', ...
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
