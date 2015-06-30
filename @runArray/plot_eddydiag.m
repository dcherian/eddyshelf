% plot eddy diagnostics

function [] = plot_eddydiag(runArray)

    annostr = 'plot_eddydiag(runArray)';

    % mark crosses on eddy track at tcen
    tcen = [0 100 200 300]; % in days

    corder_backup = runArray.sorted_colors;

    % vertical scale
    hfig1 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all

    % center track
    hfig2 = []; %figure; hold all;

    % KE, PE - normalized by initial value
    hfig3 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

    % actual KE, PE
    hfig4 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

    % Ro, L-scale, U-scale
    hfig5 = figure; subplot(3,1,1); hold all; subplot(3,1,2); hold  all;
            subplot(3,1,3); hold all;

    hfig6 = []; %figure; hold all

    hfig7 = []; %figure;
                %subplot(3,1,1); hold all;
                %subplot(3,1,2); hold all;
                %subplot(3,1,3); hold all;

    hfig8 = []; %figure; hold all;

    % x,y velocities
    hfig9 = figure; subplot(211); hold all; subplot(212); hold all;

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

        run = runArray.array(ii);
        %asp = run.eddy.hcen./(run.eddy.vor.dia/2);

        %tscale = run.tscale; find_approx(asp, 0.5, 1);
        %ndtime = run.eddy.t/tscale * 86400;
        run.fit_traj(1.0);
        [~,~,tind] = run.locate_resistance;
        ndtime = run.eddy.t*86400 ./ run.eddy.turnover;
        if isempty(runArray.name)
            name{ff} = run.name;
        else
            name{ff} = runArray.name{ii};
        end

        % center-track
        if ~isempty(hfig2)
            xnorm = run.rrdeep; (run.eddy.vor.dia(1)/2);
            ynorm = run.rrdeep; (run.eddy.vor.dia(1)/2);
            if runArray.rotate_ns && run.bathy.axis == 'x'
                y = (run.eddy.mx - run.bathy.xsb)/ynorm;
                x = (run.eddy.my - run.eddy.my(1))/xnorm;
            else
                x = (run.eddy.mx - run.eddy.mx(1))/xnorm;
                y = (run.eddy.my - run.bathy.xsb)/ynorm;
            end
            figure(hfig2)
            hgplt2(ff) = plot(x, y);
            %set(gca, 'ColorOrder', cbrewer('seq', 'Reds', len), ...
            %         'NextPlot', 'replacechildren');

            %clr = get(hgplt, 'color');
            %for ii=1:length(tcen)
            %    tind = find_approx(run.eddy.t, tcen(ii), 1);
            %    plot(x(tind), y(tind), ...
            %         'x', 'Color', clr, 'MarkerSize', 12);
            %end
        end

        % vertical scale
        if ~isempty(hfig1)
            figure(hfig1)
            colors = get(gca, 'ColorOrder');
            subplot(2,1,1)
            hgplt1(ff) = plot(ndtime, run.eddy.Lgauss, ...
                              'Color', colors(ff,:));
            plot(ndtime(run.traj.tind), run.eddy.Lgauss(run.traj.tind), ...
                 'x', 'Color', colors(ff,:));
            subplot(2,1,2)
            ind = 1; % 0vor or SSH contour ρ criterion?
            if isfield(run.eddy, 'vol')
                plot(ndtime, run.eddy.vol(:,ind), 'Color', colors(ff,:));
                plot(ndtime(run.traj.tind), run.eddy.vol(run.traj.tind,ind), ...
                     'x', 'Color', colors(ff,:));
            end
        end

        % KE, PE
        if ~isempty(hfig3)
            try
                figure(hfig3)
                subplot(2,1,1)
                plot(ndtime, run.eddy.KE./run.eddy.vol);
                subplot(2,1,2)
                hgplt3(ff) = plot(ndtime, run.eddy.PE./run.eddy.vol);
            catch ME
            end
        end
        if ~isempty(hfig4)
            try
                Escale = 1/2 * 1000 * run.eddy.V(2)^2 * run.eddy.vol(1);
                figure(hfig4)
                subplot(2,1,1)
                plot(ndtime, run.eddy.KE./Escale);
                subplot(2,1,2)
                hgplt4(ff) = plot(ndtime, run.eddy.PE./Escale);
            catch ME
            end
        end

        Bu = (sqrt(run.params.phys.N2) * run.eddy.Lgauss ./ run.eddy.Ls / ...
              run.params.phys.f0).^2;

        % Ro, Length scale, Velocity Scale
        if ~isempty(hfig5)
            figure(hfig5)
            subplot(3,1,1)
            try
                Rh = run.eddy.V./run.params.phys.beta./run.eddy.Ls.^2;
                plot(ndtime, run.eddy.rhossh.Ro);
            catch ME
                warning('didn''t plot');
            end
            subplot(3,1,2)
            hgplt5(ff) = plot(ndtime, run.eddy.rhossh.dia/2/1000);
            subplot(3,1,3)
            if ~isfield(run.params.nondim, 'S_sl')
                run.params.nondim.S_sl = 0;
            end
            plot(ndtime, run.eddy.V);
        end
        %{figure(hfig6)
        %try
        %    if isempty(run.vorsurf)
        %                            run.vorsurf = dc_roms_read_data([run.dir ...
        %                    '/ocean_vor.nc'], 'rv', ...
        %                                [], {'z' Inf Inf}, ...
        %                                [], run.rgrid, ...
        %                                'his', ...
        %                                'single');
        %         disp('reading pv...');
        %tic;
        %run.vorsurf = single(squeeze(ncread([run.dir ...
        %                    '/ocean_vor.nc'], ...
        %                                    'rv', [1 1 run.rgrid.N-1 1], ...
        %                                    [Inf Inf 1 ...
        %                    Inf])));
        %%try
        %          scale = ncreadatt([run.dir '/ocean_vor.nc'], 'rv', ...
        %                              'scale_factor');
        %            run.vorsurf = run.vorsurf .* scale;
        %        catch ME
        %        end
        %        toc;
        %    end
        %    imx = vecfind(run.rgrid.x_rho(1,:), run.eddy.mx);
        %    imy = vecfind(run.rgrid.y_rho(:,1), run.eddy.my);

        %rvcen = [];
        %    for kk=1:size(run.vorsurf, 3)
        %        rvcen(kk) = run.vorsurf(imx(kk), imy(kk), kk);
        %    end

        %   hgplt = plot(ndtime, rvcen);
        %    addlegend(hgplt, name);
        %catch ME
        %    disp(ME.cause);
        %    disp(run.name);
        %    end
        %    %}

        % Horizontal length-scale, vertical scale, depth at center
        if ~isempty(hfig7)
            figure(hfig7)
            subplot(3,1,1)
            hgplt7(ff) = plot(ndtime, run.eddy.Ls/1000);
            subplot(3,1,2)
            plot(ndtime, run.eddy.Lgauss);
            subplot(3,1,3)
            plot(ndtime, run.eddy.hcen);
        end

        % volume
        if ~isempty(hfig8)
            vol = smooth(run.eddy.vor.lmaj./run.rrdeep .* run.eddy.vor.lmin./run.rrdeep ...
                         .* run.eddy.Lgauss./run.bathy.xsl, 1);
            tindex = find_approx(run.time./run.eddy.tscale, 1, 1);
            figure(hfig8)
            hgplt8(ff) = plot(ndtime, vol);
            plot(ndtime(tindex), vol(tindex), '.', 'MarkerSize', 24, ...
                 'Color', get(hgplt8(ff), 'Color'));
        end

        % x and y velocities
        if ~isempty(hfig9)
            figure(hfig9)
            subplot(211)
            hgplt9(ff) = plot(ndtime, run.eddy.rhovor.cvx./run.eddy.V(1));
            subplot(212)
            plot(ndtime, run.eddy.rhovor.cvy./run.eddy.V(1));
        end
    end

    drawnow; pause(1);
    if ~isempty(hfig1)
        figure(hfig1);
        subplot(2,1,1)
        ylabel(['Vertical scale (m)']);
        %ylabel(['aspect ratio / Bottom slope = \alpha_{iso} ' ...
        %        '/ \alpha_{bot}']);
        subplot(2,1,2)
        ylabel('Volume/Volume(1)');
        %ylabel('Aspect ratio = \alpha_{iso}');
        xlabel('Time /(turnover time)');
        insertAnnotation(annostr);
        legend(hgplt1, name);
    end

    if ~isempty(hfig2)
        figure(hfig2)
        ylabel(['(Center - Y_{sb}) / Deep water rossby ' ...
                'radius']);
        xlabel('(Center - X(1)) /  Deep water rossby radius');
        title(['Crosses at [' num2str(tcen) '] days']);
        hleg = legend;
        set(hleg, 'Location', 'Northwest');
        limy = ylim;
        ylim([-2 max(limy)]);
        liney(0);
        beautify([18 18 20]);
        insertAnnotation(annostr);
        legend(hgplt2, name);
    end

    if ~isempty(hfig3)
        figure(hfig3);
        subplot(2,1,1); ylabel('KE/vol');
        subplot(2,1,2); ylabel('PE/vol');
        insertAnnotation(annostr);
        legend(hgplt3, name);
    end
    if ~isempty(hfig4)
        figure(hfig4);
        subplot(2,1,1); ylabel('KE');
        subplot(2,1,2); ylabel('PE');
        legend(hgplt4, name);
    end
    if ~isempty(hfig5)
        figure(hfig5);
        subplot(3,1,1); ylabel('Surface Ro = <v_x>/f');
        subplot(3,1,2); ylabel('Ls (km)');
        subplot(3,1,3); ylabel('U (m/s)');
        insertAnnotation(annostr);
        legend(hgplt5, name);
        %ylabel(['\alpha_{bot}/\alpha_{iso} ' ...
        %                    '* Bu/Ro']);
    end
    if ~isempty(hfig6)
        figure(hfig6); ylabel('max surf RV');
        insertAnnotation(annostr);
        legend(hgplt6, name);
    end
    if ~isempty(hfig7)
        figure(hfig7)
        subplot(3,1,1);
        ylabel('Radius (km)');
        subplot(3,1,2);
        ylabel('Vertical scale (m)');
        xlabel('Non-dim time');
        subplot(3,1,3);
        ylabel('Water depth at center');
        xlabel('Non-dim time');
        insertAnnotation(annostr);
        legend(hgplt7, name);
    end
    if ~isempty(hfig8)
        figure(hfig8);
        title('Dots = eddy center crosses slopebreak');
        ylabel('L_{maj} \times L_{min} \times H (m^3)');
        xlabel('Non-dim time');
        beautify;
        insertAnnotation(annostr);
        legend(hgplt8, name);
    end
    if ~isempty(hfig9)
        figure(hfig9);
        subplot(211);
        liney(0);
        ylabel('Center x-velocity / Eddy velocity scale');
        beautify;

        subplot(212);
        liney(0);
        ylabel('Center y-velocity / Eddy velocity scale');
        xlabel('Non-dimensional time');
        beautify;
        legend(hgplt9, name);
    end

    runArray.reset_colors(corder_backup);
end
