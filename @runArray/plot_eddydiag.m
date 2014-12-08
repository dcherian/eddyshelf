% plot eddy diagnostics

function [] = plot_eddydiag(runArray)
% mark crosses on eddy track at tcen
    tcen = [0 100 200 300]; % in days
    hfig1 = figure; subplot(2,1,1); hold all subplot(2,1,2); hold all

    hfig2 = figure; hold all;

    hfig3 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

    hfig4 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

    hfig5 = figure; subplot(3,1,1); hold all; subplot(3,1,2); hold  all;
    subplot(3,1,3); hold all;

    hfig6 = []; %figure; hold all

    hfig7 = []; %figure;
                %subplot(3,1,1); hold all;
                %subplot(3,1,2); hold all;
                %subplot(3,1,3); hold all;

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);

        run = runArray.array(ii);
        asp = run.eddy.Lgauss./(run.eddy.vor.dia/2);

        tscale = run.tscale; find_approx(asp, 0.5, 1);
        ndtime = run.eddy.t/tscale * 86400;
        if isempty(runArray.name)
            name = run.name;
        else
            name = runArray.name{ii};
        end

        if hfig2
            figure(hfig2)
            hgplt = plot(run.eddy.mx/1000, (run.eddy.my - run.bathy.xsb)/run.rrdeep);
            addlegend(hgplt, name);
            clr = get(hgplt, 'color');
            for ii=1:length(tcen)
                tind = find_approx(run.eddy.t, tcen(ii), 1);
                plot(run.eddy.mx(tind)/1000, ...
                     (run.eddy.my(tind)-run.bathy.xsb)/run.rrdeep, ...
                     'x', 'Color', clr, 'MarkerSize', 12);
            end
        end

        if ~isempty(hfig1)
            figure(hfig1)

            subplot(2,1,1)
            hgplt = plot(ndtime, run.eddy.Lgauss);
            addlegend(hgplt, name);
            subplot(2,1,2)
            hgplt = plot(ndtime, asp);
            addlegend(hgplt, num2str(asp(1)));
        end

        if hfig3
            try
                figure(hfig3)
                subplot(2,1,1)
                hgplt = plot(ndtime, run.eddy.KE./run.eddy.KE(1));
                addlegend(hgplt, name);
                subplot(2,1,2)
                hgplt = plot(ndtime, run.eddy.PE./run.eddy.PE(1));
                addlegend(hgplt, name);
                figure(hfig4)
                subplot(2,1,1)
                hgplt = plot(ndtime, run.eddy.KE);
                addlegend(hgplt, name);
                subplot(2,1,2)
                hgplt = plot(ndtime, run.eddy.KE./run.eddy.PE);
                addlegend(hgplt, name);
            catch ME
            end
        end
        figure(hfig5)
        Bu = (sqrt(run.params.phys.N2) * run.eddy.Lgauss ./ run.eddy.Ls / ...
              run.params.phys.f0).^2;
        subplot(3,1,1)
        try
            hgplt = plot(ndtime, run.eddy.Ro);
            addlegend(hgplt, name);
        catch ME
            warning('didn''t plot');
        end
        subplot(3,1,2)
        hgplt = plot(ndtime, run.eddy.Ls/1000);
        addlegend(hgplt, name);
        subplot(3,1,3)
        if ~isfield(run.params.nondim, 'S_sl')
            run.params.nondim.S_sl = 0;
        end
        plot(ndtime, run.eddy.V);
        %plot(ndtime, run.params.nondim.S_sl * ...
        %     sqrt(run.params.phys.N2) * run.eddy.Lgauss ./ run.eddy.V);

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

        if hfig7
            figure(hfig7)
            subplot(3,1,1)
            hgplt = plot(ndtime, run.eddy.Ls/1000);
            addlegend(hgplt, name);
            subplot(3,1,2)
            plot(ndtime, run.eddy.Lgauss);
            subplot(3,1,3)
            plot(ndtime, run.eddy.hcen);
        end
    end
    if hfig1
        figure(hfig1);
        subplot(2,1,1)
        ylabel(['aspect ratio / Bottom slope = \alpha_{iso} ' ...
                '/ \alpha_{bot}']);
        subplot(2,1,2)
        ylabel('Aspect ratio = \alpha_{iso}');
        xlabel('Non dimensional time');
    end

    if hfig2
        figure(hfig2)
        ylabel(['Center - X_{sb} / Deep water rossby ' ...
                'radius']);
        xlabel('Non dimensional time');
        title(['Crosses at [' num2str(tcen) '] days']);
        hleg = legend;
        set(hleg, 'Location', 'Northwest');
        beautify([18 18 20]);
    end

    if hfig3
        figure(hfig3);
        subplot(2,1,1); ylabel('KE/KE(1)');
        subplot(2,1,2); ylabel('PE/PE(1)');
    end
    if hfig4
        figure(hfig4);
        subplot(2,1,1); ylabel('KE');
        subplot(2,1,2); ylabel('KE/PE');
    end
    if hfig5
        figure(hfig5);
        subplot(3,1,1); ylabel('Surface Ro = <v_x>/f');
        subplot(3,1,2); ylabel('Ls (km)');
        subplot(3,1,3); ylabel('U (m/s)');
        %ylabel(['\alpha_{bot}/\alpha_{iso} ' ...
        %                    '* Bu/Ro']);
    end
    if hfig6
        figure(hfig6); ylabel('max surf RV');
    end
    if hfig7
        figure(hfig7)
        subplot(3,1,1);
        ylabel('Radius (km)');
        subplot(3,1,2);
        ylabel('Vertical scale (m)');
        xlabel('Non-dim time');
        subplot(3,1,3);
        ylabel('Water depth at center');
        xlabel('Non-dim time');
    end
end
