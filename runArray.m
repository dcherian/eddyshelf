classdef runArray < handle
    properties
        % folders
        folders;
        % array of run instances
        array;
        % description
        name;
        % length of array
        len;
    end
    methods
        % constructor
        function [runArray] = runArray(folders, name)
            if ~exist('name', 'var'), runArray.name = []; end

            runArray.array = runs.empty([length(folders) 0]);
            kk = 1;
            for ii = 1:length(folders)
                warning off;
                try
                    runArray.folders{kk} = folders{ii};
                    runArray.array(kk) = runs(runArray.folders{kk});
                    disp(runArray.array(kk).name)
                    kk = kk + 1
                catch ME
                    disp([folders{ii} ' did not work'])
                    disp(ME.message)
                    continue;
                end
            end
            runArray.len = kk-1;
        end

        function [] = plot_fluxes(runArray)
            figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all
            for ii=1:runArray.len
                run = runArray.array(ii);
                if isempty(runArray.name)
                    name = ...
                        num2str(run.eddy.Lgauss(run.eddy ...
                                                .tscaleind));
                else
                    name = runArray.name{ii};
                end
                subplot(2,1,1)
                hgplt = plot(run.csflux.time(1:end-2)/run.eddy.tscale, ...
                             smooth(run.csflux.west.shelf(1:end-2)/1e6, 3));
                addlegend(hgplt, name);

                indices = ~isnan(run.csflux.west.shelf);
                itrans = cumtrapz(run.csflux.time(indices), ...
                                  run.csflux.west.shelf(indices));
                index = find_approx(itrans, 0.05 * max(itrans));
                [mtrans,mind] = max(itrans)
                subplot(2,1,2)
                plot(run.csflux.time(indices)/run.eddy.tscale, ...
                     itrans);
            end
            subplot(2,1,1)
            ylabel('Flux (Sv)');
            subplot(2,1,2)
            ylabel('Total volume transported (m^3)');
            xlabel('Non-dimensional time');
        end

        function [] = plot_eddydiag(runArray)
            hfig1 = figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all

            hfig2 = figure;
            hold all;

            hfig3 = [];%figure;
                       %subplot(2,1,1); hold all;
                       %subplot(2,1,2); hold all;

            hfig4 = [];%figure;
                       %subplot(2,1,1); hold all;
                       %subplot(2,1,2); hold all;

            hfig5 = figure;
            hold all;

            hfig6 = figure;
            hold all

            for ii=1:runArray.len
                run = runArray.array(ii);
                ndtime = run.eddy.t/run.eddy.tscale * 86400;
                if isempty(runArray.name)
                    name = num2str(run.params.eddy.depth);
                else
                    name = runArray.name{ii};
                end
                figure(hfig1)
                subplot(2,1,1)
                hgplt = plot(ndtime, run.eddy.V./run.eddy.Lgauss./sqrt(1e-5));
                addlegend(hgplt, name);
                subplot(2,1,2)
                hgplt = plot(ndtime, run.eddy.Lgauss);

                Ro = run.eddy.V ./ run.eddy.Ls / run.params.phys.f0;

                if hfig2
                    figure(hfig2)
                    hgplt = plot(ndtime, run.eddy.prox);
                    addlegend(hgplt, name)
                end
                if hfig3
                    try
                        figure(hfig3)
                        subplot(2,1,1)
                        hgplt = plot(ndtime, run.eddy.KE./run.eddy.vol);
                        addlegend(hgplt, name);
                        subplot(2,1,2)
                        hgplt = plot(ndtime, run.eddy.PE./run.eddy.vol);
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
                %figure(hfig5)
                %hgplt = plot(ndtime, run.eddy.vol/run.eddy.vol(1));
                %addlegend(hgplt, name);

                figure(hfig6)
                try
                    if isempty(run.vorsurf)
                        %                            run.vorsurf = dc_roms_read_data([run.dir ...
                        %                    '/ocean_vor.nc'], 'rv', ...
                        %                                [], {'z' Inf Inf}, ...
                        %                                [], run.rgrid, ...
                        %                                'his', ...
                        %                                'single');
                        disp('reading pv...');
                        tic;
                        run.vorsurf = single(squeeze(ncread([run.dir ...
                                            '/ocean_vor.nc'], ...
                                                            'rv', [1 1 run.rgrid.N-1 1], ...
                                                            [Inf Inf 1 ...
                                            Inf])));
                        try
                            scale = ncreadatt([run.dir '/ocean_vor.nc'], 'rv', ...
                                              'scale_factor');
                            run.vorsurf = run.vorsurf .* scale;
                        catch ME
                        end
                        toc;
                    end
                    imx = vecfind(run.rgrid.x_rho(1,:), run.eddy.mx);
                    imy = vecfind(run.rgrid.y_rho(:,1), run.eddy.my);

                    rvcen = [];
                    for kk=1:size(run.vorsurf, 3)
                        rvcen(kk) = run.vorsurf(imx(kk), imy(kk), kk);
                    end

                    hgplt = plot(ndtime, rvcen);
                    addlegend(hgplt, name);
                catch ME
                    disp(ME.cause);
                    disp(run.name);
                end
            end
            if hfig1
                figure(hfig1);
                subplot(2,1,1)
                ylabel('Diameter (km)');
                subplot(2,1,2)
                ylabel('Isobath of center');
                xlabel('Non dimensional time');
            end

            if hfig3
                figure(hfig3);
                subplot(2,1,1); ylabel('KE/vol');
                subplot(2,1,2); ylabel('PE/vol');
            end
            if hfig4
                figure(hfig4);
                subplot(2,1,1); ylabel('KE');
                subplot(2,1,2); ylabel('KE/PE');
            end
            if hfig5
                figure(hfig5);ylabel('eddy volume (m^3)');
            end
            if hfig6
                figure(hfig6); ylabel('max surf RV');
            end
        end
            
        function [] = plot_param(runArray)
            hfig1 = figure;
            hold all
            hfig2 = figure;
            hold all
            for ii=1:runArray.len
                run = runArray.array(ii);
                if isempty(runArray.name)
                    name = run.name
                else
                    name = runArray.name{ii}
                end
                eddy_ndtime = run.eddy.t/run.eddy.tscale*86400;
                csflx_ndtime = run.csflux.time/run.eddy.tscale * 86400;
                etind = find_approx(eddy_ndtime, 1.5, 1)
                cstind = find_approx(csflx_ndtime, 1.5, 1);

                meanprox(ii) = nanmean(run.eddy.hcen(etind:end));
                meanflux(ii) = nanmean(run.csflux.west.shelf(cstind: ...
                                                             end));
                meanLz(ii) = nanmean(run.eddy.Lgauss(etind:end));
                meancy(ii) = nanmean(run.eddy.cy(etind:end));

                param(ii) = run.params.nondim.eddy.Ro/ ...
                    run.params.nondim.S_sl;

                figure(hfig1);
                hgplt = plot(param(ii), meanprox(ii), '.', 'MarkerSize', 16);
                addlegend(hgplt, name);
                disp(['run = ', run.name , ' | mean prox. = ', ...
                      num2str(meanprox(ii))]);
                %    pause;

                figure(hfig2);
                hgplt = plot(param(ii), meanflux(ii), '.', 'MarkerSize', 16);
                addlegend(hgplt, name);
                disp(['run = ', run.name , ' | mean prox. = ', ...
                      num2str(meanflux(ii))]);
            end
            figure(hfig1);
            ylabel('Proximity to shelfbreak (m)');
            xlabel('Slope parameter, Ro/S');
            figure(hfig2);
            ylabel('meandist flux');
            xlabel('Slope parameter, Ro/S');
        end

        function [] = plot_test1(runArray)
            for ii=1:runArray.len
                run = runArray.array(ii);
                figure
                ndtime = run.eddy.t ./ run.eddy.tscale * 86400;
                tind = find_approx(ndtime, 1.2, 1)
                pcolorcen(run.vorsurf(:,:,tind)'./5e-5);
                colorbar;
                caxis([-0.5 0.5]);
                title([run.name ' | ' runArray.name{ii}]);
            end

        end

        function [] = plot_test2(runArray)
            figure;
            hold all
            for ii=1:runArray.len
            end
        end
    end
end
