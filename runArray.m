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
            hold all
            for ii=1:runArray.len
                run = runArray.array(ii);
                hgplt = plot(run.csflux.time/run.eddy.tscale, ...
                             smooth(run.csflux.west.shelf/1e6, 3));
                addlegend(hgplt, run.name);
            end
        end

        function [] = plot_eddydiag(runArray)
            hfig1 = figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all

            hfig2 = figure;
            hold all;
            for ii=1:runArray.len
                run = runArray.array(ii);
                ndtime = run.eddy.t/run.eddy.tscale * 86400;
                
                figure(hfig1)
                subplot(2,1,1)
                hgplt = plot(ndtime, run.eddy.vor.dia/1000);
                addlegend(hgplt, run.name);
                subplot(2,1,2)
                hgplt = plot(ndtime, run.eddy.hcen);

                Ro = run.eddy.V ./ run.eddy.Ls / run.params.phys.f0;
                figure(hfig2)
                hgplt = plot(ndtime, Ro);
                addlegend(hgplt, run.name)
            end
            figure(hfig1);
            subplot(2,1,1)
            ylabel('Diameter (km)');
            subplot(2,1,2)
            ylabel('Isobath of center');
            xlabel('Non dimensional time');
        end

        function [] = plot_param(runArray)
            hfig1 = figure;
            hold all
            hfig2 = figure;
            hold all
            for ii=1:runArray.len
                run = runArray.array(ii);
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
                addlegend(hgplt, run.name);
                disp(['run = ', run.name , ' | mean prox. = ', ...
                      num2str(meanprox(ii))]);
                %    pause;

                figure(hfig2);
                hgplt = plot(param(ii), meanflux(ii), '.', 'MarkerSize', 16);
                addlegend(hgplt, run.name);
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
    end
end
