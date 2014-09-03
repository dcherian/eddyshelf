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

            runArray.array = runs.empty([length(folders) 0]);
            kk = 1;
            for ii = 1:length(folders)
                warning off;
                try
                    runArray.folders{kk} = folders{ii};
                    runArray.array(kk) = runs(runArray.folders{kk});
                    disp([runArray.array(kk).name ' completed'])
                    kk = kk + 1;
                catch ME
                    disp([folders{ii} ' did not work'])
                    disp(ME.message)
                    continue;
                end
                if ~exist('name', 'var') || isempty(name)
                    runArray.name{ii} = runArray.array(ii).name;
                else
                    runArray.name = name;
                end
            end
            runArray.len = kk-1;
        end

        function [] = test_hashes(runArray)
            for ii=1:runArray.len
                if ~strcmpi(runArray.array(ii).csflux.hash, ...
                    'ee34764138b91a2d150b58c7791bc60d480847e1')
                    if ~strcmpi(runArray.array(ii).csflux.hash, ...
                                '2a76dc848f7ca33a4d6953ce79451e72293c72ee')
                        warning([runArray.array(ii).name ' does not ' ...
                                 'have most recent flux ' ...
                                 'calculated']);
                    end
                end
            end
        end

        function [] = print_params(runArray, command)
            for ii=1:runArray.len
                out = eval(['runArray.array(ii).' command]);
                if ~ischar(out)
                    out = num2str(out);
                end
                disp([runArray.array(ii).name ' | ' out]);
            end
        end

        function [] = plot_fluxes(runArray)
            hfig1 = figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all

            hfig2 = figure;
            subplot(2,2,1); hold all
            subplot(2,2,2); hold all
            subplot(2,2,3); hold all
            subplot(2,2,4); hold all

            hfig3 = figure;
            hold all;
            for ii=1:runArray.len
                run = runArray.array(ii);

                if run.params.flags.flat_bottom
                    continue;
                end

                if isempty(runArray.name)
                    %                    name = ...
                    %    num2str(run.eddy.Lgauss(run.eddy ...
                    %                            .tscaleind));
                    name = run.name;
                else
                    name = runArray.name{ii};
                end

                tind = find_approx(run.eddy.t*86400 / run.eddy.tscale, ...
                                   1.5, 1);
                Ue = run.eddy.V(tind);
                He = run.bathy.hsb; %run.eddy.Lgauss(tind);
                Le = run.eddy.vor.dia(tind)/2;

                fluxscl = 1;Ue * Le * He;
                transscl = 1;fluxscl * (2*Le/Ue);

                figure(hfig1)
                subplot(2,1,1)
                hgplt = plot(run.csflux.time(1:end-2)/run.eddy.tscale, ...
                             smooth((run.csflux.west.shelf(1:end-2, ...
                                                           1))/fluxscl, 3));
                addlegend(hgplt, name, 'NorthWest');

                subplot(2,1,2)
                plot(run.csflux.time/run.eddy.tscale, ...
                     run.csflux.west.itrans.shelf(:,1)/transscl);

                % total transport
                ttrans = max(abs(run.csflux.west.itrans.shelf(:,1)));
                figure(hfig2)
                try
                    % find location of center at t=1.5 (arbitrary
                    % choice)
                    xnd = (run.eddy.vor.se(tind) - ...
                           run.csflux.west.shelfwater.bins)./(run.eddy.vor.dia(tind)/2);
                    subplot(2,2,1)
                    hgplt = plot(xnd .* run.rrshelf/1000, run.csflux.west.shelfwater.itrans);
                    addlegend(hgplt, name, 'NorthEast');

                    subplot(2,2,2)
                    hgplt = plot(xnd, run.csflux.west.shelfwater.itrans ...
                                 ./ttrans);

                    addlegend(hgplt, name, 'NorthEast');

                    subplot(2,2,3)
                    plot(run.csflux.west.shelfwater.vertitrans, ...
                         run.csflux.west.shelfwater.vertbins);

                    subplot(2,2,4)
                    plot(run.csflux.west.shelfwater.vertitrans./ttrans, ...
                         run.csflux.west.shelfwater.vertbins ./ ...
                         run.bathy.hsb);
                catch ME
                    disp(ME)
                end

                figure(hfig3)
                hgplt = plot(run.csflux.time/run.eddy.tscale, ...
                             run.csflux.west.shelfwater.envelope);
                addlegend(hgplt, name);
            end
            figure(hfig1)
            subplot(2,1,1)
            ylabel('Flux (Sv)');
            liney(0, [], [1 1 1]*0.75);
            subplot(2,1,2)
            ylabel('Total volume transported');
            xlabel('Non-dimensional time');

            figure(hfig2)
            subplot(2,2,1)
            ylabel('Total volume transported (m^3)');
            xlabel('cs location (km)');

            subplot(2,2,2)
            limy = ylim;
            ldefbt = sqrt(9.81 * run.bathy.hsb)/run.params.phys.f0;
            xnd = xnd .* run.rrshelf ./ ldefbt;
            plot(xnd(2:end-1), 1./xnd(2:end-1).^2, 'color', [1 1 ...
                                1]*1);
            ylim(limy)
            ylabel('Normalized volume transported');
            xlabel('cs location / shelfbreak rossby radius');

            subplot(2,2,3)
            xlabel('Total volume transported (m^3)');
            ylabel('Vertical bin (m)');

            subplot(2,2,4)
            xlabel('Normalized volume transported');
            ylabel('Vertical bin / Shelfbreak depth');

            figure(hfig3)
            xlabel('Non-dimensional time');
            ylabel('most on-shore source of water');

        end

        function [] = plot_eddydiag(runArray)
            hfig1 = figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all

            hfig2 = figure;
            hold all;

            hfig3 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

            hfig4 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

            hfig5 = figure;
            subplot(3,1,1); hold all;
            subplot(3,1,2); hold all;
            subplot(3,1,3); hold all;

            hfig6 = []; %figure; hold all

            hfig7 = figure;
            subplot(3,1,1); hold all;
            subplot(3,1,2); hold all;
            subplot(3,1,3); hold all;

            for ii=1:runArray.len
                run = runArray.array(ii);
                ndtime = run.eddy.t/run.eddy.tscale * 86400;
                if isempty(runArray.name)
                    name = run.name;
                else
                    name = runArray.name{ii};
                end
                figure(hfig1)
                asp = run.eddy.Lgauss./(run.eddy.vor.dia/2);
                subplot(2,1,1)
                hgplt = plot(ndtime, asp./run.bathy.sl_slope);
                addlegend(hgplt, name);
                subplot(2,1,2)
                hgplt = plot(ndtime, asp);
                addlegend(hgplt, num2str(asp(1)));

                if hfig2
                    figure(hfig2)
                    hgplt = plot(ndtime, (run.eddy.my - run.bathy.xsb)/run.rrdeep);
                    addlegend(hgplt, name)
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
                end
                subplot(3,1,2)
                hgplt = plot(ndtime, sqrt(Bu));
                addlegend(hgplt, name);
                subplot(3,1,3)
                if ~isfield(run.params.nondim, 'S_sl')
                    run.params.nondim.S_sl = 0;
                end
                plot(ndtime, run.params.nondim.S_sl * ...
                     sqrt(run.params.phys.N2) * run.eddy.Lgauss ./ run.eddy.V);

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
                    hgplt = plot(ndtime, run.eddy.Ls);
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
                subplot(3,1,1); ylabel('Ro');
                subplot(3,1,2); ylabel('ND/fL'); liney(0.5);
                subplot(3,1,3); ylabel(['\alpha_{bot}/\alpha_{iso} ' ...
                                    '* Bu/Ro']);
            end
            if hfig6
                figure(hfig6); ylabel('max surf RV');
            end
            if hfig7
                figure(hfig7)
                subplot(3,1,1);
                ylabel('Radius');
                subplot(3,1,2);
                ylabel('Vertical scale (m)');
                xlabel('Non-dim time');
                subplot(3,1,3);
                ylabel('Water depth at center');
                xlabel('Non-dim time');
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
            figure;
            subplot(2,1,1); hold all;
            subplot(2,1,2); hold all;
            for ii=1:runArray.len
                run = runArray.array(ii);
                %figure
                %ndtime = run.eddy.t ./ run.eddy.tscale * 86400;
                %tind = find_approx(ndtime, 1.2, 1)
                %pcolorcen(run.vorsurf(:,:,tind)'./5e-5);
                %colorbar;
                %caxis([-0.5 0.5]);
                %title([run.name ' | ' runArray.name{ii}]);
                tind = find_approx(run.eddy.t*86400 / run.eddy.tscale, ...
                                   1.5, 1);
                N = sqrt(run.params.phys.N2);
                f0 = run.params.phys.f0;
                Le = N/f0 * run.bathy.hsb;
                Le = run.eddy.vor.dia(tind)/2;
                xnd = (run.eddy.vor.se(tind) - ...
                       run.csflux.west.shelfwater.bins)./Le;
                ttrans = sum(run.csflux.west.shelfwater.itrans);
                subplot(2,1,1)
                hgplt = plot(xnd .* run.rrshelf/1000, ...
                             run.csflux.west.shelfwater.itrans./ttrans);
                addlegend(hgplt, run.name, 'NorthEast');

                subplot(2,1,2)
                plot(run.eddy.t * 86400/run.eddy.tscale, (run.eddy.prox)./run.rrdeep);
            end
        end

        function [] = plot_test2(runArray)
            figure;
            hold all
            for ii=1:runArray.len
                run = runArray.array(ii);

                ii

                if isempty(run.ubot) || isempty(run.vbot)
                    run.read_velbot;
                end

                % get bottom velocity scale
                if size(run.ubot,3) > length(run.eddy.V)
                    Ub = squeeze(max(max(bsxfun(@times, avg1(run.ubot(:,2: ...
                                                              end-1, 2:2: ...
                                                              end),1), ...
                                            run.eddy.vormask),[], ...
                                         1), [], 2));
                else
                    Ub = squeeze(max(max(bsxfun(@times, avg1(run.ubot(: ...
 ...
                    ,2:end-1,:)), run.eddy.vormask), [], 1), [], 2));
                end

                ndtime = run.eddy.t * 86400 / run.eddy.tscale;
                hgplt = plot(ndtime, Ub ./ run.eddy.V');
                addlegend(hgplt, run.name);
            end
        end

        function [] = plot_test3(runArray)
            figure;
            hold all
            for ii=1:runArray.len
                run = runArray.array(ii);

                ndtime = run.eddy.t * 86400 / run.eddy.tscale;
                hgplt = plot(ndtime, run.eddy.dia/1000);
                addlegend(hgplt, run.name);
            end
        end

        function [] = plot_fluxhov(runArray, axis)
        %figure;
        %    hold all

            cmap = flipud(cbrewer('div','RdBu',24));
            for ii=1:runArray.len
                run = runArray.array(ii);

                if axis == 'x'
                    tmat = repmat(run.time/86400, [size(run.csflux.shelfxt, ...
                                                        1) 1]);
                    xmat = repmat(run.rgrid.x_rho(1,2:end-1)'/1000, ...
                                  [1 length(run.time)]);
                else
                    tmat = repmat(run.time/86400, [size(run.csflux.shelfxt, ...
                                                        1) 1]);
                    xmat = repmat(run.rgrid.z_r(:,run.bathy.isb,1)', ...
                                  [1 length(run.time)]);
                end
                figure;
                hold on
                pcolorcen(xmat, tmat, run.csflux.shelfxt);
                colormap(cmap);
                colorbar; center_colorbar;
                title(run.name);
                plot(run.eddy.vor.cx/1000, run.eddy.t);
                plot(run.eddy.vor.ee/1000, run.eddy.t);
                plot(run.eddy.vor.we/1000, run.eddy.t);
                ylabel('Time (days)');
                xlabel('X (km)');
                %hgplt = plot(run.eddy.t * 86400 ./ run.eddy.tscale, ...
                %             run.eddy.prox);
                %addlegend(hgplt, run.name);
            end
        end
    end
end
