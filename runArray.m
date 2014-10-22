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
        % actual indices to plot
        filter = [];
    end
    methods
        % constructor
        function [runArray] = runArray(folders, name)

            runArray.array = runs.empty([length(folders) 0]);
            kk = 1;
            for ii = 1:length(folders)
                warning off;
                try
                    runArray.folders{kk} = ['../topoeddy/' folders{ii}];
                    runArray.array(kk) = runs(runArray.folders{kk});
                    disp([runArray.array(kk).name ' completed'])

                    if ~exist('name', 'var') || isempty(name)
                        runArray.name{kk} = runArray.array(kk).name;
                    else
                        runArray.name = name;
                    end

                    kk = kk + 1;
                catch ME
                    disp([folders{ii} ' did not work'])
                    disp(ME.message)
                    continue;
                end
            end
            runArray.len = kk-1;
        end

        function [] = print_names(runArray)
            for ii=1:runArray.len
                disp([num2str(ii) ' | ' runArray.array(ii).name]);
            end
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

        function [diags] = print_diag(runArray, name)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            diags = nan(size(runArray.filter));

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                tind = run.tscaleind;

                %%%%% slope parameter
                if strcmpi(name, 'slope param')
                    diags(ff) = run.eddy.Ro(1) ./ run.bathy.S_sl;
                    diagstr = num2str(diags(ff));
                end

                %%%%% test critical iflux hypothesis for eddy to
                %%%%% start moving northward
                if strcmpi(name, 'critical flux')
                    iflux = run.csflux.west.itrans.shelf(:,1);
                    dcy = diff(run.eddy.vor.cy);

                    ind = find(dcy(run.csflux.tscaleind:end) > 0, ...
                               1, 'first') - 1 + run.csflux.tscaleind;

                    % check ind detection
                    %figure; plot(run.eddy.vor.cy); linex(ind);

                    % get flux at ind
                    run.csflux.critrans = ...
                        run.csflux.west.itrans.shelf(ind,1);

                    diags(ff) = run.csflux.critrans;
                    diagstr = num2str(diags(ff));
                end

                % penetration
                if strcmpi(name, 'hcen')
                    hfinal = mean(run.eddy.hcen(tind:end));
                    hinit = run.eddy.hcen(1);

                    diag_h = (hinit - hfinal)./run.eddy.Lgauss(tind);

                    diag_l = (mean(run.eddy.my(tind:end) - ...
                                   run.bathy.xsb))./(run.eddy.vor.dia(tind)/2);

                    diagstr = ['h = ' num2str(diag_h) ...
                               ' | L = ' num2str(diag_l)];
                end

                %%%%% beta v/s beta_t
                if strcmpi(name, 'betas')
                    diagstr = [num2str( ...
                        run.params.phys.beta ./ run.bathy.sl_slope ...
                        ./ run.params.phys.f0 ...
                        ) ' | ' num2str(mean(run.eddy.hcen(run.tscaleind:end)))];
                end

                %%%%% shelf flux
                if strcmpi(name, 'shelf flux')
                    if ~isfield(run.csflux, 'time')
                        continue;
                    end
                    if ff == 1
                        hfig_flux = figure; hold on; hax1 = gca;
                        hfig_fluxerr = figure; hold on;
                    end

                    ind = run.eddy.tscaleind;
                    transscl = 0.075 * 9.81/run.params.phys.f0 .* ...
                               run.eddy.amp(ind).* run.bathy.hsb/1000;

                    % flux vector for applicable time
                    % convert everything to double since I'm
                    % dealing with large numbers here (time in
                    % seconds) and integrated flux (m^3)
                    fluxvec = double(smooth(run.csflux.west.shelf(run.tscaleind: ...
                                                           end,1), 6));
                    ifluxvec = double(smooth(run.csflux.west.itrans.shelf(run.tscaleind: ...
                                                           end,1), 6));
                    tvec = double(run.csflux.time(run.tscaleind:end));

                    % change origin
                    ifluxvec = (ifluxvec - ifluxvec(1));
                    tvec = (tvec - tvec(1));

                    E = [ones(size(tvec))' tvec'];

                    %%%%%%%%%%% See Wunsch(1996) pg. 116
                    % P matrix
                    x = E\ifluxvec;
                    intercept = x(1);
                    avgflux = x(2);
                    true = ifluxvec; est = intercept + avgflux .* ...
                           (tvec-tvec(1))';
                    res = true-est;
                    % (E' * E) ^-1
                    %ETEI = inv(E'*E);
                    % from http://blogs.mathworks.com/loren/2007/05/16/purpose-of-inv/
                    [Q,R] = qr(E,0);
                    S = inv(R);
                    ETEI = S*S';
                    % assuming noise vector (res) is white
                    P = ETEI * E' * var(res) * E * ETEI;
                    err = sqrt(diag(P));
                    err = err(2); % standard error

                    %%%%%%%%%%% use MATLAB regress
                    [b, bint, r, rint, stats] = ...
                        regress(ifluxvec, E);
                    avgflux = b(2);
                    err = abs(bint(2) - b(2));

                    % plot fit
                    %figure; hold all;
                    %plot(tvec/86400, true, '*');
                    %plot(tvec/86400, est); plot(tvec/86400, res); liney(0);

                    %[c,lags] = xcorr(fluxvec - mean(fluxvec(:)), 'coef');
                    %plot(lags, c); linex(0); liney(0);

                    %%%%%%%%%%% mean of instantaneous flux
                    % find number of peaks
                    %mpd = 6;
                    % crests
                    %[~,pl] = findpeaks(fluxvec, 'MinPeakDistance', mpd);
                    % troughs
                    %[~,nl] = findpeaks(-1*fluxvec, 'MinPeakDistance', mpd); ...

                    % make sure peak to trough distance is not
                    % smaller than mpd
                    %indices = sort([pl; nl]);
                    %mask = [0; diff(indices) < mpd];
                    %filtered = indices(~isnan(fillnan(indices .* ~mask,0)));
                    %dof = length(filtered) + 1; % (crude) degrees of freedom;

                    % check dof calculation
                    %figure; plot(fluxvec); linex(filtered); title(num2str(dof));pause;

                    %flx = mean(max(fluxvec/1000));
                    %flx = run.csflux.west.avgflux.shelf(1)/1000;
                    % standard deviation
                    %sdev = sqrt(1./(length(fluxvec)-1) .* sum((fluxvec - flx*1000).^2))/1000;
                    % error bounds
                    %errmean = abs(conft(0.05, dof-1) * sdev / sqrt(dof));
% $$$ % $$$
% $$$                     % check error bounds with itrans
% $$$                     hfig2 = figure;
% $$$                     set(gcf, 'renderer', 'opengl');
% $$$                     subplot(2,1,1);
% $$$                     plot(run.csflux.time/run.tscale, ...
% $$$                          run.csflux.west.shelf(:,1)/1000);
% $$$                     liney([flx-err flx flx+err]);
% $$$                     subplot(2,1,2);
% $$$                     plot(run.csflux.time/run.tscale, ...
% $$$                          run.csflux.west.itrans.shelf(:,1));
% $$$                     Ln = createLine(1, ...
% $$$                                    run.csflux.west.itrans.shelf(run.tscaleind,1), ...
% $$$                                    1, (flx-err)*1000*run.tscale);
% $$$                     L = createLine(1, ...
% $$$                                    run.csflux.west.itrans.shelf(run.tscaleind,1), ...
% $$$                                    1, flx*1000*run.tscale);
% $$$                     Lp = createLine(1, ...
% $$$                                    run.csflux.west.itrans.shelf(run.tscaleind,1), ...
% $$$                                    1, (flx+err)*1000*run.tscale);
% $$$                     hold on; drawLine(L);drawLine(Ln,'Color','g'); ...
% $$$                         drawLine(Lp,'Color','r');
% $$$                     limy = ylim; ylim([0 limy(2)]);
% $$$                     %pause;
% $$$                     try
% $$$                         close(hfig2);
% $$$                     catch ME; end
% $$$ % $$$
                    run.eddy.paramflux = avgflux;
                    diagstr = [num2str(avgflux/1000,'%.2f') '±' ...
                               num2str(err/1000,'%.2f') ' mSv | scale = ' ...
                               num2str(transscl)];

                    paramerr = avgflux/1000 - transscl;
                    figure(hfig_fluxerr)
                    plot(run.eddy.Ro(tind), paramerr, '*');

                    figure(hfig_flux);
                    errorbar(transscl, avgflux/1000, err/1000, 'x');
                    %plot(transscl, flx, 'x');
                    text(transscl, double(avgflux + err)/1000, run.name, ...
                         'Rotation', 90, 'FontSize', 12);

                    %errorbar(ff , avgflux/1000, err/1000, 'x');
                    %set(hax1, 'Xtick', [1:runArray.len]);
                    %lab = cellstr(get(hax1,'xticklabel'));
                    %lab{ff} = getname(runArray, ff);
                    %set(hax1,'xticklabel', lab);
                end

                disp([run.name ' | ' name ' = ' diagstr])
            end

            if exist('hfig_flux', 'var')
                figure(hfig_flux);
                limy = ylim;
                ylim([0 limy(2)]);
                line45; axis square;
                ylabel('Flux (mSv)');
                xlabel('Parameterization (mSv)');
                beautify([18 18 20]);
            end
        end

        function [] = plot_fluxes(runArray)
            hfig1 = figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all

            hfig2 = figure;
            hold all
            %            subplot(2,2,1); hold all
            %subplot(2,2,2); hold all
            %subplot(2,2,3); hold all
            %subplot(2,2,4); hold all

            hfig3 = figure;
            hold all;

            hfig4 = figure;
            subplot(2,1,1); hold all
            subplot(2,1,2); hold all

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);

                run = runArray.array(ii);
                name = getname(runArray, ii);

                if run.params.flags.flat_bottom
                    continue;
                end

                tind = find_approx(run.eddy.t*86400 / run.tscale, ...
                                   1.5, 1);
                Ue = run.eddy.V(tind);
                He = run.bathy.hsb; %run.eddy.Lgauss(tind);
                Le = run.eddy.vor.dia(tind)/2;

                fluxscl = 1e6;Ue * Le * He;
                transscl = 1;fluxscl * (2*Le/Ue);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHELF WATER
                figure(hfig1)
                subplot(2,1,1)
                hgplt = plot(run.csflux.time(1:end-2)/run.tscale, ...
                             smooth((run.csflux.west.shelf(1:end-2, ...
                                                           1))/fluxscl, 3));
                addlegend(hgplt, name, 'NorthWest');

                subplot(2,1,2)
                plot(run.csflux.time/run.tscale, ...
                     run.csflux.west.itrans.shelf(:,1)/transscl);

                % total transport
                ttrans = max(abs(run.csflux.west.itrans.shelf(:,1)));
                figure(hfig2)
                try
                    % find location of center at t=1.5 (arbitrary
                    % choice)
                    %                    xnd = (run.csflux.west.shelfwater.bins)./run.rrshelf;
                    %                    subplot(2,2,1)
                    %hgplt = plot(xnd, run.csflux.west.shelfwater.itrans);
                    %addlegend(hgplt, name, 'NorthEast');

                    %subplot(2,2,2)
                    %plot(xnd, run.csflux.west.shelfwater.itrans ...
                    %             ./ttrans);

                    %subplot(2,2,3)
                    %plot(run.csflux.west.shelfwater.vertitrans, ...
                    %     run.csflux.west.shelfwater.vertbins);

                    %subplot(2,2,4)
                    profile = ...
                        run.csflux.west.shelfwater.vertitrans./ ...
                        ttrans;
                    zvec = run.csflux.west.shelfwater.vertbins ./ ...
                         run.bathy.hsb;
                    bc = baroclinicity(zvec, profile);
                    hgplt = plot(profile, zvec);
                    addlegend(hgplt, [name ' | bc = ' num2str(bc,'%.3f')], 'SouthWest');
                catch ME
                    disp(ME)
                end

                figure(hfig3)
                % change from envelope to depth
                env = run.csflux.west.shelfwater.envelope;
                env(isnan(env)) = max(env);
                ind = vecfind(run.rgrid.y_rho(:,1), env);
                metric = run.bathy.h(1,ind)./run.bathy.hsb .* ...
                         (1+run.rgrid.f(run.bathy.isb,1)./run.rgrid.f(ind,1))';

                hgplt = plot(run.csflux.time/run.tscale, ...
                             (run.bathy.xsb - env)./run.rrshelf);
                addlegend(hgplt, name);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDDY WATER
                figure(hfig4)
                subplot(2,1,1)
                hgplt = plot(run.csflux.time(1:end-2)/run.tscale, ...
                             smooth((run.csflux.east.eddy(1:end-2, ...
                                                           1))/fluxscl, 3));
                addlegend(hgplt, name, 'NorthWest');

                subplot(2,1,2)
                plot(run.csflux.time/run.tscale, ...
                     run.csflux.east.itrans.eddy(:,1)/transscl);
            end

            figure(hfig1)
            subplot(2,1,1)
            ylabel('Flux of shelf water (Sv)');
            liney(0, [], [1 1 1]*0.75);
            subplot(2,1,2)
            ylabel('Total volume transported');
            xlabel('Non-dimensional time');

            figure(hfig2)
            %subplot(2,2,1)
            %ylabel('Total volume transported (m^3)');
            %xlabel('cs location (km)');

            %subplot(2,2,2)
            %limy = ylim;
            %ldefbt = sqrt(9.81 * run.bathy.hsb)/run.params.phys.f0;
            %xnd = xnd .* run.rrshelf ./ ldefbt;
            %plot(xnd(2:end-1), 1./xnd(2:end-1).^2, 'color', [1 1 ...
            %                    1]*1);
            %ylim(limy)
            %ylabel('Normalized volume transported');
            %xlabel('cs location / shelfbreak rossby radius');

            %subplot(2,2,3)
            %xlabel('Total volume transported (m^3)');
            %ylabel('Vertical bin (m)');

            %subplot(2,2,4)
            ylim([-1 0]);
            limx = xlim;
            xlim([0 limx(2)]);
            xlabel('Normalized volume transported');
            ylabel('Vertical bin / Shelfbreak depth');

            figure(hfig3)
            xlabel('Non-dimensional time');
            ylabel('isobath most on-shore source of water / h_{sb}');

            figure(hfig4)
            subplot(2,1,1);
            ylabel('Flux of eddy water (Sv)');
            subplot(2,1,2);
            ylabel('Total volume transported (m^3)');
            xlabel('Non-dimensional time');

            figure(hfig5)
            ylabel('Normalized shelf water flux');

        end

        function [] = plot_eddydiag(runArray)
        % mark crosses on eddy track at tcen
            tcen = [0 100 200 300]; % in days
            hfig1 = [];figure;
            %subplot(2,1,1); hold all
            %subplot(2,1,2); hold all

            hfig2 = figure;
            hold all;

            hfig3 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

            hfig4 = []; %figure; subplot(2,1,1); hold all; subplot(2,1,2); hold all;

            hfig5 = []; %figure;
            %subplot(3,1,1); hold all;
            %subplot(3,1,2); hold all;
            %subplot(3,1,3); hold all;

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
                continue;

                figure(hfig1)

                subplot(2,1,1)
                hgplt = plot(ndtime, asp./run.bathy.sl_slope);
                addlegend(hgplt, name);
                subplot(2,1,2)
                hgplt = plot(ndtime, asp);
                addlegend(hgplt, num2str(asp(1)));

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

        function [] = plot_param(runArray)
            hfig1 = figure;
            hold all
            %hfig2 = figure;
            %hold all

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);

                run = runArray.array(ii);
                if isempty(runArray.name)
                    name = run.name
                else
                    name = runArray.name{ii}
                end
                eddy_ndtime = run.eddy.t/run.tscale*86400;
                csflx_ndtime = run.csflux.time/run.tscale * 86400;
                etind = find_approx(eddy_ndtime, 1.0, 1);
                cstind = find_approx(csflx_ndtime, 1.0, 1);

                etind = run.tscaleind;

                meanprox(ii) = nanmean(run.eddy.hcen(etind:end));
                meanflux(ii) = nanmean(run.csflux.west.shelf(cstind: ...
                                                             end));
                meanLz(ii) = nanmean(run.eddy.Lgauss(1));
                meancy(ii) = nanmean(run.eddy.cy(etind:end));

                param(ii) = abs(log(run.eddy.Ro(1)/ ...
                                run.params.nondim.S_sl));

                x = (meanprox(ii));
                y = meanLz(ii) * sqrt(param(ii));

                figure(hfig1);
                hgplt = plot(x, y, '.', 'MarkerSize', 16);
                addlegend(hgplt, name);
                disp(['run = ', run.name , ' | mean prox. = ', ...
                      num2str(meanprox(ii))]);
                %    pause;

                %figure(hfig2);
                %hgplt = plot(param(ii), meanflux(ii), '.', 'MarkerSize', 16);
                %addlegend(hgplt, name);
                %disp(['run = ', run.name , ' | mean prox. = ', ...
                %      num2str(meanflux(ii))]);
            end
            figure(hfig1);
            ylabel('Water depth at eddy center (m)');
            xlabel('Parameterization (m) : H = D * sqrt(ln(Ro/S_\alpha))');
            axis square;
            line45;
            beautify([18 18 20]);
            %figure(hfig2);
            %ylabel('meandist flux');
            %xlabel('Slope parameter, Ro/S');

        end

        function [] = streamerstats(runArray)

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                run.plot_velsec([run.tscale/86400:50:run.time(end)/86400]);

            end
        end

        function [] = plot_test1(runArray)
        %figure;
        %    subplot(2,1,1); hold all;
        %    subplot(2,1,2); hold all;
            figure;
            %ax(1) = subplot(1,2,1); hold on; title('Flux tscale');
            % ax(2) = subplot(1,2,2); hold on; title('Eddy tscale');

            ax = gca; hold all;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                % test parameterizing max. flux
                mflux = max(run.csflux.west.shelf(:,1));

                indices = [run.eddy.tscaleind];

                for tt=1:length(indices)
                    axes(ax(tt));
                    ind = indices(tt);

                    meanprox = mean(run.eddy.vor.cy(ind:end) - ...
                                    run.eddy.vor.dia(ind:end)/2 ...
                                    - run.bathy.xsb);
                    vscale = run.eddy.V(ind) * ...
                             exp(-2/3 * (1 + meanprox/ ...
                                         (run.eddy.vor.dia(ind)/2))^3);

                    disp([run.name ' | scale = ' num2str(vscale) ' m/s | actual d.a = ' ...
                          num2str(max(run.csflux.shelfxt(:)/run.bathy.hsb)) ...
                         ' m/s']);

                    mflux_param = vscale * run.bathy.hsb * ...
                        run.eddy.vor.dia(ind)/2;
                    plot(mflux, mflux_param, '*');
                    text(mflux, mflux_param*1.1, run.name);
                    xlabel('Max flux (m^3/s)');
                    ylabel('Parameterization (m^3/s)');
                end
            end

            % draw 45 deg lines
            line45(ax);
            linkaxes(ax, 'xy');
            axes(ax(1)); beautify;
            axes(ax(2)); beautify;
        end

        function [] = plot_test2(runArray)
            figure;
            hold all

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

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
                    Ub = squeeze(max(max(bsxfun(@times, ...
                                                avg1(run.ubot(:,2:end-1,:)), ...
                                                run.eddy.vormask), [], 1), [], 2));
                end

                ndtime = run.eddy.t * 86400 / run.tscale;
                hgplt = plot(ndtime, Ub ./ run.eddy.V');
                addlegend(hgplt, run.name);
            end
        end

        function [] = plot_test3(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            figure; hold all

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = getname(runArray, ii);

                env = run.csflux.west.shelfwater.envelope;
                tind = 1;
                diagnostic = mean(run.bathy.xsb - env(tind:end));

                if run.bathy.sl_shelf ~= 0
                    beta = run.params.phys.f0 ./ max(run.bathy.h(:)) * ...
                           run.bathy.sl_shelf;
                else
                    beta = Inf; run.params.phys.beta;
                end
                param = sqrt(0.075*run.eddy.V(1)./beta);

                hgplt = plot(run.csflux.time(tind:end)/run.tscale, ...
                             (run.bathy.xsb - env(tind:end))./run.rrshelf);
                %hgplt = plot(param, diagnostic, '*');
                addlegend(hgplt, name, 'NorthWest');
           end

           for ff=1:length(runArray.filter)
               ii = runArray.filter(ff);
               run = runArray.array(ii);
               name = getname(runArray, ii);
               if run.bathy.sl_shelf ~= 0
                   beta = run.params.phys.f0 ./ max(run.bathy.h(:)) * ...
                           run.bathy.sl_shelf;
               else
                   beta = Inf; run.params.phys.beta;
               end
               Ly = sqrt(0.075*run.eddy.V(1)./beta)./run.rrshelf;
               liney(Ly, run.name);
           end
           %axis square; line45;
           beautify([18 18 20]);
        end

        function [] = plot_fluxcor(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = getname(runArray, ii);

                %vec1 = run.eddy.vor.lmaj(run.tscaleind:end)./ ...
                %       run.eddy.vor.lmin(run.tscaleind:end);
                vec1 = run.eddy.vor.lmaj(run.tscaleind:end);
                vec2 = run.csflux.west.shelf(run.tscaleind:end);

                vec1 = vec1 - mean(vec1);
                vec2 = vec2 - mean(vec2);

                [c,lags] = xcorr(vec1, vec2, 'coef');
                corrcoef(vec1, vec2)
                dt = (run.csflux.time(2)-run.csflux.time(1))/86400;

                figure;
                subplot(2,1,1)
                plot(run.eddy.t(run.tscaleind:end)*86400./run.tscale, ...
                     smooth(vec1,4)./max(vec1));
                hold on
                plot(run.csflux.time(run.tscaleind:end)/run.tscale, ...
                     vec2./max(vec2), 'Color', [1 1 1]*0.75);
                subplot(2,1,2)
                plot(lags * dt,c);
                xlabel('Lag (days)');
                linex(0); liney(0);
            end
        end

        function [] = plot_fluxhov(runArray, axis)
            if ~exist('axis', 'var') || isempty(axis)
                axis = 'z';
            end

            cmap = flipud(cbrewer('div','RdBu',24));

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                figure;
                hold on

                if axis == 'x'
                    matrix = run.csflux.shelfxt;
                    tmat = repmat(run.csflux.time/86400, [size(matrix, 1) 1]);
                    xmat = repmat(run.rgrid.x_rho(1,2:end-1)'/1000, ...
                                  [1 size(tmat, 2)]);
                else
                    matrix = run.csflux.shelfzt;
                    tmat = repmat(run.csflux.time/run.tscale, [size(matrix, 1) 1]);
                    xmat = repmat(run.rgrid.z_r(:,run.bathy.isb,1), ...
                                  [1 size(tmat, 2)]);
                    xmat = xmat./max(abs(xmat(:)));
                end

                pcolorcen(tmat, xmat, matrix./max(abs(matrix(:))));
                colormap(cmap);
                colorbar; center_colorbar;
                title(run.name);

                if axis == 'x'
                    plot(run.eddy.vor.cx/1000, run.eddy.t);
                    plot(run.eddy.vor.ee/1000, run.eddy.t);
                    plot(run.eddy.vor.we/1000, run.eddy.t);
                    ylabel('Time (days)');
                    xlabel('X (km)');
                else
                    xlim([0 3]);
                    ylim([-1 0]);
                    ylabel('z/H_sb');
                    xlabel('Time (days)');
                end
                cblabel('Normalized transport');
            end
        end

        function [] = plot_penetration(runArray)
            hfig1 = figure;
            subplot(2,1,1); hold all;
            subplot(2,1,2); hold all;

            hfig2 = figure;
            hold all;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);

                run  = runArray.array(ii);
                name = getname(runArray, ii);

                ndtime = run.eddy.t * 86400 / run.tscale;

                figure(hfig1);
                subplot(2,1,1);
                hgplt = plot(ndtime, run.eddy.vor.ne/1000 - run.bathy.xsb/1000);
                addlegend(hgplt, name);

                subplot(2,1,2);
                plot(ndtime, run.eddy.vor.se/1000 - run.bathy.xsb/1000);

                figure(hfig2)
                hgplt = plot(run.eddy.mx/1000, run.eddy.my/1000);
                addlegend(hgplt, name);
            end

            figure(hfig1);
            subplot(2,1,1)
            ylabel('Northern edge - Y_{sb}');
            subplot(2,1,2)
            ylabel('Southern edge - Y_{sb}');
            liney(0);

            figure(hfig2)
            ylabel('Y (km)');
            xlabel('X (km)');

        end
    end
end

function [name] = getname(runArray, ii)
    if isempty(runArray.name)
        name = runArray.array(ii).name;
    else
        name = runArray.name{ii};
    end
end