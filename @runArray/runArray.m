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
        function [runArray] = runArray(folders, name, reset)

            if ~exist('reset', 'var'), reset = 0; end

            runArray.array = runs.empty([length(folders) 0]);
            kk = 1;
            for ii = 1:length(folders)
                warning off;
                try
                    runArray.folders{kk} = ['../topoeddy/' folders{ii}];
                    runArray.array(kk) = runs(runArray.folders{kk}, ...
                                              reset);
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

            plots = 1;

            diags = nan(size(runArray.filter));

            if plots
                hfig = gcf;%hfig = figure;
                hold all;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                tind = run.tscaleind;

                diagstr = [];

                %%%%% slope parameter
                if strcmpi(name, 'slope param')
                    diags(ff) = run.eddy.Ro(1) ./ run.bathy.S_sl;
                end

                %%%% resistance
                if strcmpi(name, 'resistance')
                    [xx,yy,tind] = run.locate_resistance;

                    loc = 'cen';
                    %loc = 'edge';

                    % save for later
                    run.eddy.res.xx = xx;
                    run.eddy.res.yy = yy;
                    run.eddy.res.tind = tind;
                    run.eddy.res.comment = ['(xx,yy) = center ' ...
                                        'location | tind = time index'];

                    plotx = run.eddy.Ro(1)./run.bathy.S_sl;
                    %diags(ff) = yy./run.rrdeep;
                    if strcmpi(loc, 'cen')
                        hdiag = run.eddy.hcen(tind);
                        laby = 'H_{cen}./H_{eddy}';
                    else
                        if run.bathy.axis == 'y'
                            xind = find_approx(run.rgrid.y_rho(:,1), ...
                                               run.eddy.vor.se(tind), ...
                                               1);
                            hdiag = run.bathy.h(1,xind)
                        else
                            xind = find_approx(run.rgrid.x_rho(1,:), ...
                                               run.eddy.vor.we(tind), ...
                                               1);
                            hdiag = run.bathy.h(xind,1)
                        end
                        laby = 'H_{edge}./H_{eddy}';
                    end

                    diags(ff) = hdiag ./ run.eddy.Lgauss(tind);

                    % name points with run names on summary plot
                    name_points = 1;
                    labx = 'Ro/S_\alpha';
                end

                %%%%% topography parameter - taylor column
                if strcmpi(name, 'hogg')
                    diags(ff) = 1./(run.eddy.Ro(1)) ./ ...
                           (1+ run.bathy.hsb/run.eddy.Lgauss(1) * ...
                            1./run.bathy.S_sl);
                end

                %%%%% test critical iflux hypothesis for eddy to
                %%%%% start moving northward
                if strcmpi(name, 'critical flux')
                    iflux = run.csflux.west.itrans.shelf(:,1);
                    dcy = diff(run.eddy.vor.cy);

                    % make sure my time vectors are the same
                    assert(isequal(run.eddy.t*86400, ...
                                   run.csflux.time));

                    % find where velocity changes sign
                    ind = find(dcy(run.csflux.tscaleind:end) > 0, ...
                               1, 'first') - 1 + run.csflux.tscaleind;

                    % check ind detection with flux and center
                    % location plot
                    %figure;
                    %plotyy(run.eddy.t, run.eddy.vor.cy, ...
                    %       run.csflux.time/86400, ...
                    %       run.csflux.west.shelf(:,1));
                    %linex(run.eddy.t(ind));

                    %run.animate_zeta(ind, 1);

                    % Get Flux at ind
                    run.csflux.critrans = ...
                        run.csflux.west.itrans.shelf(ind,1);

                    diags(ff) = run.csflux.critrans;
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
                    betat = run.bathy.sl_slope .* run.params.phys.f0 ...
                            ./ max(run.bathy.h(:));
                    diagstr = [num2str( ...
                        run.params.phys.beta ./ betat  ...
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
                    %lab{ff} = runArray.getname( ff);
                    %set(hax1,'xticklabel', lab);
                end

                if isempty(diagstr)
                    diagstr = num2str(diags(ff));
                end

                if plots
                    figure(hfig);
                    plot(plotx, diags(ff), '*');

                    % add function call as annotation
                    insertAnnotation(['runArray.print_diag(' name ')']);
                    % add run names
                    if name_points
                        text(plotx, diags(ff), run.name, 'FontSize', ...
                             12, 'Rotation', 90);
                    end
                    xlabel(labx);
                    ylabel(laby);
                    title(name);
                end

                disp([run.name ' | ' name ' = ' diagstr])
            end

            if plots
                figure(hfig)
                beautify([18 18 20]);
            end

            if exist('hfig_flux', 'var')
                figure(hfig_flux);
                insertAnnotation(['runArray.print_diag(' name ')']);
                limy = ylim;
                ylim([0 limy(2)]);
                line45; axis square;
                ylabel('Flux (mSv)');
                xlabel('Parameterization (mSv)');
                beautify([18 18 20]);
            end
        end

        function [] = plot_param(runArray)
            hfig1 = figure;
            insertAnnotation(['runArray.plot_param']);
            hold all
            hfig2 = figure;
            insertAnnotation(['runArray.plot_param']);
            hold all

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);

                run = runArray.array(ii);
                if isempty(runArray.name)
                    name = run.name;
                else
                    name = runArray.name{ii};
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

                param(ii) = (run.eddy.Ro(1)/ run.params.nondim.S_sl);

                x = (meanprox(ii));
                y = meanLz(ii) * sqrt(abs(log(param(ii))));

                figure(hfig1);
                hgplt = plot(x, y, '.', 'MarkerSize', 16);
                addlegend(hgplt, name);
                disp(['run = ', run.name , ' | mean prox. = ', ...
                      num2str(meanprox(ii))]);
                %    pause;

                figure(hfig2);
                hgplt = plot(param(ii), meanprox(ii), '.', 'MarkerSize', ...
                             16);
                text(param(ii), meanprox(ii), run.name)
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

        function [] = plot_jetprops(runArray)
            figure;
            ax = gca; hold all;
            insertAnnotation('runArray.jetprops');
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                ndtime = run.eddy.t*86400 ./ run.csflux.tscale;

                hplot = plot(ndtime, run.jet.vscale);
                addlegend(hplot, run.name);
            end
            linkaxes(ax, 'x');
            xlim([0 4])
        end

        function [] = plot_test1(runArray)
            figure;
            ax = gca; hold all;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                index = 1;

                name = [run.name ' | ' num2str(run.asflux.x(index)/1000 ...
                                               ) ' km'];
                ndtime = run.eddy.t*86400 ./ run.csflux.tscale;
                tind = 1:length(ndtime);

                hplot = plot(ndtime(tind), run.asflux.irflux(tind,index));
                %hplot = plot(run.eddy.cx, run.eddy.cy);
                addlegend(hplot, name);
            end
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
                name = runArray.getname( ii);

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
               name = runArray.getname( ii);
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
                name = runArray.getname( ii);

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

        function [name] = getname(runArray, ii)
            if isempty(runArray.name)
                name = runArray.array(ii).name;
            else
                name = runArray.name{ii};
            end
        end
    end
end
