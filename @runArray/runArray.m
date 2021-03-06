classdef runArray < handle
    properties
        % folders
        folders;
        % array of run instances
        array;
        % description
        name;
        % rotate NS track plots to align with EW?
        rotate_ns = 0;
        % sort by this parameter?
        sort_param = []; sorted = 0;
        % flip sorted colors
        flip_colors = 0;
        % length of array
        len;
        % actual indices to plot
        filter = [];
    end
    methods
        % constructor
        function [runArray] = runArray(folders, name, reset, reduced)

            ticstart = tic;
            if ~exist('reset', 'var'), reset = 0; end
            if ~exist('reduced', 'var'), reduced = 0; end

            if ~iscell(folders), folders = cellstr(folders); end

            runArray.array = runs.empty([length(folders) 0]);
            kk = 1;
            for ii = 1:length(folders)
                warning off;
                try
                    if isempty(findstr('run', folders{ii}));
                        runArray.folders{kk} = ['../topoeddy/run' folders{ii}];
                    else
                        runArray.folders{kk} = ['../topoeddy/' folders{ii}];
                    end
                    runArray.array(kk) = runs(runArray.folders{kk}, ...
                                              reset, 0, reduced);
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
            runArray.filter = [1:runArray.len];
            runArray.test_hashes;
            toc(ticstart);
        end

        function [] = add(runArray, dir)
            if ~isempty(runArray)
                len = runArray.len;
            else
                len = 0;
            end

            try
                runArray.folders{len+1} = dir;
                runArray.array(len+1) = dir;
                runArray.name{len+1} = runArray.array(len+1).name;
                disp([runArray.array(len+1).name ' completed'])
                runArray.test_hashes;
                runArray.len = len + 1;
            catch ME
                disp([dir ' did not work'])
                disp(ME.message)
            end

            runArray.filter = 1:runArray.len;
        end

        function [] = delete(runArray, index)
            input(['Delete ' runArray.array(index).name '? ']);
            runArray.folders{index} = [];
            runArray.array(index) = [];
            runArray.name{index} = [];
            runArray.len = runArray.len - 1;

            % clean up empty arrays
            runArray.name = ...
                runArray.name(~cellfun('isempty',runArray.name));
            runArray.folders = ...
                runArray.folders(~cellfun('isempty', ...
                                          runArray.folders));

            runArray.filter = 1:runArray.len;
        end

        function [] = reload(runArray, index)
            runArray.array(index) = runs(runArray.array(index).dir);
        end

        function [] = print_names(runArray)
            for ii=1:runArray.len
                disp([num2str(ii) ' | ' runArray.array(ii).name]);
            end
        end

        % sort members of the array by runArray.sort_param;
        function [] = sort(runArray, sort_input)

            if ~exist('sort_input', 'var') || isempty(sort_input)
                error('need sort_input to sort!');
            end

            if ischar(sort_input)
                sort_input = runArray.print_params(sort_input);
            end

            if length(sort_input) < runArray.len
                error('length(sort_input) < runArray.len. DO NOT FILTER!');
            end

            [ss,ind] = sort(sort_input, 'ascend');
            runArray.sort_param = sort_input;

            % sort arrays
            runArray.array = runArray.array(ind);

            % sort names
            for ii = 1:length(ind)
                names{ii} = runArray.name{ind(ii)};
            end
            runArray.name = names;

            runArray.sorted = 1;

            disp(['runArray sorted.']);
        end

        % helper function for setting line colors when plotting
        % diagnostics from a sorted runArray object
        function [corder_backup] = sorted_colors(runArray)
            corder_backup = get(groot, 'defaultAxesColorOrder');
            if isempty(runArray.filter)
                len = runArray.len;
            else
                len = length(runArray.filter);
            end

            if runArray.sorted
                colors = brighten(cbrewer('seq','Reds',len), -0.5);
                set(groot, 'defaultAxesLineStyleorder','-');
                if runArray.flip_colors
                    set(groot, 'defaultAxesColorOrder', flip(colors));
                else
                    set(groot, 'defaultAxesColorOrder', colors);
                end
            else
                set(groot, 'defaultAxesColorOrder', linspecer(len));
            end
        end

        function [] = reset_colors(runArray, corder_backup)
            if runArray.sorted
                set(0, 'DefaultAxesColorOrder', corder_backup);
                set(0,'DefaultAxesLineStyleOrder',{'-','--','-.'});
            end
        end

        function [] = test_hashes(runArray)
            warning on;
            disp('======= eddy =======');
            for ii=1:runArray.len
                try
                    msg = [runArray.array(ii).name];
                    if isempty(strfind(runArray.array(ii).eddy.hash, ...
                                       githash('./track_eddy.m')))
                        warning(msg);
                    end
                catch ME
                    warning(msg);
                end
            end

            disp('======= csflux =======');
            for ii=1:runArray.len
                try
                    msg = [runArray.array(ii).name ' is outdated: csflux'];
                    if isempty(strfind(runArray.array(ii).csflux.hash, ...
                                       githash('./@runs/csfluxes.m')))
                        warning(msg);
                    end
                catch ME
                    warning(msg);
                end
            end

            disp('======= streamer =======');
            for ii=1:runArray.len
                try
                    msg = [runArray.array(ii).name ' is outdated: streamer'];
                    if isempty(strfind(runArray.array(ii).streamer.hash, ...
                                       githash('./@runs/avgStreamerVelSection.m')))
                        warning(msg);
                    end
                catch ME
                    warning(msg);
                 end
            end
        end

        function [out] = print_params(runArray, command)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            out = [];

            for ii=1:length(runArray.filter)
                run = runArray.array(runArray.filter(ii));
                [~,~,tres] = run.locate_resistance;
                %[~,maxloc] = run.calc_maxflux(2);

                use run.params;
                eddy = run.eddy; csflux = run.csflux; bathy = run.bathy;

                beta_sl = phys.f0 * bathy.sl_slope./bathy.hsl;
                beta_sh = phys.f0 * bathy.sl_shelf./bathy.hsb;

                try
                    [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
                    t0 = start;
                    tend = stop;
                    [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);
                catch ME
                end

                temp = [];
                try
                    temp = eval(['run.' command]);
                catch ME
                    try
                        temp = eval(command);
                    catch ME
                        if iscell(out)
                            out{ii} = NaN;
                        else
                            out(ii) = NaN;
                        end
                    end
                end

                if ~isempty(temp)
                    if ischar(temp)
                        if ~iscell(out)
                            clear out
                            out{ii-1} = '';
                        end
                        out{ii} = temp;
                    else
                        out(ii) = temp;
                    end

                end

                if ~iscell(out) && (~ischar(out(ii)))
                    outstr = num2str(out(ii));
                else
                    outstr = out{ii};
                end
                disp([num2str(runArray.filter(ii), '%02d') ' | ' ...
                      run.name ' | ' outstr]);
            end
        end

        function [] = plot_energy(runArray)

            hfig1 = figure;
            insertAnnotation(['runArray.plot_energy']);

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff = 1:length(runArray.filter)
                ii = runArray.filter(ff);

                run = runArray.array(ii);
                name = runArray.getname(ii);

                ndtime = run.eddy.t * 86400 ./ run.eddy.turnover;

                try
                    subplot(2,1,1); hold all
                    hplt = plot(ndtime, run.eddy.KE);
                    title('KE');
                    addlegend(hplt, name);

                    subplot(2,1,2); hold all
                    plot(ndtime, (run.eddy.PE - run.eddy.PE(end)));
                    title('PE');
                catch ME
                end
            end

            subplot(211); linex(1);
            subplot(212); limy = ylim; ylim([0 limy(2)]);
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
                meanflux(ii) = nanmean(run.csflux.off.shelf(cstind: ...
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

        function [] = plotSSHgradient(runArray)
            figure; hax = gca; maximize; hold all;
            for ii = 1:length(runArray.filter)
                ff = runArray.filter(ii);
                hh(ii) = runArray.array(ff).plotSSHgradient('max flux', hax);
                names{ii} = runArray.array(ff).name;
            end

            legend(hh, names);
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
            hfig = figure;
            ax1 = subplot(211); hold all;
            ax2 = subplot(212); hold all;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            if runArray.sorted
                subplot(211);
                co = runArray.sorted_colors;
                subplot(212);
                co = runArray.sorted_colors;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = run.name;

                ndtime = run.eddy.t * 86400./ (run.eddy.turnover);
                tind = 1:ceil(50*run.eddy.turnover/86400);

                beta = run.params.phys.beta;
                Ldef = run.rrdeep;

                dEdt = smooth(diff((run.eddy.KE + run.eddy.PE)./run.eddy.vol)./ ...
                       diff(run.eddy.t*86400), 14);
                ndtime1 = avg1(ndtime);

                [~,~,rest] = run.locate_resistance;
                tinds = [run.eddy.tscaleind run.eddy.edgtscaleind rest];
                axes(ax1)
                hplot = plot(ndtime1, dEdt);
                %figure;
                %hplot = plot(ndtime, run.eddy.cvx * 1000/86400);
                addlegend(hplot, name);
                plot(ndtime1(tinds), dEdt(tinds), 'k*');
                %liney(-beta .* Ldef^2);
                %linex(ndtime(run.tscaleind));
                %pause();

                axes(ax2)
                plot(ndtime, run.eddy.KE);
                plot(ndtime(tinds), run.eddy.KE(tinds), 'k*');
            end
            axes(ax1); liney(0);
            axes(ax2); liney(0);

            insertAnnotation('runArray.plot_test1');
            beautify;

            if runArray.sorted
                runArray.reset_colors(co);
            end
        end

        function [] = plot_Rosurf(runArray)
            hfig = figure;
            ax = gca; hold all;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                if isempty(run.vorsurf)
                    run.calc_vorsurf;
                end

                name = run.name;
                ndtime = run.eddy.t * 86400./ (run.eddy.turnover);
                tind = 1:length(ndtime);

                Ro = avg1(avg1(bsxfun(@rdivide, run.vorsurf, ...
                                      avg1(avg1(run.rgrid.f', 1), ...
                                           2)),1),2);

                Ro = Ro .* run.eddy.vor.mask;

                iy = vecfind(run.rgrid.y_rho(:,1), run.eddy.my);
                fcen = run.rgrid.f(iy,1);
                pv = fcen./run.eddy.Lgauss' .* (1+run.eddy.Ro');
                hplot = plot(pv./pv(1));

                Romin = squeeze(min(min(Ro,[],1),[],2));
                %hplot = plot(ndtime, Romin./Romin(1));
                addlegend(hplot, name);
            end

            ylabel('min(vorticity/f)');
            xlabel('Time / turnover time');
            insertAnnotation('runArray.plot_test1');
            beautify;
        end

        function [] = plot_test2(runArray)

            corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            hf = figure; hold all
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.name{ii};
                ndtime = run.eddy.t*86400 / run.eddy.turnover;

                clear vec vec2

                [~,~,tind] = run.locate_resistance;

                thresh = 0.2;
                Tnorm = abs(bsxfun(@rdivide, run.eddy.T, ...
                                   run.eddy.T(:,end)));
                for tt=1:size(Tnorm,1)
                    ind = find_approx(Tnorm(tt,:), thresh, 1);
                    vec2(tt) = run.eddy.zT(tt,ind) * -1;
                end
                vec = vec2;
                hplt = plot(ndtime, vec);
                plot(ndtime, run.eddy.hcen);
                plot(ndtime(tind), vec(tind), 'kx');
            end

            legend(hplt, names);
            runArray.reset_colors(corder_backup);
        end

        function [] = plot_test3(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            corder_backup = runArray.sorted_colors;

            figure; hold all;
            insertAnnotation('runArray.plot_test3');

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.getname(ii);
                %run.fit_traj;
                %tind = run.traj.tind;
                %[~,~,tind] = run.locate_resistance(10,1/2);
                tvec = run.time; %/run.eddy.turnover;

                lmaj = run.eddy.rhovor.lmaj;
                lmin = run.eddy.rhovor.lmin;
                e = lmin./lmaj;

                ind = 1;
                [~,tind] = run.calc_maxflux(ind);

                hplt(ff) = plot(tvec, smooth(e,4));
                plot(tvec(tind), e(tind), 'kx');
            end
            legend(hplt, names);

            runArray.reset_colors(corder_backup);
        end

        function [] = plot_hT(runArray)
        % Estimate vertical scale using ΔT at eddy center.
            corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            hf = figure; hold all
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.name{ii};
                ndtime = run.eddy.t*86400 / run.eddy.turnover;

                clear vec vec2

                [~,~,tind] = run.locate_resistance;

                thresh = 0.1;
                Tnorm = abs(bsxfun(@rdivide, run.eddy.T, ...
                                   run.eddy.T(:,end)));
                for tt=1:size(Tnorm,1)
                    ind = find_approx(Tnorm(tt,:), thresh, 1);
                    vec2(tt) = run.eddy.zT(tt,ind);
                end
                vec = vec2;
                hplt = plot(ndtime, vec);
                %plot(ndtime, run.eddy.hcen);
                plot(ndtime(tind), vec(tind), 'kx');
            end

            legend(hplt, names);
            runArray.reset_colors(corder_backup);
        end

        function [hax] = plot_fluxparam(runArray, str, factor, commands)

            if ~exist('commands', 'var'), commands = ''; end

            if strcmpi(str, 'max flux') | strcmpi(str, 'avg flux') ...
                    & (isempty(factor) | ~exist('factor', 'var'))
                factor = 1; % to 2xH_{sb}
            end

            isobath = 1:8;
            if isobath(1) == 1 & ...
                    ~(strcmpi(str, 'max flux') | strcmpi(str, 'avg flux'))
                isobath(1) = [];
            end

            commands = [commands '; no_name_points'];

            figure;
            insertAnnotation(['runArray.plot_fluxparam(' str ')']);
            slopeplot = 3;
            hax = packfig(3,3);
            for ii=1:length(isobath)
                iso = isobath(ii);
                if ii >= slopeplot, hh = ii+1; else hh = ii; end
                axes(hax(hh));
                [diags(ii,:), plotx(ii,:), err(ii,:), norm(ii,:), color, ...
                 ~, P, Perr] = runArray.print_diag(str, [iso factor], hax(hh), ...
                                                   commands);

                labx = hax(hh).XLabel.String;
                laby = hax(hh).YLabel.String;
                xlabel(''); ylabel(''); title('');
                ggplot;
                hax(hh).XColor = [1 1 1];
                hax(hh).YColor = [1 1 1];

                mplt(ii) = P(1);
                cplt(ii) = P(2);
                merr(ii) = Perr(1);
                cerr(ii) = Perr(2);
            end

            % link axes properly
            linkaxes(hax([1 2 4:9]), 'x');
            linkaxes(hax([2 1]), 'y');
            linkaxes(hax([6 5 4]), 'y');
            linkaxes(hax([9 8 7]), 'y');

            linkprop(hax([1 2 4:9]), 'xtick');
            linkprop(hax([2 1]), 'ytick');
            linkprop(hax([6 5 4]), 'ytick');
            linkprop(hax([9 8 7]), 'ytick');

            hax(1).YTick = hax(1).YTick(2:end);
            hax(4).YTick = hax(4).YTick(2:end);
            hax(7).YTick = hax(7).YTick(2:end);

            hax(1).YColor = [1 1 1]*0.3;
            hax(4).YColor = [1 1 1]*0.3;
            hax(7).YColor = [1 1 1]*0.3;

            hax(7).XColor = [1 1 1]*0.3;
            hax(8).XColor = [1 1 1]*0.3;
            hax(9).XColor = [1 1 1]*0.3;

            hax(8).XLabel.String = labx;
            hax(4).YLabel.String = laby;

            axes(hax(2));
            if ~isempty(factor)
                title(['Integrated to ' num2str(factor) 'H_{sb}']);
            end

            % save to file
            fname = ['./params/param_' str];
            slope = mplt;
            intercept = cplt;
            hash = githash('./@runArray/print_diag.m');
            comment = ['(slope, intercept) = straight line fit | (merr,cerr) = error in fit ' ...
                       '| isobath = location index | factor = integrate down to factor x ' ...
                       'H_sb | (diags, plotx, err, norm) = function(isobath, run)'];
            save(fname, 'isobath', 'slope', 'intercept', 'merr', 'cerr', 'hash', 'factor', ...
                 'diags', 'plotx', 'err', 'norm', 'color');

            mmagn = 10^(orderofmagn(mplt(1)));
            cmagn = 10^(orderofmagn(cplt(1)));

            %minmagn = max([mmagn cmagn]);
            cmagn = 1/(cmagn)^2;
            mmagn = 1;

            hax(slopeplot) = subplot(3,3,slopeplot);
            errorbar(runArray.array(1).csflux.ndloc(isobath), ...
                     mplt./mmagn, merr./mmagn, 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
            hold on;
            errorbar(runArray.array(1).csflux.ndloc(isobath), ...
                     cplt/cmagn, cerr/cmagn, '.-', ...
                     'Color', [1 1 1]*0.55, 'LineWidth', 2, 'MarkerSize', 20);
            if mmagn == 1
                mlegstr = 'Slope';
            else
                mlegstr = [num2str(1/mmagn) ' x Slope'];
            end
            if cmagn == 1
                clegstr = 'c';
            else
                clegstr = ['y-intercept x ' num2str(1/cmagn)];
            end
            legend(mlegstr, clegstr, 'Location', 'NorthEast', 'FontSize', 18);
            xlabel('Location (y/R)');
            xlim([-0.05 2]);
            hax(slopeplot).YTick = sort([hax(slopeplot).YTick ...
                                min(mplt)*mmagn max(mplt)*mmagn]);
            correct_ticks('y', '%.2f', []);
            beautify([18 20 22]); pbaspect([1.732 1 1]);
            liney(0);
            set(gcf, 'renderer', 'opengl');
        end

        function [] = plot_field(runArray, varname, tind)

            n = length(runArray.filter);

            if ~exist('tind', 'var'), tind = []; end

            if length(tind) == 1 && ~ischar(tind)
                tind = tind * ones([1 n]);
            end

            figure;
            a = ceil(n/2);
            for ff = 1:length(runArray.filter)
                run = runArray.array(runArray.filter(ff));

                if strcmpi(tind, 'max flux')
                    [~,tindex] = run.calc_maxflux;
                end

                if ~ischar(tind)
                    tindex = tind(ff);
                end

                hax = subplot(a,2,ff);
                run.animate_field(varname, hax, tindex, 1);
            end
        end

        function [handles] = plot_avgProfile(runArray, varname, normalize, ...
                                             axname, ix, mask, hax, reference)

            if ~exist('normalize', 'var'), normalize = 0; end
            if ~exist('axname', 'var'), axname = 'y'; end
            if ~exist('ix', 'var'), ix = 'sb'; end
            if ~exist('mask', 'var'), mask = 1; end
            if ~exist('reference', 'var'), reference = 'center'; end

            if isempty(runArray.filter)
                runArray.filter = [1:runArray.len];
            end
            if ~exist('hax', 'var') | isempty(hax)
                figure; hax = gca;
            else
                axes(hax);
            end
            handles.hax = hax;
            hold all
            insertAnnotation('runArray.plot_avgProfile');
            for ii=1:length(runArray.filter)
                ff = runArray.filter(ii);
                run = runArray.array(ff);
                names{ii} = runArray.getname(ff);

                [vm,mm,xi] = run.avgProfile(varname, axname, ix, mask);
                if strcmpi(varname, 'zeta'), vm = vm - min(vm); end
                ind = find(mm == 1, 1, 'first');
                if normalize
                    vm = vm ./ vm(ind);
                end

                if axname == 'y'
                    %hplt(ii) = plot(xi/1000, vm);
                    %plot(xi(ind)/1000, vm(ind), 'x', ...
                    hpltsh(ii) = plot(xi(1:ind)/1000, vm(1:ind));
                    hplted(ii) = plot(xi(ind:end)/1000, vm(ind:end), '--');
                    linkprop([hpltsh(ii) hplted(ii)], 'Color');
                else
                    hpltsh(ii) = plot(vm(1:ind), xi(1:ind)/1000);
                    hplted(ii) = plot(vm(ind:end), xi(ind:end)/1000, '--');
                    linkprop([hpltsh(ii) hplted(ii)], 'Color');
                end

                phys = run.params.phys;
                betash = phys.f0/run.bathy.hsb * run.bathy.sl_shelf;
                Lbeta = sqrt(run.eddy.V(1)/betash);
                %[~,maxind] = max(vm);
                %ind2 = find_approx(xi(maxind) - Lbeta, xi, 1);
                %plot(xi(ind2)/1000, vm(ind2), '.', 'Color', hplt(ii).Color, 'MarkerSize', 18);
            end

            if axname == 'y'
                ylabel(varname);
                if ~strcmpi(reference, 'center')
                    xlabel('X - X_{front} (km)');
                else
                    xlabel('X - X_{eddy} (km)');
                end
            else
                ylabel('Y - Y_{sb} (km)');
                xlabel(varname);
            end
            title(['Mean ' varname ' | ' axname ' = ' ix]);
            %if mask == 1
            %    hax = gca;
                %hax.Title.String = [hax.Title.String ...
                %                    ' | cross = eddy water boundary'];
                %end
            handles.hleg = legend(hplted, names, 'Location', 'NorthWest');
            beautify;

            hl = linex(0);
            if ~strcmpi(varname, 'zeta'), hl(2) = liney(0); end

            % handles.htxt(1) = text(0.1,0.2, 'shelf water', 'units', 'normalized');
            % handles.htxt(2) = text(0.6,0.2, 'eddy/slope water', 'units', 'normalized');

            handles.hpltsh = hpltsh;
            handles.hplted = hplted;
            handles.hl = hl;
            LowerLines;
        end

        function [] = plot_maxflux(runArray)
            figure; hax = gca;
            subplot(211); hold on;
            subplot(212); hold on;
            for ii=1:8
                [actual, param] = runArray.print_diag('max flux', ii, hax, 'no_plots')
                subplot(211); hold on;
                plot(ii, actual, 'x');
                subplot(212); hold on;
                plot(ii, param, 'x');
            end
        end

        function [] = plot_zetacyc(runArray)

            corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            hf = figure; hold all
            insertAnnotation('runArray.plot_zetacyc');
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.name{ii};
                ndtime = run.eddy.t*86400 / run.eddy.turnover;

                [~,~,tind] = run.locate_resistance;

                run.read_zeta;

                zz = nan(size(ndtime));
                for tt=1:size(run.zeta,3)
                    ix = vecfind(run.rgrid.x_rho(1,:), ...
                                 run.eddy.vor.ee(tt) + ...
                                 run.eddy.vor.dia(1));
                    iy = vecfind(run.rgrid.y_rho(:,1), ...
                                 run.eddy.my(tt));
                    zz(tt) = run.zeta(ix,iy,tt);
                end

                hplt(ff) = plot(ndtime(1:tt), zz);
            end

            legend(hplt, names);
            xlabel('Non-dimensional time');
            title('SSH of "wake cyclone"')
            runArray.reset_colors(corder_backup);
        end

        function [] = check_fittraj(runArray)

            %corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            figure; hold all
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = runArray.name{ii};
                ndtime = run.eddy.t*86400 / run.eddy.turnover;

                [~,~,tind] = run.locate_resistance;
                %run.fit_traj();
                %vel = smooth(run.eddy.mvy, 20);

                %plot(ff, vel(run.traj.tind)./max(abs(vel(:))),
                %'x');

                plot(ff, run.eddy.Lgauss(tind)./run.eddy.Lgauss(1), ...
                     'kx');
            end

            %runArray.reset_colors(corder_backup);
        end

        function [] = check_penetration(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            corder_backup = runArray.sorted_colors;

            %figure; hold all;
            %insertAnnotation('runArray.plot_test3');

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.getname(ii);
                %run.fit_traj;
                %tind = run.traj.tind;
                [~,~,tind] = run.locate_resistance(10,1/4);
                tvec = run.time/run.eddy.turnover;

                run.animate_field('eddye', [], tind, 1);
                title(names{ff});

                %plot(tvec, run.eddy.Lgauss);
                %plot(tvec(tind), run.eddy.Lgauss(tind), 'kx');
            end
            %legend(hplt, names);

            runArray.reset_colors(corder_backup);
        end

        function [] = plot_dhdt(runArray)

            corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            hf = figure; hold all
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.name{ii};
                ndtime = run.eddy.t*86400 / run.eddy.turnover;

                tsl = find_approx(run.eddy.my, run.bathy.xsl, 1);
                tse = find_approx(run.eddy.my - run.eddy.vor.dia(1)/2, ...
                                  run.bathy.xsl, 1);
                [~,~,tind] = run.locate_resistance;

                vec = run.eddy.Lgauss; %./run.eddy.hcen;
                hplt(ff) = plot(ndtime, vec);
                plot(ndtime(tind), vec(tind), 'kx');
                plot(ndtime(tsl), vec(tsl), '.', 'MarkerSize', 34, ...
                     'Color', hplt(ff).Color);
                plot(ndtime(tse), vec(tse), 'o', 'MarkerSize', 14, ...
                     'Color', hplt(ff).Color);
            end

            legend(hplt, names);
            beautify;
            runArray.reset_colors(corder_backup);
        end

        function [] = plot_section_loc(runArray, varname, loc)
        % given a location in (km) plot sections
            fold = runArray.filter;
            for ff=1:length(fold)
                run = runArray.array(fold(ff));

                tind = find_approx(run.eddy.my, loc);

                runArray.filter = fold(ff);
                hh(ff) = runArray.plot_sections(varname, tind);
                hh(ff).LevelList = hh(1).LevelList;
                hf(ff) = gcf;
            end

            linkfig(hf, 'xy');
            runArray.filter = fold;
        end

        function [handles] = plot_ts(runArray, tsname, ax, tname)
        % plot time series
            corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            if ~exist('ax', 'var') | isempty(ax)
                hf = figure; hold all;
                ax = gca();
            else
                axes(ax); hold all;
            end

            if ~exist('tname', 'var')
                tname = [];
            end

            insertAnnotation(['runArray.plot_ts(' tsname ')']);
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                names{ff} = runArray.name{ii};

                use run.params;
                bathy = run.bathy;
                eddy = run.eddy;

                try
                    [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
                catch ME
                    disp(ME);
                end

                eval(['vec = ' tsname ';']);

                if isempty(tname)
                    ndtime = run.eddy.t*86400 / run.eddy.turnover;
                else
                    eval(['ndtime = ' tname ';']);
                end

                try
                    tind = run.FitCenterVelocity; tind = tind(1);
                    if isfield(run.csflux, 'off')
                        [~,maxloc] = run.calc_maxflux(1);
                    end

                    ax.ColorOrderIndex = ff-1;
                    handles.hplt(ff) = plot(ndtime, vec(1:run.eddy.tend));
                    handles.htind(ff) = plot(ndtime(tind), vec(tind), 'kx', ...
                                             'MarkerSize', 16, ...
                                             'HandleVisibility', 'off');
                    if exist('maxloc', 'var')
                        handles.hmaxloc(ff) = plot(ndtime(maxloc), ...
                                                   vec(maxloc), ...
                                                   'ko', ...
                                                   'HandleVisibility', 'off');
                    end
                catch ME
                    disp(ME);
                    try
                        handles.hplt(ff) = plot(ndtime, vec);
                    catch ME
                        if ~exist('handles', 'var') | length(handles.hplt) < ff
                            handles.hplt(ff) = plot(vec);
                        end
                    end
                end
                try
                    handles.hplt(ff).DisplayName = names{ff};
                catch ME
                end
            end

            legend(handles.hplt, names);
            title(tsname);
            xlabel('Time / Turnover time');
            beautify;
            runArray.reset_colors(corder_backup);
        end

        function [handles] = plot_dEdt(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            corder_backup = runArray.sorted_colors;

            figure; subplot(121); hold all; subplot(122); hold all;
            insertAnnotation('runArray.plot_dEdt');

            kk=1;
            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                if ~isfield(run.eddy, 'KE'), continue; end
                %run.fit_traj();
                run.traj = [];

                loc = 1;

                names{kk} = runArray.getname(ii);
                [~,~,tind] = run.locate_resistance();
                tvec = run.time/run.eddy.turnover;

                subplot(121)
                vec = run.eddy.KE(:,loc);
                handles.hke(kk) = plot(tvec, vec./vec(1)); %, 'Color', [1 1 1]*0.5);
                plot(tvec(tind), vec(tind)./vec(1), 'kx', 'MarkerSize', 18);

                subplot(122)
                vec = run.eddy.PE(:,loc);
                vec = vec./vec(1);
                handles.hpe(kk) = plot(tvec, vec./max(vec)); %, 'Color', [1 1 1]*0.5);
                plot(tvec(tind), vec(tind)./max(vec), 'kx', 'MarkerSize', 18);

                kk = kk+1;
            end
            legend(handles.hke, names);
            subplot(121);
            ylim([0 1.1]);
            ylabel('KE / KE_0');
            title(['Crosses at traj.tind. Values normalized by initial ' ...
                   'value']);
            xlabel('Time / Turnover time');
            beautify;

            subplot(122);
            ylim([0 1.1]);
            ylabel('PE / PE_0');
            xlabel('Time / Turnover time');
            beautify;

            %packrows;
            runArray.reset_colors(corder_backup);
        end

        function [] = plot_spectra(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            figure; hold all

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = runArray.getname(ii);


                % find sponge edges
                sz = size(run.sponge);
                sx1 = find(run.sponge(1:sz(1)/2,sz(2)/2) == 0, 1, 'first');
                sx2 = sz(1)/2 + find(run.sponge(sz(1)/2:end,sz(2)/2) == 1, 1, ...
                                     'first') - 2;
                sy2 = find(run.sponge(sz(1)/2, :) == 1, 1, 'first') - 1;

                % indices to look at
                % west, then east
                ix = [sx1 + 10; sx2 - 10];
                iy = [sy2 - 10; sy2 - 10];

                for nn=1:length(ix)
                    u = dc_roms_read_data(run.dir, 'u', [], {'x' ix(nn) ix(nn); ...
                                        'y' iy(nn) iy(nn); 'z' 72 72}, ...
                                          [], run.rgrid, 'his', 'single');

                    v = dc_roms_read_data(run.dir, 'v', [], {'x' ix(nn) ix(nn); ...
                                        'y' iy(nn) iy(nn); 'z' 72 72}, ...
                                          [], run.rgrid, 'his', 'single');

                    [t,iu,~] = unique(run.time, 'stable');

                    inertialfreq = run.params.phys.f0/2*pi;
                    [psi, lambda] = sleptap(length(iu));
                    [F,S] = mspec(1/2*(u(:,iu).^2 + v(:,iu).^2)', ...
                                  psi);
                    % fourier frequencies so that Nyquist is at 1
                    f = fourier(86400, length(iu))/2*pi;

                    hgplt = loglog(f./inertialfreq, S(:,end));
                    set(gca, 'XScale', 'log');
                    set(gca, 'YScale', 'log');
                    addlegend(hgplt, [name ' | ' num2str(ix(nn))]);
                end
            end
            linex(1);
        end

        function [] = plot_envelope(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            figure; hold all

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = runArray.getname( ii);

                env = run.csflux.off.shelfwater.envelope;
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

        function [] = plot_enflux(runArray)

            corder_backup = runArray.sorted_colors;

            figure;
            hold all

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                ndtime = run.asflux.time / run.eddy.turnover;

                filter = 'topo.';
                eval(['teflux = run.asflux.' filter 'ikeflux(:,3) + ' ...
                      'run.asflux.' filter 'ipeflux(:,3)' ...
                      '- run.asflux.' filter 'ikeflux(:,2) - ' ...
                      'run.asflux.' filter 'ipeflux(:,2);']);
                hgplt = plot(ndtime, teflux);
                addlegend(hgplt, run.name);
            end
            liney(0);
            runArray.reset_colors(corder_backup);
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
                vec2 = run.csflux.off.gshelf(run.tscaleind:end);

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
