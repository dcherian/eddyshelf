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

        % sort members of the array by runArray.sort_param;
        function [] = sort(runArray, sort_input)

            if ~exist('sort_input', 'var') || isempty(sort_input)
                error('need sort_input to sort!');
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
            corder_backup = get(0, 'DefaultAxesColorOrder');
            if runArray.sorted
                if isempty(runArray.filter)
                    len = runArray.len;
                else
                    len = length(runArray.filter);
                end

                set(0, 'DefaultAxesLineStyleorder','-');
                set(0, 'DefaultAxesColorOrder', brighten(cbrewer('seq','Reds',len), ...
                                                         -0.5));
            end
        end

        function [] = reset_colors(runArray, corder_backup)
            if runArray.sorted
                set(0, 'DefaultAxesColorOrder', corder_backup);
                set(0,'DefaultAxesLineStyleOrder',{'-','--','-.'});
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
                run = runArray.array(ii);
                out = eval(['run.' command]);
                if ~ischar(out)
                    out = num2str(out);
                end
                disp([runArray.array(ii).name ' | ' out]);
            end
        end

        function [diags, plotx] = print_diag(runArray, name)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            plots = 1;

            diags = nan(size(runArray.filter));

            if plots
                hfig = figure;
                hold all;
                name_points = 1; % name points by default
                labx = ' '; laby = ' ';
                plotx = [];
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                runName = runArray.getname(ii);

                % some commonly used variables
                tind = run.tscaleind;
                ndtime = run.eddy.t * 86400 / run.eddy.turnover;
                Lx = run.eddy.vor.lmaj;
                Ly = run.eddy.vor.lmin;
                Lz = run.eddy.Lgauss;
                Ro = run.eddy.Ro;
                V = run.eddy.V;
                if isfield(run.eddy, 'Vb'), Vb = run.eddy.Vb; end
                hcen = run.eddy.hcen;

                A = run.eddy.amp(1);

                Lr = run.rrdeep;
                beta = run.params.phys.beta;
                use run.params.phys

                alpha = run.bathy.sl_slope;
                hsb = run.bathy.hsb;
                hsl = run.bathy.hsl;
                xsb = run.bathy.xsb;
                xsl = run.bathy.xsl;
                Sa = run.bathy.S_sl;
                diagstr = [];

                %%%%% dummy
                if strcmpi(name, 'dummy')

                    % for local plots
                    %figure; hold all
                    %cb = runArray.sorted_colors;
                    sortedflag = 0;

                    diags(ff) = [];
                    % x-axis variable for plots
                    plotx(ff) = [];
                    % label points with run-name?
                    name_points = 0;
                    % x,y axes labels
                    labx = [];
                    laby = [];
                end

                %%%%% cross-isobath translation velocity
                if strcmpi(name, 'mvy')

                    local_plot = 0;

                    if local_plot
                        % for local plots
                        if ff == 1
                            hfig2 = figure; hold all
                            cb = runArray.sorted_colors;
                            sortedflag = 1;
                        end
                    end

                    Ro = Ro(1);

                    % convert to m/s
                    mvy = smooth(run.eddy.mvy, 28) * 1000/86400;
                    [diags(ff), ind] = min(mvy);
                    diags(ff) = abs(diags(ff));

                    if local_plot
                        ndtime = run.eddy.t * 86400 ./ ...
                                 run.eddy.turnover;
                        figure(hfig2);
                        hgplt = plot(ndtime, mvy);
                        addlegend(hgplt, runName);
                        plot(ndtime(ind), mvy(ind), 'k*');
                    end

                    % make normal summary plot?
                    plots = 1;
                    % x-axis variable for plots
                    plotx(ff) = 5*(1+Sa)^(-1).*(beta.*Lr^2).^2 .* f0*Lr./9.81./A;
                    % label points with run-name?
                    name_points = 1;
                    % x,y axes labels
                    labx = 'Parameterization';
                    laby = 'max |cross-isobath velocity|';
                end

                %%%%% slope parameter
                if strcmpi(name, 'slope param')
                    plots = 0;
                    diags(ff) = Ro(1)./Sa;
                end

                %%%% resistance
                if strcmpi(name, 'resistance')
                    [xx,yy,tind] = run.locate_resistance;

                    loc = 'cen';
                    loc = 'edge';

                    % save for later
                    run.eddy.res.xx = xx;
                    run.eddy.res.yy = yy;
                    run.eddy.res.tind = tind;
                    run.eddy.res.comment = ['(xx,yy) = center ' ...
                                        'location | tind = time index'];

                    plotx(ff) = 1./(1+Ro(1));
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

                    diags(ff) = hdiag ./ Lz(1);

                    % name points with run names on summary plot
                    name_points = 1;
                    labx = 'Ro/S_\alpha';
                end

                %%%%% energy loss
                if strcmpi(name, 'energy loss')

                    % algorithm:
                    %  1. Detect the minimum in vertical scale.
                    %  2. approximate dh/dt as Δh/Δt
                    %  3. compare with cg estimate
                    %
                    % Non-dimensionalizations:
                    %  1. time → eddy turnover time
                    %  2. height → initial vertical scale of eddy

                    time =  run.eddy.t * 86400;
                    ndtime = time ./ (run.eddy.vor.lmaj(1)./run.eddy.V(1));

                    vec = smooth(run.eddy.energy.intTE, 30);
                    [xmax, imax, xmin, imin] = extrema(vec);

                    % find the first minima in the height time series
                    % presumably, this gives me a good slope
                    imins = sort(imin, 'ascend');
                    index = imins(find(imins > 30, 1, 'first'))

                    try
                        intTE = run.eddy.energy.intTE;
                    catch ME
                        intTE = run.eddy_energy_ideal;
                    end

                    fig = 0;
                    if fig
                        field = intTE; run.eddy.Lgauss;
                        figure; plot(ndtime, field./field(1)); hold all
                        plot(ndtime, vec./vec(1));
                        linex(ndtime(index));
                        %plot(ndtime, vec);
                        %plot(ndtime(index), field(index), 'x', ...
                        %     'MarkerSize', 12);
                        title(run.name);
                    end

                    dt = ndtime(index) - ndtime(1);
                    dhdt = (run.eddy.Lgauss(1) - run.eddy.Lgauss(index))./dt;

                    dEdt = (intTE(1) - intTE(index))./intTE(1)./dt;
                    diags(ff) = dEdt; dhdt;

                    plotx(ff) = run.topowaves; labx = 'cg (m/s)';
                    %plotx(ff) = run.eddy.Ro(1)./run.bathy.S_sl; labx = 'Ro/S_\alpha';
                    laby = 'dE/dt';
                end

                %%%%% topography parameter - taylor column
                if strcmpi(name, 'hogg')
                    diags(ff) = 1./Ro(1) ./ ...
                           (1+hsb/Lz(1) * 1/Sa);
                end

                %%%%% rest latitude hypothesis?
                if strcmpi(name, 'rest')

                end

                %%%%% estimate slope for bottom torque balance
                if strcmpi(name, 'slope est')
                    slope = syms_angmom(run);

                    sltind = find_approx(2*pi*slope, alpha);

                    diags(ff) = run.traj.H;
                    plotx(ff) = run.eddy.hcen(sltind);
                end

                %%%%% Flierl (1987) bottom torque hypothesis.
                %% estimate ∫∫ψ
                if strcmpi(name, 'btrq est')
                    %plot(run.eddy.btrq(:,1)/1025);

                    figure;
                    subplot(211); hold all
                    hplt = plot(ndtime, beta .* Lz ./f0 .* (V/Vb)*1./alpha);
                    addlegend(hplt, runName);
                    %try
                    %    plot(ndtime, Vb./V /6);
                    %catch ME
                    %    plot(ndtime, Vb./V'/6);
                    %end

                    title(runName);
                    % rhoamp = rho0 * TCOEF * run.eddy.T(:,end)';
                    % %plot((run.eddy.mass-1000*run.eddy.vol) .* g .* alpha./rho0)
                    % %uvec = f0 .* bsxfun(@times, run.eddy.vol, V');
                    % %plot(uvec);

                    %hold all
                    %plot(ndtime, g*alpha*rhoamp/rho0);
                    %plot(ndtime, beta.*run.eddy.Ls.*V);

                    subplot(212); hold all;
                    plot(ndtime, run.eddy.hcen);
                    %subplot(211);
                    %plot(beta .* run.eddy.voltrans)
                    %hold all
                    %subplot(212);
                    %plot(run.eddy.btrq/1025); hold all
                    %legend('Transport 1', 'Transport 2', ...
                    %       'bottom torque 1', 'bottom torque 2');
                    continue;
                end
                %%%%% Flierl (1987) bottom torque hypothesis.
                if strcmpi(name, 'bottom torque')

                    Ro = Ro(1);
                    Lx = run.eddy.vor.dia(1)/2;
                    H0 = hcen(1);
                    Tamp = run.params.eddy.tamp;
                    TCOEF = run.params.phys.TCOEF;

                    tanhfitflag = 0;
                    if tanhfitflag
                        run.fit_traj;

                        H = run.traj.H;
                        Y = run.traj.Y;
                        tind = run.traj.tind;
                    else
                        [~,~,tind] = run.locate_resistance;
                        H = (run.eddy.hcen(tind));
                        Y = run.eddy.my(tind) - run.bathy.xsb; run.eddy.my(tind);
                    end
                    if isempty(tind), continue; end

                    %if H./run.bathy.hsl > 0.9
                    %    H = (run.eddy.hedge(tind));
                    %    disp(runName);
                    %end

                    %figure;
                    %hplt = plot(ndtime, run.eddy.hcen); hold on
                    %addlegend(hplt, runName);
                    %plot(ndtime(tind), run.eddy.hcen(tind), 'k*');
                    %%liney(yscl); linex([1 2 3]*tscl);

                    diags(ff) = H./Lz(1);
                    %plotx(ff) = (beta*Lx)/alpha * V(1)/(Tamp*TCOEF);
                    %diags(ff) = Y./(run.rrdeep);
                    %plotx(ff) = Ro./Sa;
                    plotx(ff) = beta*Lz(1)/alpha/f0;
                    name_points = 1;
                    laby = 'H_{cen}/L_z^0';
                    labx = '\beta L_z^0 / (\alpha f_0)';

                    % there is an older version that didn't really work
                    %  ndtime = run.eddy.t*86400 ./ run.eddy.turnover;
                    % Lx = sqrt(run.eddy.vor.lmaj .* run.eddy.vor.lmin);
                    % c = smooth(run.eddy.mvx, 10)';

                    % Vb = run.eddy.Vb;

                    % hcen = run.eddy.hcen';

                    % % last 'n' turnover periods
                    % n = 30;
                    % ind = find_approx(run.eddy.t, run.eddy.t(end) - ...
                    %                   n * run.eddy.turnover/86400);

                    % % average quantities
                    % Vm = nanmean(V(ind:end));
                    % Vbm = nanmean(Vb(ind:end));
                    % cm = nanmean(c(ind:end));
                    % Lxm = nanmean(Lx(ind:end));
                    % Am = nanmean(run.eddy.vor.amp(ind:end));

                    % h = run.bathy.h(1,:);
                    % seind = vecfind(run.eddy.yr(1,:), run.eddy.vor.se);
                    % hedge = h(seind);

                    % %figure;
                    % %plot(run.eddy.mx, run.eddy.my);
                    % %hold on;
                    % %plot(run.eddy.mx(ind), run.eddy.my(ind), 'k*');
                    % %title(run.name);

                    % %c = runArray.sorted_colors;

                    % %num = alpha * Lx;
                    % %deno = c./Vb +  V./Vb;

                    % d1 = alpha .* f0./beta .* Vbm./Vm;
                    % d2 = cm./Vm .* f0./beta./Lxm .* Am;
                    % diags(ff) = d1/100;

                    % plotx(ff) = mean((hedge(ind:end) + hcen(ind:end)) /2);
                    % %hold all;
                    % %hplt = plot(ndtime, vec);
                    % %addlegend(hplt, run.name);
                    % %plot(ndtime, run.eddy.Lgauss, 'Color', get(hplt, ...
                    % %                                            'Color'))
                    % %runArray.reset_colors(c);
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
                    hfinal = mean(hcen(tind:end));
                    hinit = hcen(1);

                    diag_h = (hinit - hfinal)./Lz(tind);

                    diag_l = (mean(run.eddy.my(tind:end) - ...
                                   xsb))./(run.eddy.vor.dia(tind)/2);

                    diagstr = ['h = ' num2str(diag_h) ...
                               ' | L = ' num2str(diag_l)];
                end

                %%%%% beta v/s beta_t
                if strcmpi(name, 'betas')
                    plots = 0;
                    betat = Sa .* f0 ./ max(run.bathy.h(:));
                    diagstr = [num2str( ...
                        betabeta ./ betat  ...
                        ) ' | ' num2str(mean(hcen(run.tscaleind:end)))];
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
                    transscl = 0.075 * 9.81/f0 .* run.eddy.amp(ind).* hsb/1000;

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
                    plot(plotx(ff), diags(ff), '*');

                    % add function call as annotation
                    insertAnnotation(['runArray.print_diag(' name ')']);
                    % add run names
                    if name_points
                        text(plotx(ff), diags(ff), runName, 'FontSize', ...
                             12, 'Rotation', 90);
                    end
                    xlabel(labx);
                    ylabel(laby);
                    title(name);
                end

                disp([run.name ' | ' name ' = ' diagstr ' | plotx ' ...
                                    '= ' num2str(plotx(ff))]);
            end

            if exist('sortedflag', 'var') & sortedflag
                runArray.reset_colors(cb);
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

                Ro = Ro .* run.eddy.vormask;

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

            figure;
            hold all;

            %corder_backup = runArray.sorted_colors;

            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);

                ndtime = run.eddy.t*86400 / run.eddy.turnover;

                tind = run.traj.tind;
                slope = syms_angmom(run);
                alpha = run.bathy.sl_slope;

                hplt = plot(ndtime, slope);
                addlegend(hplt, run.name);
                ind = find_approx(slope, alpha);
                plot(ndtime(ind), slope(ind), 'k*');
                plot(ndtime(tind), slope(tind), 'ko');
                %addlegend(hgplt, run.name);
                %plot(ndtime(run.tscaleind), run.eddy.vol(run.tscaleind), 'k*');
            end

            %runArray.reset_colors(corder_backup);
        end

        function [] = plot_test3(runArray)
            if isempty(runArray.filter)
                runArray.filter = 1:runArray.len;
            end

            figure; hold all

            for ff=1:length(runArray.filter)
                ii = runArray.filter(ff);
                run = runArray.array(ii);
                name = runArray.getname(ii);

                hplt = plot(run.ndtime, ...
                            (run.eddy.Ro));

                addlegend(hplt, name);
            end
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
