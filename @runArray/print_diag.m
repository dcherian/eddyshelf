% [diags, plotx, err, norm, color, rmse, P, Perr] = ...
%        print_diag(runArray, name, args, hax, commands)
% commands: 'no_name_points'; 'no_sloping_shelf'; 'small_points'

function [diags, plotx, err, norm, color, rmse, P, Perr, handles] = ...
        print_diag(runArray, name, args, hax, commands)

    if ~exist('args', 'var'), args = []; end
    if ~exist('hax', 'var'), hax = []; end
    if ~exist('commands', 'var'), commands = ''; end

    if ischar(hax)
        commands = hax;
        hax = [];
    end

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    if ~strcmpi(name, 'params table') && ~strcmpi(name, 'nondim') ...
            && ~strcmpi(name, 'beta gyre') ...
            && ~strcmpi(name, 'bfric') && ~strcmpi(name, 'hbl') ...
            && ~strcmpi(name, 'arrest') && isempty(strfind(commands, 'no_plots'))
        plots = 1;
    else
        plots = 0;
    end

    annostr = ['runArray.print_diag(' name ')'];
    diags = nan(size(runArray.filter));
    plotx = diags;
    err = repmat(diags, [2 1]);
    norm = ones(size(diags));
    P = nan([2 1]); % slope and intercept
    Perr = nan([2 1]); % error bounds on slope and intercept
    save_diags = 0;

    parameterize = 0;
    if plots
        if isempty(hax)
            hfig = figure;
            hax = gca; % default axes
        else
            hfig = gcf;
            axes(hax); hold on;
        end

        % add function call as annotation
        insertAnnotation(['runArray.print_diag(' name ')']);
        hold all;
        name_points = 1; % name points by default
        line_45 = 0; %no 45° line by default
        errorbarflag = 0; % use errorbar() instead of plot()
        logscale = 0; strip_ew = 1;
        force_0intercept = 0;
        kozak = 0; % fancy Kozak scatterplot
        labx = ' '; laby = ' ';
        clr = [0 0 0];
        mark_outliers = 0;
        titlestr = name;
        if ~isempty(args), titlestr = [titlestr ' | args = ' num2str(args)]; end
        rmse = 0;
    else
        titlestr = [];
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);
        run = runArray.array(ii);
        runName = runArray.getname(ii);

        ptName = runName;
        clr = [0 0 0]; % reset color
        marker = '.'; % reset marker

        % some commonly used variables
        tind = run.tscaleind;
        [itsl, itse, tsl, tse] = run.getEddyCenterTimeScales;
        [~,uind,~] = unique(run.time, 'stable');
        ndtime = run.eddy.t * 86400 / run.eddy.turnover;
        Lx = run.eddy.vor.lmaj/2;
        Ly = run.eddy.vor.lmin/2;
        Lz = run.eddy.Lgauss;
        Ls = run.eddy.Ls;
        if isfield(run.eddy, 'rhovor')
            Ro = run.eddy.rhovor.Ro;
        end
        V = run.eddy.V;
        if isfield(run.eddy, 'Vb'), Vb = run.eddy.Vb; end

        hcen = run.eddy.hcen;
        fcen = run.eddy.fcen;
        my = run.eddy.my;
        mx = run.eddy.mx;

        A = run.eddy.amp(1);

        Lr = run.rrdeep;
        beta = run.params.phys.beta;
        use run.params.phys
        f0 = abs(f0);

        sgn = sign(run.params.phys.f0) * sign(run.params.eddy.tamp);

        alpha = run.bathy.sl_slope;
        hsb = run.bathy.hsb;
        hsl = run.bathy.hsl;
        xsb = run.bathy.xsb;
        xsl = run.bathy.xsl;
        bathy = run.bathy;
        eddy = run.eddy;

        if run.bathy.axis == 'y'
            fsb = run.rgrid.f(run.bathy.isb,1);
        else
            fsb = run.rgrid.f(1,run.bathy.isb);
        end
        Sa = run.bathy.S_sl;
        Ssh = run.bathy.S_sh;
        beta_t = f0 * alpha./Lz(1);
        N = sqrt(run.params.phys.N2);
        ash = run.bathy.sl_shelf;
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

        %%%%% shelfbc
        if strcmpi(name, 'shelfbc')
            if isempty(args)
                args(1) = 1;
            end

            % for local plots
            %figure; hold all
            %cb = runArray.sorted_colors;
            sortedflag = 0;

            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;
            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

            phi = hsb./(V0./bathy.S_sh/N);

            [clr, ptName, marker] = colorize(run, ptName);
            if run.params.misc.rdrg == 0
                clr = [0 0 0];
            end
            diags(ff) = nanmedian(run.shelfbc.shelf(:,args(1)));
            % x-axis variable for plots
            plotx(ff) = phi;
            % label points with run-name?
            name_points = 1;
            strip_ew = 1;
            % x,y axes labels
            labx = '$$\phi_o$$';
            laby = ['BC_{' num2str(run.shelfbc.thresh(args)) '}'];
        end


        %%%%% sbreakbc
        if strcmpi(name, 'sbreakbc')

            if isempty(args)
                args(1) = 1;
            end

            % for local plots
            %figure; hold all
            %cb = runArray.sorted_colors;
            sortedflag = 0;

            if ~isfield(run.shelfbc, 'sbreak')
                continue;
            end

            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;
            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

            phi = hsb./(V0./bathy.S_sh/N);

            [clr, ptName, marker] = colorize(run, ptName);
            if run.params.misc.rdrg == 0
                clr = [0 0 0];
            end
            diags(ff) = nanmedian(run.shelfbc.sbreak.shelf(:,args(1)));
            % x-axis variable for plots
            plotx(ff) = phi;
            % label points with run-name?
            name_points = 1;
            % x,y axes labels
            labx = '$$\phi_o$$';
            laby = ['BC_{' num2str(run.shelfbc.thresh(args)) '}'];
        end

        % shelfbreak SSH decay scale
        if strcmpi(name, 'sbssh')

            if isinf(bathy.Lbetash), continue; end
            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;

            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);
            [supply, errsupp, eddyonshelf] = run.SupplyJetEddyonShelf;

            diags(ff) = 2*abs(run.sbssh.X/1000);
            plotx(ff) = supply/1000;

            conf = confint(run.sbssh.fitobj);
            err(1,ff) = abs(abs(2*conf(1,2)/1000)-diags(ff));

            parameterize = 1;
            errorbarflag = 1;
            name_points = 1;

            labx = 'Eddy radius (km)';
            laby = {'Shelf-water SSH along-shelf'; 'decay scale at shelfbreak'};
        end

        %%%%% std(offshore flux profile)

        if strcmpi(name, 'fluxprofilestd')

            profile = run.FluxVertProfile(1);

            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;
            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

            phi = hsb./(V0./bathy.S_sh/N);

            diags(ff) = std(profile);
            plotx(ff) = phi;

        end

        %%%%% bc(offshore flux profile)

        if strcmpi(name, 'fluxprofilebc')

            profile = run.FluxVertProfile(1);

            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;
            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

            phi = hsb./(V0./bathy.S_sh/N);

            [~, imax] = max(profile);
            [~, imin] = min(profile);

            diags(ff) = abs(profile(imax)-profile(imin))/abs(profile(imax));
            plotx(ff) = phi;

        end

        %%%%% rhines scale
        if strcmpi(name, 'supply') | strcmpi(name, 'eddyonshelf')

            if isempty(run.supply), continue; end
            if run.params.misc.rdrg ~= 0, continue; end

            sortedflag = 0;

            % averaging between [start,stop] works better because of Rh = 3 runs.
            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;

            Ro = mean(Ro(t0:tend));
            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);
            V0 = eddy.V(1);

            betash = fsb/hsb * bathy.sl_shelf;

            Lbetash = sqrt(V0/(betash-beta));
            Lctw = V0/bathy.sl_shelf/N;
            Ldef = N*hsb/fsb;

            phi = hsb./(V0./bathy.S_sh/N);
            chi = (2/sqrt(pi)*exp(-(hsb/Lz0)^2) * V0/Lz0)/(bathy.S_sl*N);
            lambda = hsb/Lz0;

            % if strcmpi(run.name, 'ew-8392') | strcmpi(run.name, 'ew-8151')
            %     i0 = 10;
            % else
            %     i0 = 1;
            % end
            % vavg = run.supply.shsl.vavg;
            % [vmin, imin] = min(vavg(i0:end));
            % imin = imin + i0-1;
            % ind = find_approx(vavg(i0:imin), 0.2 * vmin, 1) + i0-1;
            % figure;
            % plot(run.supply.shelf.vavg);
            % linex([imin ind]);
            % title(run.name)

            [supply, errsupp, eddyonshelf] = run.SupplyJetEddyonShelf;

            % if strfind(run.name, '8234');
            %     supply = supply - eddyonshelf;
            % end

            % sqrt(U_mean/v_bot) factor for rhines scale
            % zvec = -hsb:0.1:0;
            % vel = 1 - erf((abs(zvec)./D));
            % vfactor = sqrt(trapz(zvec, vel)./hsb/vel(1));

            if strcmpi(name, 'supply')
                if lambda >= 0.35
                    continue;
                end
                diags(ff) = supply;
                err(1,ff) = errsupp;

                % if strcmpi(run.name, 'ew-8342-2')
                %     exceptions = ff;
                % end

                errorbarflag = 0;
                %diags(ff) = max(smooth(env,10))/1000;
                %diags(ff) = abs(ind-imin); % in km
                laby = {'Width of shelf water';  'supply jet (km)'};

                plotx(ff) = Lbetash/1000;
                labx = 'Rhines scale, L_\beta (km)';
            else % eddy water on shelf
                diags(ff) = eddyonshelf;
                laby = {'Penetration of';  'eddy water on shelf (km)'};
                plotx(ff) = Ldef/1000; %plotx(ff) / sqrt(1-Ro);
                labx = 'Shelf Rossby radius (km)';
            end

            ptName = [num2str(run.bathy.S_sh) ' | ' run.name(4:end)];
            strip_ew = 0;
            name_points = 1;
            parameterize = 1;
            clr = 'k';
            %line_45 = 1;
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

        %%%%% cross-shelf translation time scale
        if strcmpi(name, 'timescale')
            run.fit_traj;

            vel = smooth(run.eddy.cvy, 20);

            % stalling isobath
            h0 = Lz(1) * erfinv(1 - beta/beta_t);

            t0 = run.eddy.tscaleind;

            vel = smooth(run.eddy.mvy, 20);

            %figure;
            %plot(vel);hold all
            %plot(smooth(vel, 20));
            diags(ff) = run.traj.tscl;
            plotx(ff) = beta/f0 .* Lz(t0).^2 ./ alpha^2 .* ...
                exp( (hcen(1)/Lz(1)) ^2) ./ max(-1*vel(1:100));
            %name_points = 0;
            %diags(ff) = -1 * min(vel)./run.eddy.V(1);
            %plotx(ff) = run.params.nondim.eddy.Rh;
        end

        %%%%% slope parameter
        if strcmpi(name, 'slope param')
            plots = 0;
            diags(ff) = beta/beta_t;
            plotx(ff) = NaN;
        end

        %%%%% energy loss - as vertical scale
        if strcmpi(name, 'dhdt')
            use_traj = 0;
            if use_traj
                run.fit_traj(1.1);
                tind = run.traj.tind;
            else
                [~,~,tind] = run.locate_resistance;
            end

            diags(ff) = (Lz(tind)-Lz(1))./Lz(1);
            plotx(ff) = V(1)./beta./Ls(1).^2;

            if plotx(ff) > 25
                plotx(ff) = NaN;
                diags(ff) = NaN;
            end

            laby = '\Delta L_z / L_z^0';
            labx = 'U/\beta L_s^2';
            %plotx(ff) = beta*Ls(1)/f0;% - much worse
        end

        % can run cross shelfbreak?
        if strcmpi(name, 'cross sb')
            tind = find_approx((run.eddy.vor.se - xsb), 0);

            % figure;
            % plot(Lz); linex(tind); liney(hsb);
            % title(run.name);
            % continue;

            diags(ff) = Lz(tind)/hsb;
            plotx(ff) = V(1)./beta./Ls(1).^2;

            laby = '\Delta L_z / L_z^0';
            labx = 'U/\beta L_s^2';
        end


        %%%%% energy loss
        if strcmpi(name, 'dEdt')
            [~,~,tind] = run.locate_resistance;
            try
                PE = abs(run.eddy.PE(1:tind,1));
            catch ME
                run.eddy_bulkproperties;
            end
            KE = run.eddy.KE(1:tind,1);

            tvec = ndtime(1:tind)';
            %dPEdt = nanmean(fillnan(smooth(diff(PE./PE(1))./ ...
            %                    diff(tvec), 30), Inf));
            %dKEdt = nanmean(fillnan(smooth(diff(KE./KE(1))./ ...
            %                    diff(tvec), 30), Inf));


            [~,pind] = min(PE(5:end));
            [~,kind] = min(KE(5:end));

            use_polyfit = 0;
            if use_polyfit
                PEfit = polyfit(tvec(1:pind), PE(1:pind)./PE(1), 1);
                KEfit = polyfit(tvec(1:kind), KE(1:kind)./KE(1), 1);
            else
                [PEfit,bint,r,rint,stats] = regress(PE(1:pind)./PE(1), ...
                                                    [tvec(1:pind) ...
                                    ones([pind 1])]);
                error_PE = bint(1,2) - PEfit(1);

                [KEfit,bint,r,rint,stats] = regress(KE(1:kind)./KE(1), ...
                                                    [tvec(1:kind) ...
                                    ones([kind 1])]);
                error_KE = bint(1,2) - KEfit(1);
            end

            dPEdt = (PE(pind))./PE(1) - 1; %PEfit(1);
            dKEdt = (KE(kind))./KE(1) - 1; %KEfit(1);

            en = 'PE';

            %laby = ['|d' en '/dt|'];
            laby = ['\Delta' en '/' en '(0)'];

            %figure; hold all;
            %plot(tvec, PE./PE(1)); plot(tvec, tvec*PEfit(1) + ...
            %                            PEfit(2));
            %title(run.name);
            %plot(tvec, KE./KE(1)); plot(tvec, tvec*KEfit(1) + KEfit(2));

            eval(['diags(ff) = abs(d' en 'dt);']);
            plotx(ff) = beta * Lx(1) ./ f0; %...
            %(beta * Lz(1) - alpha * f0 * (1-erf(hedge(1)./Lz(1))));
            eval(['err(1,ff) = abs(error_' en ');']);

            errorbarflag = 1;
            name_points = 1;
            labx = '\beta L / f_0';
            %labx = ['$$\beta L_z - ' ...
            %        'f_0 \alpha_{bot} (1 - \mathrm{erf} (h_{edge}/L_z) )  $$'];
            titlestr = '';
        end

        if strcmpi(name, 'energy loss old')

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

            try
                if any(size(run.eddy.KE) == 1)
                    TE = run.eddy.KE + run.eddy.PE;
                else
                    TE = run.eddy.KE(:,1) + run.eddy.PE(:,1);
                end
            catch ME
                disp(run.name);
                continue;
            end

            vec = smooth(TE, 30);
            [xmax, imax, xmin, imin] = extrema(vec);

            % find the first minima in the height time series
            % presumably, this gives me a good slope
            imins = sort(imin, 'ascend');
            index = imins(find(imins > 30, 1, 'first'))

            fig = 0;
            if fig
                field = run.eddy.Lgauss;
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

            dEdt = (TE(1) - TE(index))./TE(1)./dt;
            diags(ff) = dEdt; dhdt;

            plotx(ff) = run.topowaves; labx = 'cg (m/s)';
            %plotx(ff) = run.eddy.Ro(1)./run.bathy.S_sl; labx = 'Ro/S_\alpha';
            laby = 'dE/dt';
        end

        %%%%% buoyancy arrest
        if strcmpi(name, 'arrest')
            try
                Q = V(1) * Lx(1) * Lz(1);
                r = run.params.misc.rdrg;
                c5 = 39; % 39 (inflow) or 96 (outflow)
                d = sqrt(N2)/f0 * r / V(1); % C_D N/f

                diags(ff) = c5/d * (1+Sa^2)/Sa^2 * ...
                    (abs(Q)/hcen(tind)/f0)/3000;
            catch ME
                diags(ff) = NaN;
            end
        end

        %%%%% bottom friction
        if strcmpi(name, 'bfric')
            rdrg = run.params.misc.rdrg;
            tscale = 1./(beta.*Lx(1)); % [s]
            ftscale = hcen(1)./rdrg; % [s]
            diags(ff) = tscale./ftscale;
            % this is correct.
            % if tscale ≫  frictional tscale,
            % friction will kill the wave.
        end

        if strcmpi(name, 'hbl')
            vi = run.eddy.V(1);
            s = run.bathy.S_sl;
            N = sqrt(N2);
            Ri = 0.4;

            L = 1/2 * (-1 + sqrt(1 + 4*Ri*s^2));
            diags(ff) = -vi/(s*N) * L;
        end

        %%%%% β gyre timescale
        if strcmpi(name, 'beta gyre')
            diags(ff) = (1./(beta .* Lx(1)))/run.eddy.turnover;
            plotx(ff) = 0;
        end

        %%%%% estimate slope for bottom torque balance
        if strcmpi(name, 'slope est')
            [slope,pbot,angmom] = syms_angmom(run);

            sltind = find_approx(beta*angmom - pbot* alpha, 0);

            figure;
            plot(ndtime, hcen);
            hold on; plot(ndtime(sltind), hcen(sltind), 'k*');
            diags(ff) = run.traj.H;
            plotx(ff) = run.eddy.hcen(sltind);

            laby = 'Actual water depth';
            labx = 'Water depth where angmom/pbot = \alpha';
            name_points = 1;
            line_45 = 1;
        end

        %%%%% Flierl (1987) bottom torque hypothesis.
        if strcmpi(name, 'bottom torque')
            Ro = Ro(1);
            %Lx = run.eddy.vor.dia(1)/2;
            H0 = hcen(1);
            Tamp = run.params.eddy.tamp;
            TCOEF = run.params.phys.TCOEF;

            if run.params.nondim.eddy.Rh(1) > 50, continue; end

            avgresflag = 0; % this is crap. keep it 0
            if avgresflag %|| run.params.eddy.tamp < 0
                nsmooth = 10;
                [~,~,tind,H,erres] = run.averageResistance(nsmooth);
                titlestr = ['Using average resistance'];
                errorbarflag = 1;

                %factor = 1/4;
                %[~,~,tind] = run.locate_resistance(nsmooth, factor);
                %Y = run.eddy.my(tind) - run.bathy.xsb; ...
                %    run.eddy.my(tind);
            else
                [tind,~,~,BadFitFlag] = run.FitCenterVelocity;
                if BadFitFlag
                    warning(['Bad Fit: Skipping ' run.name '.'])
                    continue;
                end

                titlestr = ['Using fits to translation velocity'];
                errorbarflag = 1;
            end

            if isempty(tind) | isnan(tind)
                warning(['Skipping  ' run.name '. locate_resistance did not work.']);
                continue;
            end

            [itsl,itse,tsl,tse] = run.getEddyCenterTimeScales;
            t0 = tind; run.eddy.tscaleind;
            titlestr = [titlestr ' | t0 = ' num2str(t0)];

            %Lz = run.eddy.Lzfit;
            %Lz(Lz > 1e3) = nan;

            trange = 1:itsl; %tind(1); %1:ceil(mean([itsl tind(1)]));
            Lz0 = nanmean(run.eddy.Lgauss(trange)); %Lz(1,tind(1));
            z0  = 0; %nanmean(run.eddy.z0fit(trange,1));
            % sometimes the fits are trash.
            % go back one time instant at a time till something reasonable is found.
            % while abs(z0) > 1e3 | abs(Lz0) > 1e3
            %     tind = tind - 1;
            %     z0 = run.eddy.z0fit(tind(1),1);
            %     Lz0 = run.eddy.Lzfit(1,tind(1));
            %end
            %z0 = 0;
            % [estimate, lower bound, upper bound]
            % Since H decreases with time, using upper bound on tfit
            % gives lower bound on H
            H = run.eddy.hcen([tind(1) tind(3) tind(2)]);

            beta_z = (alpha*abs(f0)/Lz0);

            [clr, ptName, marker] = colorize(run, ptName);
            clr = 'k';

            zz = (H(1) - z0)./Lz0;
            if strcmpi(args, 'regression') | isempty(args)
                diags(ff) = (1 - erf(zz))./(1-erf(-z0/Lz0));
                plotx(ff) = beta ./ (beta_z - beta);

                parameterize = 1;
                laby = '$$\frac{U_b}{U_s} = 1 - \mathrm{erf}(\frac{H}{L_z^0})$$';
                labx = '$$\frac{\beta}{\beta_z - \beta}$$';
            end

            if strcmpi(args, 'functional')
                diags(ff) = zz;
                plotx(ff) = beta/beta_z;

                parameterize = 0;
                laby = '$$\frac{H_{arr}}{L_z}$$';
                labx = '$$\frac{\beta}{\beta_z}$$';
            end

            if isnan(diags(ff)), plotx(ff) = NaN; end

            % if avgresflag
            %     % these are not justifiable error bounds
            %     zzlo = (-erres(4))./Lz0;
            %     zzhi = (+erres(4))./Lz0;

            %     err(1,ff) = (1 - erf(H./Lz0+zzlo)) - diags(ff);
            %     err(2,ff) = (1 - erf(H./Lz0+zzhi)) - diags(ff);

            %     err(1,ff) = exp(-(H./Lz0)) - diags(ff);
            %     err(2,ff) = exp(-(H./Lz0+zzhi)) - diags(ff);
            % else
            %     %err(1,ff) = (1 - erf(H(2)./Lz0)) - diags(ff);
            %     %err(2,ff) = (1 - erf(H(3)./Lz0)) - diags(ff);
            %     % F = 1 - erf( (H-z0)/Lz )
            %     % error bounds are hopeless. z0, Lz are correlated.
            %     % and then numerator and denominator of y-axis are correlated.
            %     F0 = 2/sqrt(pi) * exp( - ((H(1)-z0)./Lz0)^2 );
            %     dFdH = F0./Lz0;
            %     dFdz0 = F0./Lz0;
            %     dFdLz = F0 * (H(1) - z0)./Lz0^2;

            %     dH = mean(abs(diff(H([2 1 3]))));
            %     dLz = mean([abs(Lz(2,tind(1)) - Lz0), abs(Lz(3,tind(1)) - Lz0)]);
            %     dz0 = 0; %mean([abs(z0(tind(1),2) - z0), abs(z0(tind(1),3) - z0)]);

            %     dnum = sqrt( dFdH^2 * dH^2 + dFdLz^2 * dLz^2 + dFdz0^2 * dz0^2 );
            %     dden = 0; sqrt( dFdLz^2 * dLz^2 + dFdz0^2 * dz0^2 );

            %     err(1,ff) = sqrt(dnum^2 + dden^2);
            %     err(2,ff) = sqrt(dnum^2 + dden^2);
            % end

            kozak = 1;
            name_points = 1; line_45 = 0;
            slope = num2str(round(alpha*1e3));
            errorbarflag = 0;

            % mark NS isobath runs with gray points
            if run.bathy.axis == 'x'
                clr = [1 1 1]*0.75;
            end
        end

        if strcmpi(name, 'btrq')
            % β ∫∫∫ yu = αg (∫∫η + 1/ρ0 ∫∫∫ρ' )

            [angmom,pbot] = run.estimate_bottrq;

            diags(ff) = angmom;
            plotx(ff) = pbot;
        end

        if strcmpi(name, 'velprofile')
            it = [1 tind];
            ix = vecfind(run.rgrid.x_rho(1,:), run.eddy.mx(it));
            iy = vecfind(run.rgrid.y_rho(:,1), ...
                         run.eddy.my(it) - ...
                         run.eddy.vor.lmin(it)/3);

            for tt=1:length(it)
                u(:,tt) = ...
                    dc_roms_read_data(run.dir, 'u', it(tt), {'x' ix(tt) ...
                                    ix(tt); 'y' iy(tt) iy(tt)}, [], ...
                                      run.rgrid, 'his');
                [u0(tt),Z1(tt)] = ...
                    erf_fit((double(run.rgrid.z_u(:,iy(tt),ix(tt)))), ...
                            (u(:,tt)./u(end,tt)), 0);

            end

            % figure out factor to account for gradient wind
            Lz = run.eddy.Lgauss(it);
            H = run.eddy.hcen(tind);
            run.eddy.grfactor = Z1./Lz*-1;

            plotx(ff) = beta/beta_t;
            diags(ff) = run.eddy.grfactor(2);
        end

        if strcmpi(name, 'rh')
            [~,~,tind] = run.locate_resistance;

            plotx(ff) = ff;
            diags(ff) = V(tind)./beta/Ls(tind).^2;
        end

        if strcmpi(name, 'deccel')
            [itsl,itse,tsl,tse] = run.getEddyCenterTimeScales;
            [~,~,restind] = run.averageResistance;
            if isnan(restind), continue; end
            tres = run.time(restind);

            [~,tfit,tscl] = run.FitCenterVelocity;
            tfit = tfit(1);
            V0 = run.eddy.V(itse);
            Lrh = sqrt(V0/beta);

            diags(ff) = tscl; %(tres-tse);
            plotx(ff) = 1./beta ./abs(run.rrdeep);

            laby = 't_{res} - t_{SE}';
            name_points = 1;
        end

        if strcmpi(name, 'initial distance')
            [~,tfit,tscl] = run.FitCenterVelocity;

            diags(ff) = tscl(1);
            plotx(ff) = sgn * (run.params.eddy.cy - run.bathy.xsl)/abs(run.rrdeep);
        end

        % compare against flat bottom rest latitude
        if strcmpi(name, 'rest latitude')
            [~,~,tind] = run.locate_resistance();

            Yb2 = max(run.rgrid.y_rho(:))/2;
            y = Yb2 + 1./beta * (fcen(1)*(1-sgn*Ro(tind)) - f0);

            dy = (y - sgn*run.eddy.my(tind));

            diags(ff) = dy/1000;
            plotx(ff) = ff;
        end

        % compare rest latitude for shallow water PV
        if strcmpi(name, 'shallow PV')
            [~,~,tind] = run.locate_resistance();
            diags(ff) = hcen(1) * fcen(tind)./fcen(1) .* ...
                (1 - Ro(tind))./(1 - Ro(1)) ./ hcen(tind);
            plotx(ff) = ff;
        end

        % nondim parameters to compare runs
        if strcmpi(name, 'diff') || strcmpi(name, 'nondim')
            [~,~,tres] = run.locate_resistance;

            [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
            t0 = start;
            tend = stop;
            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

            Ldefsh = N * hsb / fsb;
            betash = f0 * run.bathy.sl_shelf./run.bathy.hsb;
            Lbetash = sqrt(V(1)/(betash-beta));

            phi = hsb./(V0./bathy.S_sh/N);
            chi = (2/sqrt(pi)*exp(-(hsb/Lz0)^2) * V0/Lz0)/(bathy.S_sl*N);

            rdrg = run.params.misc.rdrg;
            Lf = sqrt(rdrg/fsb/bathy.sl_shelf * 50e3);
            Lf = bathy.L_shelf.^2 / (rdrg/fsb/bathy.sl_shelf);

            F = rdrg /2/betash/hsb/1000;
            dL = F * (1+F/2/Lbetash);
            if ff == 1
                close;
                disp(sprintf('| %10s | %5s | %4s | %4s | %8s | %8s | %5s | %6s | %6s | %4s | %5s | %5s |', ...
                             '', 'Rh', 'Ro', 'bl/f', 'beta_sh', ...
                             'L_rh/Lsh', 'phi', 'xi', 'lambda', 'S_sh', 'L_rh', 'L_def/L_rh', ...
                             'Latw'));
                disp(repmat('-', [1 100]));
            end
            disp(sprintf('| %10s | %5.2f | %4.2f | %4.2f | %.2e | %8.2f | %5.2f | %7.2f |  %5.2f | %0.2f | %5.2f | %5.2f | %5.2f |' , ...
                         run.name, run.params.nondim.eddy.Rh, ...
                         Ro(1), beta*Lx(1)/f0, betash, ...
                         Lbetash./run.bathy.L_shelf, ...
                         phi, chi, hsb./Lz(itsl), bathy.S_sh, Lbetash/1000, Ldefsh/Lbetash, Lf/1000));
            continue;
        end

        if strcmpi(name, 'params table')
            [~,~,tres] = run.locate_resistance;

            % cross-shelf flux runs
            if strfind(args, 'csf')
                if ff == 1
                    fid = fopen(args, 'w');
                    fprintf(fid, ['|- \n | Name | $L$ (km) | $L_z$ (m) | $\\Ro$ | $\\lambda$ | $L\\sh$ (km) | ' ...
                                  '$L\\sl$ (km) | $S\\sh$ | $S\\sl$ | $N^2$ | $f_0$ | $\\beta$' ...
                                  '| $r_f$ (m/s) | \n |- \n']);
                    fclose(fid);
                end
                fid = fopen(args, 'a');
                fprintf(fid, ['| %s | %2.0f | %3.0f | %.2f | %.2f | %3d | %3d | %.2f | %.2f ' ...
                              '| %1.1e | %1.1e | %1.1e | %1.1e |\n'], ...
                        runName, Lx(itsl)/1000, Lz(itsl), ...
                        Ro(1), hsb./Lz(itsl), bathy.L_shelf/1000, ...
                        bathy.L_slope/1000, bathy.S_sh, bathy.S_sl, ...
                        N^2, f0, beta, run.params.misc.rdrg);
            end

            % sl runs
            if strfind(args, 'sl')
                [tind,~,~,BadFitFlag] = run.FitCenterVelocity;
                if BadFitFlag
                    H = 0;
                else
                    H = run.eddy.hcen(tind(1));
                end

                if ff == 1
                    fid = fopen(args, 'w');
                    fprintf(fid, ['|- \n ' ...
                                  '|  | $L$ (km) | $L\\slo$ | $L_z$ (m) | $H_\\text{max}$ (m)| $H$ (m)| ' ...
                                  ' $\\Ro$ | $\\Rh$ | ' ...
                                  '$L\\slo$ (km) | $N^2$ | $f_0$ | $\\beta$ | $\\alpha_\\text{sl}$ ' ...
                                  '| $r_f$ (m/s) | \n |- \n']);
                    fclose(fid);
                end
                fid = fopen(args, 'a');
                fprintf(fid, ['| %s | %2.0f | %3.0f | %3.0f | %3.0f |  %4.0f | ' ...
                              ' %.2f | %2.0f | ' ...
                              ' %1.1e | %1.1e | %1.1e | %1.1e | %1.1e |\n'], ...
                        runName, Lx(1)/1000, bathy.L_slope/1000, Lz(1), max(bathy.h(:)), H, ...
                        Ro(1), run.params.nondim.eddy.Rh, ...
                        N^2, run.params.phys.f0, beta, alpha, run.params.misc.rdrg);
            end

            % sloping shelf runs
            if strfind(args, 'sh')

                Lf = bathy.L_shelf.^2 / (run.params.misc.rdrg/fsb/bathy.sl_shelf);

                [start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
                t0 = start;
                tend = stop;
                [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

                Ldefsh = N * hsb / fsb;
                betash = f0 * run.bathy.sl_shelf./run.bathy.hsb;
                Lbetash = sqrt(V(1)/(betash-beta));

                phi = hsb./(V0./bathy.S_sh/N);
                %chi = (2/sqrt(pi)*exp(-(hsb/Lz0)^2) * V0/Lz0)/(bathy.S_sl*N);
                if ff == 1
                    fid = fopen(args, 'w');
                    fprintf(fid, ['|-\n| %10s | %5s | %4s | %4s | %8s | %5s |' ...
                                  '%6s | %6s | %4s | %5s | %5s | %3s | %5s | \n|-\n'], ...
                            '', '$\Rh$', '$\Ro$', '$\beta$', '$\beta\shl$', ...
                            '$\phi_o$', '$\lambda$','$\Ssl$', '$\Ssh$', '$L_\text{rh}$', ...
                            '$L_\text{def}$', '$r_f$', '$L\atw$');
                    fclose(fid);
                end
                fid = fopen(args, 'a');
                fprintf(fid, ['| %10s | %5.2f | %4.2f | %1.1e |' ...
                              '%.2e | %5.2f | %7.2f |  %5.2f | %0.2f | %5.0f |' ...
                              '%5.0f | %1.1e | %5.0f | \n'], ...
                        run.name, run.params.nondim.eddy.Rh, ...
                        Ro(1), beta, betash, ...
                        phi, hsb./Lz(itsl), bathy.S_sl, bathy.S_sh, ...
                        Lbetash/1000, Ldefsh/1000, ...
                        run.params.misc.rdrg, Lf/1000);
            end

            fclose(fid);
            continue;
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
        if strcmpi(name, 'check max')
            isobath = args(1);

            % m^2/s
            [start,stop] = ...
                run.flux_tindices(run.csflux.off.slope(:,isobath, isobath));

            keyboard;
            vvec = trapz(run.csflux.time(start:stop), ....
                         repnan(run.csflux.off.slopezt(:, start:stop, ...
                                                        isobath, isobath),0), 2);
        end

        if strcmpi(name, 'max flux') | strcmpi(name, 'avg flux')
            default_factor = 1; % integrate to 2xHsb

            if isempty(args)
                isobath = 1;
                factor = default_factor;
            else
                isobath = args(1);
                if length(args) == 2
                    factor = args(2);
                else
                    factor = default_factor;
                end
            end

            if strfind(commands, 'no_sloping_shelf') & (isobath == 1)
                if run.bathy.sl_shelf ~= 0
                    disp(['Skipping ' run.name ': sloping shelf.']);
                    continue;
                end
            end

            if ~isfield(run.csflux, 'off') | (isobath > size(run.csflux.vertbins, 2))
                disp(['Skipping ' run.name]);
                continue;
            end

            source = isobath;
            if ~isinf(factor)
                integrate_zlimit = factor*hsb; % calculate flux above this depth
                fluxvec = run.recalculateFlux(integrate_zlimit, isobath, source);
            else
                integrate_zlimit = run.csflux.h(isobath);
                fluxvec = run.csflux.off.slope(:, isobath, source);
            end

            [maxflux, maxloc, errflx] = run.calc_maxflux(fluxvec,isobath);
            if strfind(name, 'max')
                flux = maxflux;

                t0 = 1;
                tend = maxloc;
            else
                %flux_tref = run.recalculateFlux(integrate_zlimit, 1,1);
                [start,stop] = run.flux_tindices(fluxvec);
                [flux, errflx] = run.calc_avgflux(fluxvec, 0);

                t0 = 1; %start;
                tend = ceil((start+stop)/2);
            end

            [V0, L0, Lz0] = run.EddyScalesForFlux(t0, tend);

            if hsb/Lz0 > 0.5  ...
                    | (strcmpi(run.name, 'ew-2041') & (isobath > 3)) ...
                    | strcmpi(run.name, 'ew-2043') ...
                    | strcmpi(run.name, 'ew-8383')
                warning('skipping because splitting is probably happening.');
                continue;
            end

            eddyscl = V0 * L0 * Lz0;

            zvec = run.csflux.vertbins(:, isobath);
            xvec = run.rgrid.x_rho(1,2:end-1);

            [v,mask] = run.makeStreamerSection(isobath, maxloc, V0, L0, Lz0);
            zind = find_approx(zvec, -abs(integrate_zlimit));
            vmask = v .* mask;

            fluxscl = trapz(zvec(zind:end), ...
                            trapz(xvec, vmask(:,zind:end), 1), 2);

            Lsupp = run.bathy.Lbetash*1.01;
            alpha = run.bathy.sl_shelf;
            if alpha == 0
                % Lsupp = run.bathy.L_shelf;
                continue;
            end
            supplyscl = V0 * Lsupp * (hsb - alpha*Lsupp/2);
            fluxscl = supplyscl;

            debug_fluxes = 0;
            if debug_fluxes
                disp(sprintf([name ' = %.2f mSv, fluxscl = %.2f mSv | ' ...
                             'V0 = %.2f m/s, L0 = %.2f km, Lz0 = %.2f m'], ...
                             flux/1000, fluxscl/1000, V0, L0/1000, Lz0));
            end

            norm(ff) = 1000; eddyscl;
            plotnorm = 1000;

            diags(ff) = flux/plotnorm;
            plotx(ff) = double(fluxscl)/plotnorm;
            err(1,ff) = errflx/plotnorm;

            diagstr = [num2str(flux/1000,'%.2f') '±' ...
                       num2str(errflx/1000,'%.2f') ' mSv'];

            if err(1,ff) ~= 0
                errorbarflag = 1;
            else
                errorbarflag = 0;
            end

            [clr, ptName] = colorize(run, ptName);

            if plotnorm == 1000
                normstr = '(mSv)';
            else
                normstr = '/ Eddy volume flux';
            end

            save_diags = 0;
            mark_outliers = 0; name_points = 1;
            parameterize = 1; logscale = 0;
            force_0intercept = 0;
            line_45 = 0;
            laby = [upper(name(1)) name(2:end) ' at isobath ' normstr];
            labx = ['Predicted Flux Q ' normstr];
            titlestr = [titlestr ' | ND isobath = ' ...
                        num2str(run.csflux.ndloc(:,isobath))];
        end

        if strcmpi(name, 'xpeak')
            iso = args(1);

            if iso > length(run.csflux.x), continue; end
            if ~isfield(run.csflux, 'slopex')
                run.streamerstruct;
            end

            flux = run.csflux.slopex.flux(:, iso, iso);
            xvec = run.csflux.slopex.xi;

            [~,xloc] = max(flux);

            R = run.csflux.R;
            L = run.eddy.vor.dia(1)/2;
            yoR = run.csflux.ndloc(iso); % y/R - used in csflux
            y0oL =  R/L * (1 - yoR); % y0/L - used in derivation
            xfrac = sqrt(1 - y0oL^2);

            name_points = 0;

            diags(ff) = abs(xvec(xloc)) - L;
            plotx(ff) = N/f0 * run.predict_zpeak(iso, 'use');
        end

        if strcmpi(name, 'zpeak')
            debug = 0;

            if debug, figure; title(run.name); hold on; end

            parameterize = 1;
            nsmooth = 3;
            zmax = [];
            kk = 1;
            for iso=args(1)%length(run.csflux.x)
                vec = smooth(run.csflux.off.slopewater.vertitrans(:,iso,iso), ...
                             nsmooth);
                zvec = run.csflux.vertbins(:,iso);

                [~,imax] = max(abs(vec));
                zmax(kk) = zvec(imax);
                if debug, plot(vec, zvec); end
                kk = kk + 1;
            end

            if length(zmax) > 1
                E = [run.csflux.ndloc(2:end)' ones(size(run.csflux.ndloc(2:end)))'];
                [P,Pint,R,Rint,stats] = regress(zmax', E, 0.05);

                errorbarflag = 0;
                err(1,ff) = Pint(2) - P(2);
            else
                P(2) = zmax;
            end

            [~,maxloc] = run.calc_maxflux(args(1));
            fvec = run.rgrid.f(:,1);
            fratio = fvec(run.csflux.ix(iso))./fvec(run.bathy.isb);

            [~,zpeak] = run.predict_zpeak(iso, []);
            diags(ff) = abs(P(2));
            plotx(ff) = abs(zpeak);
            % plotx(ff) = (Ro(1)) * Lz(maxloc)/hsb * run.csflux.ndloc(iso);

            name_points = 1;
            labx = 'Ro H_{sb}';
            laby = 'z_{peak} - H_{sb}/2';

            if debug
                liney(-diags(ff));title(run.name);
            end
        end

        if strcmpi(name, 'streamerwidth')
            % use run.streamer.xprof and estimate peak width at half maximum

            xprof = run.streamer.xprof(:,1);
            xivec = run.streamer.xivec;
            [xmax,imax] = max(xprof);

            % find half-maximum
            hm(1) = find_approx(xprof(1:imax), xmax/2, 1);
            hm(2) = imax-1 + find(xprof(imax:end) <= xmax/2, 1, 'first');

            % peak width at half maximum
            pwhm = abs(diff(xivec(hm)));

            [start,stop] = run.flux_tindices(1);
            [~,L0,~] = run.EddyScalesForFlux(start,stop);

            diags(ff) = run.supply.zeta.xscale*2/1000;
            plotx(ff) = run.sbssh.X/1000*2;
            %run.supply.xscale/1000 - run.supply.IntersectScale/1000; bathy.Lbetash/1000; abs(L0)/1000;

            parameterize = 1;
        end

        if strcmpi(name, 'zpeakwidth')
            parameterize = 1;
            iso = args(1);

            if iso > size(run.csflux.off.slope, 3), continue; end

            [~,~,zwidth] = run.streamer_peak(iso);
            prediction = run.predict_zpeak(iso, []);

            name_points = 1;

            diags(ff) = abs(zwidth(1));
            plotx(ff) = abs(prediction);

            labx = 'Prediction (m)';
            laby = 'Observed peak width (m)';
            force_0intercept = 0;
            % continuity
            %Lz0 = Lz(maxloc);
            %fn = hsb * (1-erf(hsb/Lz0)) + Lz0/sqrt(pi)*(1 - exp(-(hsb/Lz0)^2));
            %fn = hsb + Lz0/sqrt(pi);
            %plotx(ff) = Ro(1) * fn/hsb * R/Ly;
        end

        if strcmpi(name, 'streamervel')
            [avgflux, ~] = run.calc_avgflux;
            [maxflux, maxloc] = run.calc_maxflux;

            vxt = run.csflux.shelfxt(:,:,1) .* ...
                  run.csflux.offmask./run.bathy.hsb;

            diags(ff) = mean(vxt(:,maxloc) .* (vxt(:,maxloc) > 0));
            plotx(ff) = V(tind);

            parameterize = 1;
            laby = 'Mean streamer velocity (m/s)';
            labx = 'V(tind)';
        end

        % width of outflow supply jet
        if strcmpi(name, 'outflow')
            if run.bathy.sl_shelf == 0, continue; end

            [~,maxloc] = run.calc_maxflux(1);
            V = run.eddy.rhovor.Vke(maxloc);
            beta_t = run.bathy.sl_shelf * f0/run.bathy.hsb;

            [start,stop]  = run.flux_tindices(run.csflux.off.slope(:,1,1));
            start = ceil((start+stop)/2);
            width = run.bathy.xsb - ...
                    nanmedian(smooth(run.csflux.off.slopewater.envelope(start:stop,1), 10));

            diags(ff) = width/1000;
            plotx(ff) = sqrt(V/beta_t)/1000;

            laby = 'Source (km)';
            labx = '$$\sqrt{U/\beta_t}$$ (km)';
            name_points = 1;
        end

        % width of inflow supply jet
        if strcmpi(name, 'inflow')
            if run.bathy.sl_shelf == 0, continue; end

            [~,maxloc] = run.calc_maxflux(1);
            V = run.eddy.rhovor.Vke(maxloc);
            beta_t = run.bathy.sl_shelf * f0/run.bathy.hsb;

            [start,stop]  = run.flux_tindices(run.csflux.off.slope(:,1,1));
            start = ceil((start+stop)/2);
            width = nanmean(smooth(run.onshelf.edd.env(start:stop), 10));

            diags(ff) = width/1000;
            plotx(ff) = sqrt(V/beta_t)/1000;

            laby = 'Source (km)';
            labx = '$$\sqrt{U/\beta_t}$$ (km)';
            name_points = 1;
        end


        %%%%%%%%%%%%%%%%%%%% deprecated

        %%%%%
        if strcmpi(name, 'aspect')
            % critical aspect ratio hypothesis
            run.fit_traj(1.1);

            tind = run.traj.tind;

            H = hcen(tind);

            diags(ff) = H./Lz(tind);
            plotx(ff) = beta/beta_t;

            laby = 'H/Lz';
            labx = 'index';
        end

        %%%%% topography parameter - taylor column
        if strcmpi(name, 'hogg')
            diags(ff) = 1./Ro(1) ./ ...
                (1+hsb/Lz(1) * 1/Sa);
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

        %%%%% test critical iflux hypothesis for eddy to
        %%%%% start moving northward
        if strcmpi(name, 'critical flux')
            iflux = run.csflux.off.itrans.shelf(:,1);
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
            %       run.csflux.off.shelf(:,1));
            %linex(run.eddy.t(ind));

            %run.animate_zeta(ind, 1);

            % Get Flux at ind
            run.csflux.critrans = ...
                run.csflux.off.itrans.shelf(ind,1);

            diags(ff) = run.csflux.critrans;
        end

        %%%%% resistance
        %%%%% old version of bottom torque
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

        if isempty(diagstr)
            diagstr = num2str(diags(ff));
            if ~isnan(err(1,ff))
                diagstr = [diagstr '±' num2str(err(1,ff))];
            end
        end

        if plots
            figure(hfig);
            axes(hax);

            %[clr,ptName] = colorize(run, ptName);

            if marker == '.'
                markersize = 28;
                if strfind(commands, 'small_points')
                    markersize = 14;
                end
            else
                markersize = 18;
            end

            if errorbarflag
                if ~isnan(err(1,ff)) && isnan(err(2,ff))
                    err(2,ff) = err(1,ff);
                end
                errorbar(plotx(ff), diags(ff), abs(err(1,ff)), abs(err(2,ff)), ...
                         'x', 'LineWidth', 2, 'Color', clr, ...
                         'Marker', marker, 'MarkerSize', markersize);
            else
                plot(plotx(ff), diags(ff), 'Marker', marker, 'Color', clr, ...
                     'MarkerSize', markersize);
            end

            % add run names
            if (name_points | strfind(commands, 'name_points')) ...
                    & isempty(strfind(commands, 'no_name_points'))
                if strip_ew
                    ll = strfind(ptName, 'ew-');
                    if ~isempty(ll)
                        ptName = ptName(ll+3:end);
                    end
                end

                text(plotx(ff), diags(ff), ptName, 'FontSize', ...
                     12, 'Rotation', 0, 'Color', clr, ...
                     'VerticalAlignment','Bottom');
            end
            if ff == 1 | isnan(diags(1))
                if ~iscell(labx) & strfind(labx, '$')
                    xlabel(labx, 'interpreter', 'latex');
                else
                    xlabel(labx);
                end
                if ~iscell(laby) & strfind(laby, '$')
                    ylabel(laby, 'interpreter', 'latex')
                else
                    ylabel(laby);
                end
                htitle = title(titlestr);
            end
        end

        color(ff,:) = clr;
        disp([num2str(ii,'%02d') ' | ' run.name ' | ' name ' = ' diagstr ' | plotx ' ...
              '= ' num2str(plotx(ff))]);
    end

    if exist('sortedflag', 'var') & sortedflag
        runArray.reset_colors(cb);
    end

    if save_diags
        savename = name;
        savename(isspace(name)) = [];
        save(['./diags/diags_' savename '.mat']);
    end

    if parameterize && plots
        limx = xlim;
        limy = ylim;

        xvec = linspace(0.5*nanmin(plotx(:)), 1*nanmax(plotx), 100);

        err = max(abs(err), [], 1);
        % if strcmpi(name, 'bottom torque')
        %     err = fillnan(err, 0); % 6362-1
        %     diags(isnan(err)) = NaN;
        %     plotx(isnan(err)) = NaN;
        % end

        if exist('exceptions', 'var')
            diags(exceptions) = NaN;
            plotx(exceptions) = NaN;
            err(exceptions) = NaN;
        end

        if ~force_0intercept
            E = [cut_nan(plotx') ones(size(cut_nan(plotx')))];
        else
            E = cut_nan(plotx');
        end

        if ~isempty(cut_nan(err)) & ~strcmpi(name, 'max flux')
            disp('Using weighted least squares');
            disp(['% error: ' num2str(round(100*err./diags))]);

            [P,stderror,MSE] = lscov(E, cut_nan(diags'), 1./cut_nan(err));
            % The weights are *not* the exact covariance matric of B
            %stderror = stderror * sqrt(1/MSE);
            tval = abs(conft(0.05,length(diags)));
            Perr = tval * stderror;

            Pint = [P-Perr, P+Perr];
        else
            [P,Pint,R,Rint,stats] = regress(cut_nan(diags'), E, 0.05);

            Perr(1) = P(1) - Pint(1,1);
            Perr(2) = P(2) - Pint(2,1);

            % diagnose outliers
            outind = logical(zeros(size(plotx)));
            for ii=1:size(Rint,1)
                if ~((Rint(ii,1) <= 0) && (Rint(ii,2) >= 0))
                    outind(ii) = true;
                end
            end

            % mark outliers
            if mark_outliers
                plot(plotx(outind), diags(outind), 'ro', 'MarkerSize', 13);
            end
        end

        % R - value @ 99%
        [r,Pcorr,rlo,rhi] = corrcoef(cut_nan(diags), cut_nan(plotx), 'alpha', 0.01);
        stats(1) = r(2)^2; % return R² always.
        rerr = rhi(2)-r(2);

        if rlo * sign(r) < 0
            warning('Correlation not significant at 99%');
        end

        if force_0intercept
            P(2) = 0;
            Pint(2,1:2) = 0;
        end

        % [Q,R] = qr(E,0);
        % S = inv(R);
        % ETEI = S*S';
        % % assuming noise vector (res) is white
        % P2 = ETEI * E' * var(res) * E * ETEI;
        % errfit = sqrt(diag(P2));
        % errfit = errfit(2); % standard error

        c = P(1);
        rmse = sqrt(nanmean((diags - c*plotx - P(2)).^2));

        hparam = plot(xvec, c*xvec + P(2), '--', ...
                      'Color', [1 1 1]*0.5);
        slopestr = num2str(c, '%.2f');
        intstr = num2str(P(2), '%.2f');

        if exist('Pint', 'var')
            hparam(2) = plot(xvec, Pint(1,1) * xvec + Pint(2,1), ...
                             'Color', [1 1 1]*0.5, 'LineStyle', '--');
            hparam(3) = plot(xvec, Pint(1,2) * xvec + Pint(2,2), ...
                             'Color', [1 1 1]*0.5, 'LineStyle', '--');

            slopestr = [slopestr '\pm' ...
                        num2str(abs(Pint(1,1)-P(1)), '%.2f')];
            intstr = [intstr '\pm' ...
                      num2str(abs(Pint(2,1)-P(2)), '%.2f')];
        end

        if force_0intercept, intstr = '0'; end
        %hleg = legend(hparam(1), ['y = (' slopestr ') x + (' ...
        %                    intstr '); rmse = ' num2str(rmse,3)], ...
        %              'Location', 'SouthEast');
        uistack(hparam, 'bottom');

        tlen = 0;
        if exist('stats', 'var')
            disp(['stats = ' num2str(stats)]);
            textstr{1} = ['C_r = ' num2str(sqrt(stats(1)), '%.2f')];
            %textstr{2} = ['rmse = ' num2str(rmse, '%.2f')];
            tlen = length(textstr);
        end
        textstr{tlen+1} = ['m = ' slopestr];
        textstr{tlen+2} = ['c = ' intstr];

        htext = text(0.03,0.92,textstr, ...
                     'FontSize', 16, ...
                     'Units', 'normalized', ...
                     'HorizontalAlignment', 'left', ...
                     'VerticalAlignment', 'top');
    end

    if plots
        figure(hfig);
        maximize(); pause(0.05);
        set(gcf, 'renderer', 'painters');
        if line_45, line45(hax); end

        if logscale
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
        end

        pbaspect([1.618 1 1]);

        if strcmpi(name, 'bottom torque') & strcmpi(args, 'functional')
            c = 0.01;
            x = linspace(4e-3, 0.15, 100);
            y = erfinv(1 - 1.45* (x./(1-x)) - 0.02);
            y0 = erfinv(1 - (1.45-0.15)*(x./(1-x)) - 0.01);
            y1 = erfinv(1 - (1.45+0.15)*(x./(1-x)) - 0.03);
            hold on;
            hparam(1) = plot(x, y, '--','Color', [1 1 1]*0.3);
            hparam(2) = plot(x, y0, '--', 'Color', [1 1 1]*0.3);
            hparam(3) = plot(x, y1, '--', 'Color', [1 1 1]*0.3);

            uistack(hparam, 'bottom');
            hleg = legend(hparam, ['$$1 - \mathrm{erf}(\frac{H}{L_z^0}) ' ...
                                '= (1.45 \pm 0.15)  \frac{\beta}{\beta_z - \beta}' ...
                                '+ (0.02\pm0.01)$$'], ...
                          'Location', 'NorthEast');
            hleg.FontSize = 28;
            set(hleg, 'interpreter', 'latex');
        end

        if strcmpi(name, 'bottom torque') ...
                & (strcmpi(args, 'regression') | isempty(args))
            hleg = legend(hparam, ['$$ \frac{U_b}{U_s} = ' ...
                            '1 - \mathrm{erf}(\frac{H}{L_z^0}) ' ...
                                '= ' num2str(c,3) ' \frac{\beta}{\beta_z - \beta}' ...
                                '+ ' num2str(P(2),2) '$$'], ...
                      'Location', 'SouthEast');
            set(hleg, 'interpreter', 'latex');
        end

        if strcmpi(name, 'avg flux') || strcmpi(name, 'max flux') ...
                || strcmpi(name, 'zpeakwidth');
            figure(hfig);
            % axis square;
            xlim([0 max(plotx)]);
            ylim([0 max(ylim)]);

            hiso = text(0.68,0.15, ...
                        ['y/R = ' ...
                         num2str(run.csflux.ndloc(isobath), '%.2f')], ...
                        'Units', 'normalized', ...
                        'HorizontalAlignment', 'left', ...
                        'VerticalAlignment', 'top', 'FontSize', 16);
        end

        beautify([24 26 28]);

        if kozak
            % call after beautify
            MakeKozak;
        end
    end

    if exist('htext', 'var')
        handles.htext = htext;
    end
    if exist('hiso','var')
        handles.hiso = hiso;
    end
    if exist('hparam','var')
        handles.hparam = hparam;
    end
    if exist('hleg','var')
        handles.hleg = hleg;
    end

    if strcmpi(name, 'params table')
        diags = '|- \n';

        fid = fopen(args, 'a');
        fprintf(fid, diags);
        fclose(fid);

        type(args);
    end
end

% colorize points
function [clr, ptName, marker] = colorize(run, ptName)
    clr = [0 0 0]; %[96 125 139]/255;
    marker = '.';

    try
        if run.params.flags.conststrat == 0
            ptName = ' N^2';
            clr = [231,41,138]/255;
            return;
        end
    catch ME
    end

    if run.bathy.sl_shelf ~= 0
        clr = [255 87 34]/255;
    end

    if run.params.misc.rdrg ~= 0
        clr = [63 81 181]/255;
        marker = '+';
    end

    if (run.bathy.sl_shelf ~=0) ...
            & (run.params.misc.rdrg ~= 0)
        clr = [255 160 0]/255;
    end

    if run.params.eddy.tamp < 0
        ptName = ' C';
        clr = [117,112,179]/255;
        marker = 'o';
    end

    % if run.params.phys.f0 < 0
    %     ptName = ' f^{-}';
    %     clr = [27,158,119]/255;
    % end

    if ~isempty(strfind(run.name, 'ew-64461-9')) | ...
            ~isempty(strfind(run.name, 'ew-64461-8'))
        marker = 'o';
    end
end