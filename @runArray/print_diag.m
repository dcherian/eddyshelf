function [diags, plotx] = print_diag(runArray, name, args, hax)

    if ~exist('args', 'var'), args = []; end
    if ~exist('hax', 'var'), hax = []; end

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    if ~strcmpi(name, 'nondim') && ~strcmpi(name, 'beta gyre') ...
            && ~strcmpi(name, 'bfric') && ~strcmpi(name, 'hbl') ...
            && ~strcmpi(name, 'arrest')
        plots = 1;
    else
        plots = 0;
    end

    annostr = ['runArray.print_diag(' name ')'];
    diags = nan(size(runArray.filter));
    plotx = diags;

    if plots
        if isempty(hax)
            hfig = figure;
            hax = gca; % default axes
        else
            hfig = gcf;
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
        clr = 'k';
        error = [];
        parameterize = 0;
        titlestr = name;
    end

    for ff=1:length(runArray.filter)
        ii = runArray.filter(ff);
        run = runArray.array(ii);
        runName = runArray.getname(ii);

        ptName = runName;
        clr = 'k'; % reset color

        % some commonly used variables
        tind = run.tscaleind;
        [~,uind,~] = unique(run.time, 'stable');
        ndtime = run.eddy.t * 86400 / run.eddy.turnover;
        Lx = run.eddy.vor.lmaj;
        Ly = run.eddy.vor.lmin;
        Lz = run.eddy.Lgauss;
        Ls = run.eddy.Ls;
        Ro = run.eddy.Ro;
        V = run.eddy.V;
        if isfield(run.eddy, 'Vb'), Vb = run.eddy.Vb; end

        hedge = run.eddy.hedge;
        hcen = run.eddy.hcen;
        fcen = run.eddy.fcen;
        my = run.eddy.my;
        mx = run.eddy.mx;

        A = run.eddy.amp(1);

        Lr = run.rrdeep;
        beta = run.params.phys.beta;
        use run.params.phys
        f0 = abs(f0);

        alpha = run.bathy.sl_slope;
        hsb = run.bathy.hsb;
        hsl = run.bathy.hsl;
        xsb = run.bathy.xsb;
        xsl = run.bathy.xsl;
        Sa = run.bathy.S_sl;
        beta_t = f0 * alpha./Lz(1);
        N = sqrt(run.params.phys.N2);
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
            eval(['error(ff) = abs(error_' en ');']);

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

            tanhfitflag = 0;
            if tanhfitflag %|| run.params.eddy.tamp < 0
                run.fit_traj()

                tind = run.traj.tind;
                Y = run.eddy.my(tind); %run.traj.Y;

                titlestr = ['tanh(t/T) at t = ' ...
                            num2str(run.traj.tcrit) 'T'];
            else
                nsmooth = 10;
                factor = 1/3;
                [~,~,tind] = run.locate_resistance(nsmooth, factor);
                Y = run.eddy.my(tind) - run.bathy.xsb; ...
                    run.eddy.my(tind);
                titlestr = ['Water depth at which cross-isobath translation ' ...
                            'velocity drops by ' num2str(1-factor,2)];
            end
            if isempty(tind), continue; end

            H = run.eddy.hcen(tind);
            t0 = 1; %run.eddy.tscaleind;
            titlestr = [titlestr ' | t0 = ' num2str(t0)];

            Lz0 = Lz(1); %*run.eddy.grfactor(1);
            beta_t = (alpha*abs(fcen(tind))/Lz(tind));

            zz = H./Lz0;
            diags(ff) = (1-erf(zz)); %./((zz*(1-erf(zz)) + 1/sqrt(pi) * ...
                                     %(1 - exp(-zz^2))));
            plotx(ff) = beta/beta_t./(1-beta/beta_t);
            errorbarflag = 0; kozak = 1;
            name_points = 1; line_45 = 0;
            parameterize = 1;
            slope = num2str(round(alpha*1e3));

            if errorbarflag
                error(ff) = 2/sqrt(pi) * exp(-(H/Lz(t0))^2) * ...
                    (alpha * run.traj.yerr)/Lz(t0);
            end

            try
                if run.params.flags.conststrat == 0
                    ptName = ' N^2';
                    clr = [231,41,138]/255;
                end
            catch ME
                if run.params.eddy.tamp < 0
                    ptName = ' C';
                    clr = [117,112,179]/255;
                else
                    if run.params.phys.f0 < 0
                        ptName = ' f^{-}';
                        clr = [27,158,119]/255;
                    else
                        clr = 'k';
                        ptName = ''; %[runName];
                    end
                end
                %ptName = num2str(err,2);
            end
            % mark NS isobath runs with gray points
            if run.bathy.axis == 'x'
                clr = [1 1 1]*0.75;
            end

            laby = '$$\frac{U_b}{U_s} = 1 - \mathrm{erf}(\frac{H}{L_z^0})$$';
            labx = '$$\beta/\beta_t$$';
            %if ff == 1, hax = subplot(121); hold all; end
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

        % compare against flat bottom rest latitude
        if strcmpi(name, 'rest latitude')
            [~,~,tind] = run.locate_resistance();

            sgn = sign(f0) * sign(run.params.eddy.tamp);
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

            if ff == 1
                close;
                disp(sprintf('| %18s | %5s | %4s | %4s | %8s | %5s | %6s |', ...
                             'name', 'Rh', 'Ro', 'bl/f', 'beta_t', ...
                             'L/Lsl', 'lambda'));
            end
            disp(sprintf('| %18s | %5.2f | %4.2f | %4.2f | %.2e | %5.2f | %5.2f |', ...
                         run.name, run.params.nondim.eddy.Rh, ...
                         Ro(1), beta*Lx(1)/f0, beta_t, ...
                         Lx(1)./run.bathy.L_slope, ...
                         hsb./Lz(1)));
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
        if strcmpi(name, 'max flux')
            if isempty(args)
                isobath = 3;
            else
                isobath = args(1);
            end

            if ~isfield(run.csflux, 'time') || ...
                    isobath > size(run.csflux.west.slope, 2)
                disp(['Skipping ' run.name]);
                continue;
            end

            [maxflux, maxloc] = ...
                run.calc_maxflux(run.csflux.west.slope(:,isobath, isobath));

            if isnan(maxflux), continue; end

            %if ~isfield(run.eddy, 'rhovor'), continue; end

            % maxflux = max(run.csflux.west.slope(:,isobath, isobath));
            H = run.csflux.h(isobath);

            tind = 1; maxloc;
            %[V1, L, x0] = run.fit_vel(tind);
            % V1/V0 = 2.3 quite dependably
            %V0 = V1*2.3;
            V0 = run.eddy.V(tind) * 2.33;
            % distance of eddy center from shelfbreak
            R = run.csflux.R;
            % vor.dia is radius to max vel =
            L = run.eddy.Ls(1)/2;
            Lz0 = Lz(maxloc);

            % I can make slope ~ 1 by
            % V0 * 0.1
            % xfrac * 2
            yoR = run.csflux.ndloc(isobath); % y/R - used in csflux
            y0oL =  R/L * (1 - yoR); % y0/L - used in derivation
            xfrac = sqrt(1 - y0oL^2);

            syms x z;
            %            fluxscl = -1 * abs(V0) * ...
            %          int(x/L*exp(-(x/L)^2), -Inf, -xfrac*L) * exp(-y0oL^2) ...
            %          * int(1 - erf(-z/Lz0), z, -H, 0);
            fluxscl = abs(V0) * L/2 * exp(-xfrac^2) * exp(-y0oL^2) ...
                      * int(1 - erf(-z/Lz0), z, -H, 0);
            %fluxscl = abs(V0) * L * exp(-xfrac)* exp(-y0oL^2) ...
            %          * int(1 - erf(-z/Lz0), z, -H, 0);

            norm = 1e3;
            transscl = double(fluxscl)/norm;

            % colorize
            crit = 0.40;
            if round(hsb/Lz(1),2) > crit
                clr = 'r';
            end
            if round(hsb/Lz(1),2) < crit
                clr = 'b';
            end

            parameterize = 1; logscale = 0;
            force_0intercept = 0;
            errorbarflag = 0; name_points = 1; line_45 = 0;
            laby = 'Slope water max flux (mSv)';
            labx = 'Volume flux in eddy (mSv)';
            titlestr = [titlestr ' | ND isobath = ' ...
                        num2str(run.csflux.ndloc(:,isobath))];

            diags(ff) = maxflux/norm;
            plotx(ff) = transscl;
            error(ff) = 0.10*diags(ff);
        end

        if strcmpi(name, 'avg flux')
            isobath = 3;

            if ~isfield(run.csflux, 'time') || ...
                    isobath > size(run.csflux.west.slope, 2)
                continue;
            end

            [avgflux, err] = ...
                run.calc_avgflux(run.csflux.west.slope(:,isobath, isobath));

            H = run.csflux.h(isobath);
            R = run.csflux.R;
            V0 = V(1)/ 0.43;
            L = run.eddy.vor.dia(1)/2;
            Lz0 = Lz(1);

            yoR = run.csflux.ndloc(isobath); % y/R - used in csflux
            y0oL = R/L * (1 - yoR); % y0/L - used in derivation
            xfrac = sqrt(1 - y0oL^2);

            syms x z;
            fluxscl = -1 * V0 * ...
                      int(x/L*exp(-(x/L)^2), -Inf, -xfrac*L) * exp(-y0oL^2) ...
                      * int(1 - erf(z/Lz0), z, -H, 0);

            %transscl = 0.075 * g/f0 .* run.eddy.amp(ind) .* hsb/1000;
            %transscl = 0.023 * Lx(1) * V(tind) * hsb / 1000 + 0.0;
            transscl = double(fluxscl); %V(1) * H * Lx(1);

            diagstr = [num2str(avgflux/1000,'%.2f') '±' ...
                       num2str(err/1000,'%.2f') ' mSv | scale = ' ...
                       num2str(transscl)];

            % colorize
            if run.bathy.S_sh ~= 0
                % ptName = 'sh';
                clr = 'g';
            end
            %if run.params.misc.rdrg ~= 0
            %    % ptName = 'f';
            %    clr = 'r';
            %end

            if round(hsb/Lz(1),2) > 0.13
                clr = 'r';
            end
            if round(hsb/Lz(1),2) < 0.13
                clr = 'b';
            end

            %if run.params.nondim.eddy.Rh < 10
            %    clr = 'b';
            %end

            % plot error
            %if ff == 1 || ~exist('hfig_fluxerr', 'var')
            %    hfig_fluxerr = figure; hold on;
            %end
            % paramerr = avgflux/1000 - transscl;
            % figure(hfig_fluxerr)
            % plot(run.eddy.Ro(tind), paramerr, '*');
            % xlabel('Ro'); ylabel('paramerr');

            parameterize = 1;
            errorbarflag = 1; name_points = 1; line_45 = 0;
            laby = 'Flux (mSv)';
            labx = 'Parameterization (mSv)';
            titlestr = [titlestr ' | ND isobath = ' ...
                        num2str(run.csflux.ndloc(:,isobath))];

            diags(ff) = avgflux/1000;
            error(ff) = err/1000;
            plotx(ff) = transscl;
        end

        if strcmpi(name, 'streamervel')
            [avgflux, ~] = run.calc_avgflux;
            [maxflux, maxloc] = run.calc_maxflux;

            vxt = run.csflux.shelfxt(:,:,1) .* ...
                  run.csflux.westmask./run.bathy.hsb;

            diags(ff) = mean(vxt(:,maxloc) .* (vxt(:,maxloc) > 0));
            plotx(ff) = V(tind);

            parameterize = 1;
            laby = 'Mean streamer velocity (m/s)';
            labx = 'V(tind)';
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
        end

        if plots
            figure(hfig);
            axes(hax);
            if errorbarflag
                errorbar(plotx(ff), diags(ff), error(ff), 'x', ...
                         'LineWidth', 2, 'Color', clr);
            else
                plot(plotx(ff), diags(ff), '.', 'Color', clr, ...
                     'MarkerSize', 20);
            end

            % add run names
            if name_points
                if strip_ew, ptName = ptName(4:end); end

                text(plotx(ff), diags(ff), ptName, 'FontSize', ...
                     12, 'Rotation', 0, 'Color', clr, ...
                     'VerticalAlignment','Bottom');
            end
            if ff == 1
                if strfind(labx, '$$')
                    xlabel(labx, 'interpreter', 'latex');
                else
                    xlabel(labx);
                end
                if strfind(laby, '$$')
                    ylabel(laby, 'interpreter', 'latex')
                else
                    ylabel(laby);
                end
                htitle = title(titlestr);
            end
        end

        disp([num2str(ii,'%02d') ' | ' run.name ' | ' name ' = ' diagstr ' | plotx ' ...
              '= ' num2str(plotx(ff))]);
    end

    if exist('sortedflag', 'var') & sortedflag
        runArray.reset_colors(cb);
    end

    if parameterize && plots
        limx = xlim;
        limy = ylim;

        xvec = linspace(0.5*min(plotx(:)), 1*max(plotx), 100);

        if ~force_0intercept
            E = [plotx' ones(size(plotx'))];

            %P = polyfit(cut_nan(plotx), cut_nan(diags), 1);
            [P,Pint,R,Rint,stats] = regress(diags', E, 0.05);
        else
            E = cut_nan(plotx');
            P(1) = E\cut_nan(diags');
            P(2) = 0;
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
                      'Color', [1 1 1]*0.75);
        slopestr = num2str(c, '%.2f');
        intstr = num2str(P(2), '%.2f');

        if exist('Pint', 'var')
            hparam(2) = plot(xvec, Pint(1,1) * xvec + Pint(2,1), ...
                             'Color', [1 1 1]*0.75);
            hparam(3) = plot(xvec, Pint(1,2) * xvec + Pint(2,2), ...
                             'Color', [1 1 1]*0.75);

            slopestr = [slopestr '\pm' ...
                        num2str(abs(Pint(1,1)-P(1)), '%.2f')];
            intstr = [intstr '\pm' ...
                      num2str(abs(Pint(2,1)-P(2)), '%.2f')];
        end

        hleg = legend(hparam(1), ['y = (' slopestr ') x + (' ...
                            intstr '); rmse = ' num2str(rmse,3)], ...
                      'Location', 'SouthEast');

        uistack(hparam, 'bottom');

        if exist('stats', 'var')
            disp(['stats = ' num2str(stats)]);
    end

    if plots
        figure(hfig);
        maximize(); drawnow; pause(1);
        set(gcf, 'renderer', 'painters');
        if line_45, line45(hax); end

        if logscale
            set(gca, 'XScale', 'log');
            set(gca, 'YScale', 'log');
        end

        if kozak
            ax2 = kozakscatterplot(hax, [min(plotx) max(plotx)], ...
                             [min(diags) max(diags)]);
            grid off;
            %hleg = legend;
            %pos = hleg.Position;
            %hleg.Position = [pos(1)-0.1 pos(2:end)];
        else
            pbaspect([1.618 1 1]);
        end

        if strcmpi(name, 'bottom torque')
            hleg = legend(hparam, ['$$ \frac{U_b}{U_s} = ' ...
                            '1 - \mathrm{erf}(\frac{H}{L_z^0}) ' ...
                            '= ' num2str(c,3) ' \frac{\beta}{\' ...
                            'beta_t} + ' num2str(P(2),2) '$$ ' ...
                            '; rmse = ' num2str(rmse,3)], ...
                      'Location', 'SouthEast');
            set(hleg, 'interpreter', 'latex');
        end

        if strcmpi(name, 'avg flux') || strcmpi(name, 'max flux')
            figure(hfig);
            axis square;
            xlim([0 max(xlim)]);
            ylim([0 max(ylim)]);
        end

        beautify([20 24 30]);
    end
end