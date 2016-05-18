%  animate 2D field at the surface presumably.
function [handles] = animate_field(runs, name, hax, t0, ntimes, opt)
    if ~exist('hax','var'), hax = []; end

    runs.video_init(name);
    titlestr = [];

    fontsize = [];

    % zoom in axis?
    limx = [165 410]; limy = [0 120];
    limx = []; limy = [];

    commands = '';

    dt = 1; dx = 0; dy = 0;
    factor = 1; % scale variable (1 by default)

    csfluxplot = 0; % 0 = no flux plot
                    % 1 = instantaneous x-profile;
                    % 2 = time series plot
    csfluxIsobath = 4; % for profile
    csfluxFinalize = 0; % for csfluxplot = 2

    asfluxplot = 0; % 0 = no flux plot
                    % 1 = instantaneous y-profile;
                    % 2 = time series plot : left, right and total
    if (asfluxplot == 1) || (asfluxplot == 2) % which location for asflux plot?
        asindex = [1 2];
    end

    % time series?
    enfluxplot = 0; % plot AS energy flux ?
    vecplot = 0; % plot some time vector (assign tvec and vec);
    pointplot = 0; % mark some point on the map

    % eddy contours?
    rhocontourplot = 1; % plot eddy drho contour too?
    vorcontourplot = 0; % vorticity contour
    sshplot = 0; % plot ssh-contour too?
    dyeplot = 0; % plot eddye contour too?

    % extra contours
    addcsdye = 0; % add csdye to eddye plot?
    addzeta = 0; % overlay zeta contours
    stopzeta = Inf; % stop display zeta

    csdcontourplot = 1; % contour csd contours
    try
        isobath = [1 7];
        isobath(isobath > length(runs.csflux.x)) = [];
        csdcontours = runs.csflux.x(isobath);
    catch ME
        if runs.params.flags.flat_bottom
            csdcontours = [0 15 30 45 60] * 1000;
        else
            csdcontourplot = 0;
        end
    end
    modify_bathy = 1; % modify bathymetric contours to show
                      % isobaths corresponding to csdcontours?
    bathycolor = [1 1 1]*0.45;
    addvelquiver = 0; % velocity vectors?
    dxi = 5; dyi = 5; % decimate vectors
    uref = 1; vref = uref; % scale vectors
    scale = 3;
    normquiver = 0; % normalize so that magnitude is one everywhere
    quivercolor = 'k';

    % graphics
    nocolorbar = 0;
    AnimateZoom = 0;
    ZoomStart = 1;
    ZoomEnd = 1;
    ZoomXLimStart = [];
    ZoomXLimEnd = [];
    ZoomYLimStart = [];
    ZoomYLimEnd = [];

    % eddy diagnostics
    drawtrack = 1; % plot eddy track?
    drawedgetrack = 0; % plot eddy's edge tracks?
    drawcenter = 1; % mark eddy center location.

    % grid diagnostics
    telesplot = 0;  % plot lines where grid stretching starts
                    % and ends

    if exist('opt', 'var') & ~isempty(opt)
        fields = fieldnames(opt);
        for ff = 1:length(fields)
            eval([fields{ff} ' = opt.' fields{ff} ';']);
        end
    end

    if ~isfield(runs.eddy, 'rhovor')
        rhocontourplot = 0;
    end

    if csfluxplot == 2
        vecplot = 1;
        tvec = runs.csflux.time/86400;
        for kkk=1:length(csfluxIsobath)
            if csfluxIsobath(kkk) == 1
                isobathNames{kkk} = [runs.bathy.axis ' = shelfbreak'];
            else
                isobathNames{kkk} = [runs.bathy.axis ' = ' ...
                                    num2str(round(runs.csflux.x(csfluxIsobath(kkk))/1000)) ...
                    ' km'];
            end
            vec(:,kkk) = runs.csflux.off.slope(:,csfluxIsobath(kkk), ...
                                               csfluxIsobath(kkk));
        end
        if kkk == 1 & csfluxIsobath == 1
            laby = {'Shelf water flux'; '(m^3/s)'};
        else
            laby = {'Shelf-slope water'; 'flux (m^3/s)'};
        end
        locy = runs.bathy.xsb/1000; locx = [];
    end
    %%%%%%%%%%%%%%%%%%% options processing ends

    if pointplot
        px = (runs.eddy.mx + runs.eddy.vor.dia(1))/1000;
        py = runs.eddy.my/1000;
    end

    if vecplot & (csfluxplot ~= 2)
        try
            %%% integrated energy asflux
            %tvec = runs.eddy.t;
            %vec = (runs.asflux.ikeflux(:,3) + runs.asflux.ipeflux(:,3) ...
            %      - runs.asflux.ikeflux(:,2) - runs.asflux.ipeflux(:,2));
            %laby = 'Integrated energy flux';
            %locx = runs.asflux.x(2:3)/1000; locy = [];

            %%% dE/dt
            %tvec = avg1(runs.eddy.t);
            %vec = smooth(diff(runs.eddy.KE + runs.eddy.PE)./ ...
            %             diff(runs.eddy.t*86400), 4);
            %laby = 'dE/dt';
            %locx = []; locy = [];

            %tvec = runs.eddy.t;
            %vec = runs.eddy.cvy;
            %laby = 'cvy'
            %locx = []; locy = [];

            %%% asflux time series
            %tvec = runs.asflux.time/runs.tscale;
            %vec = runs.asflux.ikeflux(:,2);
            %locx = runs.asflux.x(2); locy = [];

            %%% area plot
            %vec = runs.eddy.vor.lmin .* runs.eddy.vor.lmaj;
            %tvec = runs.eddy.t;
            %laby = 'Surface area (m^2)';
            %locx = []; locy = [];

            %%% length scales
            vec = runs.eddy.rhovor.dia;
            tvec = runs.eddy.t;
            laby = 'Diameter';
            locx = []; locy = [];
        catch ME
            disp('vecplot failed!');
            vecplot = 0;
        end
    end

    % over flat bottom, addcsdye makes no sense
    if runs.params.flags.flat_bottom
        addcsdye = 0;
    end

    if ~exist('t0', 'var'), t0 = 1; end
    if isempty(t0), [~,~,t0] = runs.locate_resistance; end
    t0 = runs.process_time(t0);
    if ~exist('ntimes', 'var') | isempty(ntimes)
        ntimes = length(runs.time) - t0 + 1;
    end

    ix = runs.spng.sx1:runs.spng.sx2;
    iy = runs.spng.sy1:runs.spng.sy2;

    % read zeta if required
    if strcmpi(name, 'zeta') || addzeta
        runs.read_zeta(t0, ntimes);
        varname = 'zeta';
        titlestr = 'SSH (m)';

        %if strcmpi(name, 'zeta'), addzeta = 0; end
        if ~addzeta, addcsdye = 0; end
    end

    if strcmpi(name, 'pbot')
        if isempty(runs.pbot)
            runs.read_pbot;
        end
        varname = 'pbot';
        titlestr = 'Bottom pressure';
        addcsdye = 0;
    end

    if strcmpi(name, 'betagyre')
        if ~isfield(runs.eddy, 'betagyre')
            runs.betagyre;
        end
        varname = 'eddy.betagyre';
        titlestr = 'Asymmetric streamfn';
        addcsdye = 0;
    end

    if strcmpi(name, 'vor')
        if isempty(runs.vorsurf)
            runs.calc_vorsurf;
        end
        varname = 'vorsurf';
        titlestr = 'Vorticity/f_0';
        addcsdye = 0;
        factor = runs.params.phys.f0;
    end

    % read eddye if required
    if (dyeplot && isempty(runs.eddsurf)) || strcmpi(name, 'eddye') ...
            || strcmpi(name, runs.eddname);
        runs.read_eddsurf(t0, ntimes);
        if (strcmpi(name, 'eddye') || strcmpi(name, runs.eddname)) ...
                && (addcsdye == 0)
            varname = 'eddsurf';
            titlestr = 'Eddy dye';
        end
        if addcsdye == 1
            varname = 'edcsdyesurf';
            titlestr = 'Dyes';
        end
    end

    % read csdye if required
    if (dyeplot && isempty(runs.csdsurf)) || strcmpi(name, 'csdye') ...
            || (strcmpi(name, 'eddye') && addcsdye == 1) ...
            || csdcontourplot || strcmpi(name, runs.csdname)
        runs.read_csdsurf(t0, ntimes);
        if strcmpi(name, 'csdye') || strcmpi(name, runs.csdname)
            varname = 'csdsurf';
            titlestr = 'Cross-shelf dye';
            factor = 1000;
            addcsdye = 0;
        end
    end

    % read rhosurf if required
    if strcmpi(name, 'rho')
        runs.read_rhosurf(t0, ntimes);
        varname = 'rhosurf';
        titlestr = 'Surface \rho';
    end

    if strcmpi(name, 'pv')
        runs.read_pvsurf(t0, ntimes);
        varname = 'pvsurf';
        titlestr = 'Surface PV';
    end

    % bottom velocities
    if strcmpi(name, 'ubot')
        if isempty(runs.ubot), runs.read_velbot; end
        varname = 'ubot';
        titlestr = 'Bottom u';
    end
    if strcmpi(name, 'vbot')
        if isempty(runs.vbot), runs.read_velbot; end
        varname = 'vbot';
        titlestr = 'Bottom v';
    end

    % surface velocities
    if strcmpi(name, 'u')
        varname = 'usurf';
        runs.read_velsurf(t0, ntimes);
    end
    if strcmpi(name, 'v')
        varname = 'vsurf';
        runs.read_velsurf(t0, ntimes);
    end

    if addvelquiver
        runs.read_velsurf(t0, ntimes);
    end

    if ~strcmpi(name, 'eddye')
        addcsdye = 0;
    end

    if isempty(titlestr)
        titlestr = ['Surface ' name];
    end

    % which subplot do I need?
    if (~isempty(runs.csflux) && csfluxplot ~= 0) || ...
            enfluxplot || vecplot
        subplots_flag = 'x';
    else
        if (~isempty(runs.asflux) && asfluxplot == 1)
            subplots_flag = 'y';
        else
            if (~isempty(runs.asflux) && asfluxplot == 2)
                subplots_flag = 'x';
            else
                subplots_flag = [];
            end
        end
    end

    % if provided with axis handle, then assume no subplots
    if ~isempty(hax) & length(hax) == 1
        subplots_flag = [];
    end

    % process addcsdye
    if strcmpi(name, 'eddye') && addcsdye
        tind = t0:t0+ntimes-1;
        if isempty(runs.csflux)
            csfluxThresh = runs.bathy.xsb + 10e3;
        else
            csfluxThresh = runs.csflux.x(csfluxIsobath);
        end

        runs.edcsdyesurf(:,:,tind) = (runs.csdsurf(:,:,tind) < csfluxThresh)*-1;
        runs.edcsdyesurf(:,:,tind) = runs.edcsdyesurf(:,:,tind) + ...
            runs.eddsurf(:,:,tind) .* (runs.edcsdyesurf(:,:,tind) == 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actually plot
    if isempty(hax)
        handles.hfig = figure;
        insertAnnotation([runs.name '.animate_field(' name ')']);
        if subplots_flag == 'x'
            handles.subax(1) = subplot(5,1,[1 2 3]);
            handles.subax(2) = subplot(5,1,[4 5]);
        else
            if subplots_flag == 'y'
                handles.subax(1) = subplot(1,3,[1 2]);
                handles.subax(2) = subplot(1,3,3);
            else
                % no subplots
                handles.subax(1) = gca;
            end
        end
    else
        handles.hfig = gcf;
        if length(hax) == 2
            handles.subax(1) = hax(1);
            handles.subax(2) = hax(2);
        else
            handles.subax(1) = hax;
            axes(hax);
        end
    end

    ii=t0;

    % plot field
    axes(handles.subax(1));
    handles.hfield = runs.plot_surf(varname, 'pcolor', ii);
    handles.hfield.YData = handles.hfield.YData - dy;
    handles.hfield.XData = handles.hfield.XData - dx;
    handles.hfield.CData = handles.hfield.CData / factor;
    hold on;
    handles.hcb = colorbar;
    if AnimateZoom
        if isempty(ZoomXLimStart)
            ZoomXLimStart = xlim;
        end
        if isempty(ZoomXLimEnd)
            ZoomXLimEnd = xlim;
        end
        if isempty(ZoomYLimStart)
            ZoomYLimStart = ylim;
        end
        if isempty(ZoomYLimEnd)
            ZoomYLimEnd = ylim;
        end
    end
    if ~strcmpi(name, 'csdye') && ~strcmpi(name, 'rho') ...
            && ~strcmpi(name, 'dye_01')
        center_colorbar;
    end
    if strcmpi(name, 'pv')
        caxis([min(min(runs.pvsurf(:,:,1))) ...
               max(max(runs.pvsurf(:,:,1)))]);
    end
    if strcmpi(varname, 'eddsurf')
        colormap(runs.eddyeColormap);
        caxis([0 1]);
    end
    if strcmpi(varname, 'rhosurf')
        colormap(flip(colormap));
    end
    clim = caxis;

    % clim = [-1 1]*3e-3;

    if addvelquiver
        % get on interior RHO points
        u = avg1(runs.usurf(:,2:end-1,:),1);
        v = avg1(runs.vsurf(2:end-1,:,:),2);

        if normquiver
            % make all arrows equal length.
            V = hypot(u,v);
            u = u./V;
            v = v./V;

            uref = 1; vref = 1; scale = 1;
        end

        rangex = runs.spng.sx1:dxi:runs.spng.sx2;
        rangey = runs.spng.sy1:dyi:runs.spng.sy2;

        handles.hquiv = quiver(runs.eddy.xr(rangex, rangey)/1000 - dx, ...
                               runs.eddy.yr(rangex, rangey)/1000 - dy, ...
                               u(rangex, rangey, ii)./uref, ...
                               v(rangex, rangey, ii)./vref, scale, ...
                               'Color', quivercolor, 'LineWidth', 1);
    end

    if csdcontourplot
        handles.hcsd = runs.plot_surf('csdsurf', 'contour', ii);
        handles.hcsd.LevelList = csdcontours;
        handles.hcsd.Color = runs.shelfSlopeColor('dark');
        handles.hcsd.LineWidth = 1;
        handles.hcsd.XData = handles.hcsd.XData - dx;
        handles.hcsd.YData = handles.hcsd.YData - dy;
    end

    if strcmpi(name, 'ubot') || strcmpi(name, 'vbot')
        eval(['cmax = max(abs(runs.' varname '(:)));']);
        caxis([-1 1].*cmax);
    end

    % bathy
    handles.hbathy = runs.plot_bathy('contour', bathycolor);
    for bb=1:3
        handles.hbathy{bb}.XData = handles.hbathy{bb}.XData - dx;
        handles.hbathy{bb}.YData = handles.hbathy{bb}.YData - dy;
    end
    if csdcontourplot && modify_bathy
        if runs.bathy.axis == 'y'
            ax = runs.rgrid.y_rho(:,1) - dy;
            hvec = runs.bathy.h(1,:);
        else
            ax = runs.rgrid.x_rho(1,:) - dx;
            hvec = runs.bathy.h(:,1);
        end
        ind = vecfind(ax, csdcontours);
        handles.hbathy{1}.LevelList = round(hvec(ind));
    end

    % plot track
    if drawtrack
        handles.htrack = plot(runs.eddy.mx/1000 - dx, runs.eddy.my/1000 - dy, ...
                              'Color', [152 78 163]/255);
        if drawedgetrack
            if runs.bathy.axis == 'y'
                handles.hedgetrack(1) = plot(runs.eddy.mx/1000 - dx, ...
                                             runs.eddy.vor.ne/1000 - dy);
                handles.hedgetrack(2) = plot(runs.eddy.mx/1000 - dx, ...
                                             runs.eddy.vor.se/1000 - dy);
            else
                handles.hedgetrack(1) = plot(runs.eddy.vor.ee/1000 - dx, ...
                                             runs.eddy.my/1000 - dy);
                handles.hedgetrack(2) = plot(runs.eddy.vor.we/1000 - dx, ...
                                             runs.eddy.my/1000 - dy);
            end
        end
    end

    if drawcenter
        handles.hcen = plot(runs.eddy.mx(ii)/1000 - dx, runs.eddy.my(ii)/1000 - dy, 'kx');
    end

    % zeta too?
    if addzeta & (ii ~= stopzeta);
        handles.hzeta = runs.plot_surf('zeta','contour', ii);
        handles.hzeta.LevelList(handles.hzeta.LevelList < 0) = [];
        handles.hzeta.LevelList = handles.hzeta.LevelList([1:2:end]);
        handles.hzeta.LineWidth = 1;
        handles.hzeta.XData = handles.hzeta.XData - dx;
        handles.hzeta.YData = handles.hzeta.YData - dy;
        set(handles.hzeta, 'Color', 'k', 'LineWidth', 1);

        handles.hzetaneg = runs.plot_surf('zeta','contour', ii);
        handles.hzetaneg.LevelList(handles.hzetaneg.LevelList > 0) = [];
        handles.hzetaneg.LineWidth = 1;
        handles.hzetaneg.XData = handles.hzeta.XData - dx;
        handles.hzetaneg.YData = handles.hzeta.YData - dy;
        handles.hzetaneg.Color = 'k';
        handles.hzetaneg.LineWidth = 1;
        handles.hzetaneg.LineStyle = '--';
    end

    % plot eddy contours
    if vorcontourplot
        handles.he = runs.plot_eddy_contour('contour',ii);
        handles.he.XData = handles.he.XData - dx;
        handles.he.YData = handles.he.YData - dy;
        handles.he.LineWidth = 1;
    end
    if sshplot
        handles.hssh = runs.plot_eddy_sshcontour('contour',ii);
        handles.hssh.XData = handles.hssh.XData - dx;
        handles.hssh.YData = handles.hssh.YData - dy;
        handles.hssh.LineWidth = 1;
    end
    if rhocontourplot
        handles.hrho = runs.plot_rho_contour('contour', ii);
        handles.hrho.XData = handles.hrho.XData - dx;
        handles.hrho.YData = handles.hrho.YData - dy;
        handles.hrho.LineWidth = 1;
    end

    % telescoping lines
    if runs.params.flags.telescoping && telesplot
        linex([runs.params.grid.ixn runs.params.grid.ixp], 'telescope','w');
        liney([runs.params.grid.iyp],'telescope','w');
    end

    % mark point
    if pointplot
        hpt = plot(px(ii) - dx, py(ii) - dy, 'kx', 'MarkerSize', 22);
    end

    % misc stuff
    handles.htitle = runs.set_title(titlestr,ii);
    xlabel('X (km)');ylabel('Y (km)');

    axis image;

    if strcmpi(name, 'zeta')
        if dyeplot
            [~,handles.hdye] = contour(runs.rgrid.x_rho(ix,iy)/1000 - dx, ...
                                runs.rgrid.y_rho(ix,iy)/1000 - dy, ...
                                runs.eddsurf(ix,iy,ii)', ...
                                [1 1]*0.95, 'LineWidth', ...
                                2, 'Color', 'r');
        end
    end

    % draw angle
    %L = createLine(runs.eddy.vor.cx(ii)/1000, runs.eddy.vor.cy(ii)/1000, ...
    %               1, -1*runs.eddy.vor.angle(ii)*pi/180);
    %       hline = drawLine(L);
    if vecplot
        linex(locx - dx); liney(locy - dy);
    end
    caxis(clim); % restore colorbar limits
    maximize(gcf); % pause(0.2);
    beautify(fontsize);
    set(gcf, 'renderer', 'painters');
    handles.htlabel = runs.add_timelabel(ii);

    if addcsdye
        handles.hshelflabel = text(0.1,0.1,'shelf-slope water', ...
                                   'Units', 'normalized', 'Color', 'w');
        caxis([-1 1]*1.5);
        colorbar('hide');
    end

    if nocolorbar
        colorbar('hide');
    end

    if AnimateZoom
        ZoomDeltaT = ceil((ZoomEnd - ZoomStart + 1)/dt);
        ZoomDeltaXmin = (ZoomXLimEnd(1) - ZoomXLimStart(1));
        ZoomDeltaYmin = (ZoomYLimEnd(1) - ZoomYLimStart(1));
        ZoomDeltaXmax = (ZoomXLimEnd(2) - ZoomXLimStart(2));
        ZoomDeltaYmax = (ZoomYLimEnd(2) - ZoomYLimStart(2));

        if ii >= ZoomStart
            xlim(ZoomXLimStart);
            ylim(ZoomYLimStart);
        end
    else
        if ~isempty(limx), xlim(limx); end
        if ~isempty(limy), ylim(limy); end
    end

    % second plot
    if subplots_flag == 'x'
        axes(handles.subax(2));
        if vecplot
            hvec = plot(tvec, vec);
            htime = linex(tvec(ii));
            %ylim([-1 1] .* max(vec)/3);
            % if vector goes through zero, mark zero
            if min(vec(:)) .* max(vec(:)) < 0
                liney(0);
            end
            if csfluxplot == 2
                handles.hleg2 = legend(hvec, isobathNames, 'Location', 'NorthWest');
                handles.hleg2.Position(1) = 0.15;
                handles.hleg2.Position(2) = 0.29;
                if csfluxFinalize
                    handles.hleg2.delete;
                end
            end
            ylabel(laby);
            xlabel('Time (days)');
            if csfluxplot == 2 & csfluxFinalize
                handles.subax(2).YAxis.Exponent = 3;
                ylim([0 max(ylim)]);
            end
        end

        if asfluxplot == 2
            filter = '';

            eval(['left = runs.asflux.' filter 'ikeflux(:,asindex(1)) + '...
                  'runs.asflux.' filter 'ipeflux(:,asindex(1))']);
            eval(['right = runs.asflux.' filter 'ikeflux(:,asindex(2)) + ' ...
                  'runs.asflux.' filter 'ipeflux(:,asindex(2))']);

            total = right - left;

            tvec = runs.asflux.time/86400;

            hold all
            plot(tvec, left);
            plot(tvec, right);
            plot(tvec, total, 'k');
            hleg = legend('west', 'east', 'total (into domain > 0)', ...
                          'Location' ,'NorthWest');
            set(hleg, 'box', 'off');
            liney(0);
            htime = linex(tvec(ii));

            axes(ax);
            linex(runs.asflux.x(asindex)/1000 - dx);
        end

        if csfluxplot == 1
            slopext = runs.csflux.slopext(ix-1, :, csfluxIsobath, csfluxIsobath);
            hflux = plot(runs.rgrid.xr(ix,1)/1000 - dx, ...
                         slopext(:, ii));
            ylim([-1 1] *max(slopext(:)));
            if ~nocolorbar
                oldpos = handles.subax(1).Position;
                newpos = handles.subax(2).Position;
                newpos(3) = oldpos(3);
                handles.subax(2).Position = newpos;
            end
            linkaxes(handles.subax, 'x');
            linkprop(handles.subax, 'XTick');
            hee = linex(runs.eddy.mx(ii)/1000 - dx);
            liney(0);
            xlabel('X (km)');
            handles.htext2 = text(0.05, 0.15, 'Depth Integrated Transport (m^2/s)', ...
                                  'Units', 'normalized');
        end

        beautify(fontsize); % - slows everything down for some reason
        axes(handles.subax(1)); % bring xlabel up
    end
    if subplots_flag == 'y'
        handles.subax(2) = subplot(1,3,3);

        % AS fluxes
        if asfluxplot == 1
            cla
            matrix1 = runs.asflux.ipefluxyt(:,:,asindex(1)) + ...
                      runs.asflux.ikefluxyt(:,:,asindex(1));
            yl = size(matrix1,1);
            isb = runs.bathy.isb;
            hflux1 = plot(matrix1(:,ii), ...
                          runs.rgrid.yr(1,isb:isb+yl-1)/1000 - dy);
            hleg = addlegend(hflux1, ...
                             ['x = ' ...
                              num2str(runs.asflux.x(asindex(1))/1000) ...
                              ' km']);
            if length(asindex) > 1
                matrix2 = runs.asflux.ipefluxyt(:,:,asindex(2)) ...
                          + runs.asflux.ikefluxyt(:,:, asindex(2));
                hold all
                hflux2 = plot(matrix2(:,ii), ...
                              runs.rgrid.yr(1,isb:isb+yl-1)/1000 - dy);
                hleg = addlegend(hflux2, ...
                                 [' x = ' ...
                                  num2str(runs.asflux.x(asindex(2))/1000) ...
                                  ' km']);
            else
                matrix2 = [];
            end

            set(hleg, 'Location', 'NorthWest', 'Box' ,'off');
            xlim([-1 1]*max(abs([matrix1(:); matrix2(:)])));
            liney([runs.bathy.xsb runs.bathy.xsl]/1000 - dy, [], 'k');
            xlabel('Along-isobath depth-integrated energy flux');
            ylabel('Y(km)');
            linex(0); beautify(fontsize);

            axes(ax)
            linex(runs.asflux.x(asindex)/1000 - dx);

            linkaxes(ax, 'y');
        end

        beautify(fontsize);
    end

    eval(commands);
    if ntimes > 1
        runs.video_update();
        for ii = t0+1:dt:t0+ntimes*dt-1
            %L = createLine(runs.eddy.vor.cx(ii)/1000, runs.eddy.vor.cy(ii)/1000, ...
            %           1, -1*runs.eddy.vor.angle(ii)*pi/180);
            %delete(hline);
            if ~isempty(runs.csflux) && csfluxplot > 0
                %set(hee_zeta, 'XData', [1 1]* runs.eddy.vor.ee(ii)/ ...
                %              1000);
                %axes(ax);
                %hline = drawLine(L);
            end

            runs.update_surf(varname, handles.hfield,ii);
            handles.hfield.CData = handles.hfield.CData / factor;
            if vorcontourplot
                runs.update_eddy_contour(he,ii);
            end
            if addvelquiver
                set(handles.hquiv,'UData',u(rangex, rangey, ii)/uref);
                set(handles.hquiv,'VData',v(rangex, rangey, ii)/vref);
            end

            if addzeta & (ii >= stopzeta)
                handles.hzeta.delete;
                addzeta = 0;
            end

            if addzeta
                runs.update_surf('zeta', handles.hzeta, ii);
                runs.update_surf('zeta', handles.hzetaneg, ii);
                if any(handles.hzeta.LevelList < 0)
                    handles.hzeta.LevelList(handles.hzeta.LevelList < 0) = [];
                end
                if any(handles.hzetaneg.LevelList > 0)
                    handles.hzetaneg.LevelList(handles.hzeta.LevelList > 0) = [];
                end
            end

            if csdcontourplot
                runs.update_surf('csdsurf', handles.hcsd, ii);
            end

            if pointplot
                hpt.XData = px(ii) - dx;
                hpt.YData = py(ii) - dy;
            end

            if drawcenter
                handles.hcen.XData = runs.eddy.mx(ii)/1000 - dx;
                handles.hcen.YData = runs.eddy.my(ii)/1000 - dy;
            end

            if AnimateZoom & ii>=ZoomStart & ii<=ZoomEnd
                olimx = handles.subax(1).XLim;
                olimy = handles.subax(1).YLim;
                zlimx(1) = olimx(1) + ZoomDeltaXmin/ZoomDeltaT;
                zlimx(2) = olimx(2) + ZoomDeltaXmax/ZoomDeltaT;
                zlimy(1) = olimy(1) + ZoomDeltaYmin/ZoomDeltaT;
                zlimy(2) = olimy(2) + ZoomDeltaYmax/ZoomDeltaT;

                handles.subax(1).XLim = zlimx;
                handles.subax(1).YLim = zlimy;
            end

            % ssh contour
            if sshplot
                runs.update_eddy_sshcontour(handles.hssh,ii);
            end
            if rhocontourplot
                runs.update_rho_contour(handles.hrho,ii);
            end
            runs.update_title(handles.htitle,titlestr,ii);

            % eddye contours
            if dyeplot
                set(handles.hdye, 'ZData', runs.eddsurf(ix,iy,ii)');
            end

            % slopewater flux plots
            if ~isempty(runs.csflux) && csfluxplot == 1
                axis(handles.subax(2));
                if exist('htime', 'var')
                    set(htime, 'XData', [1 1]*runs.csflux.time(ii)/ ...
                               86400);
                else
                    set(hflux, 'YData', slopext(:,ii));
                    set(hee, 'XData', [1 1]*runs.eddy.vor.ee(ii)/1000);
                end
            end

            % AS eddy water flux plots
            if ~isempty(runs.asflux) && asfluxplot == 1
                axis(handles.subax(2));
                if exist('htime', 'var')
                    set(htime, 'XData', [1 1]*runs.asflux.time(ii)/ ...
                               86400);
                else
                    set(hflux1, 'XData', matrix1(:,ii));
                    if length(asindex) > 1
                        set(hflux2, 'XData', matrix2(:,ii));
                    end
                end
            end

            % mark time for time-series plots
            if vecplot || (asfluxplot == 2)
                set(htime, 'XData', [1 1]*tvec(ii));
            end

            % update time label
            runs.update_timelabel(handles.htlabel, ii);
            runs.video_update();
            pause(1);
        end
        runs.video_write();
    end
end