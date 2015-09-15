%  animate 2D field at the surface presumably.
function [] = animate_field(runs, name, hax, t0, ntimes)
    if ~exist('hax','var'), hax = []; end

    runs.video_init(name);
    titlestr = [];

    dt = 3;
    factor = 1; % scale variable (1 by default)

    csfluxplot = 0; % 0 = no flux plot
                    % 1 = instantaneous x-profile;
    asfluxplot = 0; % 0 = no flux plot
                    % 1 = instantaneous y-profile;
                    % 2 = time series plot : left, right and total
    if (asfluxplot == 1) || (asfluxplot == 2) % which location for asflux plot?
        asindex = [1 2];
    end

    % time series?
    enfluxplot = 0; % plot AS energy flux ?
    vecplot = 1; % plot some time vector (assign tvec and vec);
    pointplot = 0; % mark some point on the map

    % eddy contours?
    rhocontourplot = 0; % plot eddy drho contour too?
    vorcontourplot = 1; % vorticity contour
    sshplot = 0; % plot ssh-contour too?
    dyeplot = 0; % plot eddye contour too?

    % extra contours
    addcsdye = 0; % add csdye to eddye plot?
    addzeta = 0; % overlay zeta contours

    csdcontourplot = 1; % contour csd contours
    try
        isobath = [3 5 7];
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

    % eddy diagnostics
    drawtrack = 1; % plot eddy track?
    drawedgetrack = 0; % plot eddy's edge tracks?
    drawcenter = 1; % mark eddy center location.

    % grid diagnostics
    telesplot = 0;  % plot lines where grid stretching starts
                    % and ends

    if pointplot
        px = (runs.eddy.mx + runs.eddy.vor.dia(1))/1000;
        py = runs.eddy.my/1000;
    end

    if vecplot
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

            %%% cross-sb shelf water flux
            tvec = runs.csflux.time/86400;
            for kkk=1:length(isobath)
                vec(:,kkk) = runs.csflux.west.slope(:,isobath(kkk), ...
                                                    isobath(kkk));
            end

            laby = 'Slope water flux (m^3/s)';
            locy = runs.bathy.xsb/1000; locx = [];

            %%% area plot
            %vec = runs.eddy.vor.lmin .* runs.eddy.vor.lmaj;
            %tvec = runs.eddy.t;
            %laby = 'Surface area (m^2)';
            %locx = []; locy = [];
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
    if ~exist('ntimes', 'var')
        ntimes = length(runs.time) - t0 + 1;
    end

    ix = runs.spng.sx1:runs.spng.sx2;
    iy = runs.spng.sy1:runs.spng.sy2;

    % read zeta if required
    if strcmpi(name, 'zeta') || addzeta
        runs.read_zeta(t0, ntimes);
        varname = 'zeta';
        titlestr = 'SSH (m)';

        if strcmpi(name, 'zeta'), addzeta = 0; end
        addcsdye = 0;
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
            titlestr = 'Red = eddy water | Blue = shelf/slope water';
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

    if ~strcmpi(name, 'eddye')
        addcsdye = 0;
    end

    if isempty(titlestr)
        titlestr = ['Surface ' name];
    end

    % which subplot do I need?
    if (~isempty(runs.csflux) && csfluxplot == 1) || ...
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
    if ~isempty(hax)
        subplots_flag = [];
    end

    % process addcsdye
    if strcmpi(name, 'eddye') && addcsdye
        tind = t0:t0+ntimes-1;
        runs.edcsdyesurf(:,:,tind) = (runs.csdsurf(:,:,tind) < ...
                                      (runs.bathy.xsb+10e3))*-1;
        runs.edcsdyesurf(:,:,tind) = runs.edcsdyesurf(:,:,tind) + ...
            runs.eddsurf(:,:,tind) .* (runs.edcsdyesurf(:,:,tind) == 0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actually plot
    if isempty(hax)
        figure;
        insertAnnotation([runs.name '.animate_field(' name ')']);
    else
        axis(hax);
    end

    if subplots_flag == 'x'
        ax(1) = subplot(5,1,[1 2 3]);
    else
        if subplots_flag == 'y'
            ax(1) = subplot(1,3,[1 2]);
        end
    end
    ii=t0;

    % plot field
    hz = runs.plot_surf(varname, 'pcolor', ii);
    hz.CData = hz.CData / factor;
    hold on;
    colorbar;
    if ~strcmpi(name, 'csdye') && ~strcmpi(name, 'rho')
        center_colorbar;
    end
    clim = caxis;

    if csdcontourplot
        hcsd = runs.plot_surf('csdsurf', 'contour', ii);
        hcsd.LevelList = csdcontours;
        hcsd.Color = [1 1 1]*0.5;
        hcsd.LineWidth = 2;
    end

    if strcmpi(name, 'ubot') || strcmpi(name, 'vbot')
        eval(['cmax = max(abs(runs.' varname '(:)));']);
        caxis([-1 1].*cmax);
    end

    % bathy
    hbathy = runs.plot_bathy('contour', [1 1 1]*0.7);
    if csdcontourplot && modify_bathy
        if runs.bathy.axis == 'y'
            ax = runs.rgrid.y_rho(:,1);
            hvec = runs.bathy.h(1,:);
        else
            ax = runs.rgrid.x_rho(1,:);
            hvec = runs.bathy.h(:,1);
        end
        ind = vecfind(ax, csdcontours);
        hbathy{1}.LevelList = round(hvec(ind));
    end

    % plot track
    if drawtrack
        plot(runs.eddy.mx/1000, runs.eddy.my/1000);
        if drawedgetrack
            if runs.bathy.axis == 'y'
                plot(runs.eddy.mx/1000, runs.eddy.vor.ne/1000);
                plot(runs.eddy.mx/1000, runs.eddy.vor.se/1000);
            else
                plot(runs.eddy.vor.ee/1000, runs.eddy.my/1000);
                plot(runs.eddy.vor.we/1000, runs.eddy.my/1000);
            end
        end
    end

    if drawcenter
        hcen = plot(runs.eddy.mx(ii)/1000, runs.eddy.my(ii)/1000, 'kx');
    end

    % zeta too?
    if addzeta
        hzeta = runs.plot_surf('zeta','contour', ii);
        set(hzeta, 'Color', 'k', 'LineWidth', 2);
    end

    % plot eddy contours
    if vorcontourplot
        he = runs.plot_eddy_contour('contour',ii);
    end
    if sshplot
        he2 = runs.plot_eddy_sshcontour('contour',ii);
    end
    if rhocontourplot
        he3 = runs.plot_rho_contour('contour', ii);
    end

    % telescoping lines
    if runs.params.flags.telescoping && telesplot
        linex([runs.params.grid.ixn runs.params.grid.ixp], 'telescope','w');
        liney([runs.params.grid.iyp],'telescope','w');
    end

    % mark point
    if pointplot
        hpt = plot(px(ii), py(ii), 'kx', 'MarkerSize', 22);
    end

    % misc stuff
    ht = runs.set_title(titlestr,ii);
    xlabel('X (km)');ylabel('Y (km)');
    axis image;

    if strcmpi(name, 'zeta')
        if dyeplot
            [~,hedd2] = contour(runs.rgrid.x_rho(ix,iy)/1000, ...
                                runs.rgrid.y_rho(ix,iy)/1000, ...
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
        linex(locx); liney(locy);
    end
    caxis(clim); % restore colorbar limits
    maximize(gcf); pause(0.2);
    beautify([16 16 18]);
    ax(1) = gca;
    htext = runs.add_timelabel(ii);

    if addcsdye
        caxis([-1 1]*1.5);
        colorbar('hide');
        correct_ticks('y',[],[2 6]);
    end

    % second plot
    if subplots_flag == 'x'
        ax(2) = subplot(5,1,[4 5]);
        if vecplot
            hvec = plot(tvec, vec);
            htime = linex(tvec(ii));
            %ylim([-1 1] .* max(vec)/3);
            % if vector goes through zero, mark zero
            if min(vec(:)) .* max(vec(:)) < 0
                liney(0);
            end

            ylabel(laby);
            xlabel('Time (days)');
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
            beautify;

            axes(ax);
            linex(runs.asflux.x(asindex)/1000);
        end

        if csfluxplot == 1
            slopext = runs.csflux.slopext(ix-1, :, index);
            hflux = plot(runs.rgrid.xr(ix,1)/1000, ...
                         slopext(:, ii));
            ylim([min(slopext(:)) max(slopext(:))]);
            hee = linex(runs.eddy.vor.ee(ii)/1000);
            oldpos = get(ax(1), 'Position');
            newpos = get(ax(2), 'Position');
            newpos(3) = oldpos(3);
            set(ax(2), 'Position', newpos);
            linkaxes(ax, 'x');
            liney(0);
            xlabel('X (km)');
            title('\int v(x,z,t)dz (m^2/s)');
        end
    end
    if subplots_flag == 'y'
        ax(2) = subplot(1,3,3);

        % AS fluxes
        if asfluxplot == 1
            cla
            matrix1 = runs.asflux.ipefluxyt(:,:,asindex(1)) + ...
                      runs.asflux.ikefluxyt(:,:,asindex(1));
            yl = size(matrix1,1);
            isb = runs.bathy.isb;
            hflux1 = plot(matrix1(:,ii), ...
                          runs.rgrid.yr(1,isb:isb+yl-1)/1000);
            hleg = addlegend(hflux1, ...
                             ['x = ' ...
                              num2str(runs.asflux.x(asindex(1))/1000) ...
                              ' km']);
            if length(asindex) > 1
                matrix2 = runs.asflux.ipefluxyt(:,:,asindex(2)) ...
                          + runs.asflux.ikefluxyt(:,:, asindex(2));
                hold all
                hflux2 = plot(matrix2(:,ii), ...
                              runs.rgrid.yr(1,isb:isb+yl-1)/1000);
                hleg = addlegend(hflux2, ...
                                 [' x = ' ...
                                  num2str(runs.asflux.x(asindex(2))/1000) ...
                                  ' km']);
            else
                matrix2 = [];
            end

            set(hleg, 'Location', 'NorthWest', 'Box' ,'off');
            xlim([-1 1]*max(abs([matrix1(:); matrix2(:)])));
            liney([runs.bathy.xsb runs.bathy.xsl]/1000, [], 'k');
            xlabel('Along-isobath depth-integrated energy flux');
            ylabel('Y(km)');
            linex(0); beautify;

            axes(ax)
            linex(runs.asflux.x(asindex)/1000);

            linkaxes(ax, 'y');
        end

        beautify;
    end

    if ntimes > 1
        runs.video_update();
        for ii = t0+1:dt:t0+ntimes-1
            %L = createLine(runs.eddy.vor.cx(ii)/1000, runs.eddy.vor.cy(ii)/1000, ...
            %           1, -1*runs.eddy.vor.angle(ii)*pi/180);
            %delete(hline);
            if ~isempty(runs.csflux) && csfluxplot > 0
                %set(hee_zeta, 'XData', [1 1]* runs.eddy.vor.ee(ii)/ ...
                %              1000);
                %axes(ax);
                %hline = drawLine(L);
            end

            runs.update_surf(varname, hz,ii);
            hz.CData = hz.CData / factor;
            if vorcontourplot
                runs.update_eddy_contour(he,ii);
            end

            if addzeta
                runs.update_surf('zeta', hzeta, ii);
            end

            if csdcontourplot
                runs.update_surf('csdsurf', hcsd, ii);
            end

            if pointplot
                hpt.XData = px(ii);
                hpt.YData = py(ii);
            end

            if drawcenter
                hcen.XData = runs.eddy.mx(ii)/1000;
                hcen.YData = runs.eddy.my(ii)/1000;
            end

            % ssh contour
            if sshplot
                runs.update_eddy_sshcontour(he2,ii);
            end
            if rhocontourplot
                runs.update_rho_contour(he3,ii);
            end
            runs.update_title(ht,titlestr,ii);

            % eddye contours
            if dyeplot
                set(hedd2, 'ZData', runs.eddsurf(ix,iy,ii)');
            end

            % slopewater flux plots
            if ~isempty(runs.csflux) && csfluxplot == 1
                axis(ax(2));
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
                axis(ax(2));
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
            runs.update_timelabel(htext, ii);
            runs.video_update();
            pause(1);
        end
        runs.video_write();
    end
end