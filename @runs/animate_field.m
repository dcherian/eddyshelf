%  animate 2D field at the surface presumably.
function [] = animate_field(runs, name, t0, ntimes)
    runs.video_init(name);

    titlestr = [];

    dt = 10;

    csfluxplot = 0; % 0 = no flux plot
                    % 1 = instantaneous x-profile;
    asfluxplot = 0; % 0 = no flux plot
                    % 1 = instantaneous y-profile;
                    % 2 = time series plot : left, right and total
    if (asfluxplot == 1) || (asfluxplot == 2) % which location for asflux plot?
        asindex = [1 2];
    end

    enfluxplot = 0; % plot AS energy flux ?
    sshplot = 0; % plot ssh-contour too?
    dyeplot = 0; % plot eddye contour too?
    telesplot = 0;  % plot lines where grid stretching starts
                    % and ends

    vecplot = 1; % plot some time vector (assign tvec and vec);
    if vecplot
        %%% integrated energy asflux
        %tvec = runs.eddy.t;
        %vec = (runs.asflux.ikeflux(:,3) + runs.asflux.ipeflux(:,3) ...
        %      - runs.asflux.ikeflux(:,2) - runs.asflux.ipeflux(:,2));
        %laby = 'Integrated energy flux';
        %locx = runs.asflux.x(2:3)/1000; locy = [];

        %%% dE/dt
        tvec = avg1(runs.eddy.t);
        vec = smooth(diff(runs.eddy.KE + runs.eddy.PE)./ ...
                     diff(runs.eddy.t*86400), 4);
        laby = 'dE/dt';
        locx = []; locy = [];

        %%% asflux time series
        %tvec = runs.asflux.time/runs.tscale;
        %vec = runs.asflux.ikeflux(:,2);
        %locx = runs.asflux.x(2); locy = [];

        %%% cross-sb shelf water flux
        %tvec = runs.csflux.time/86400;
        %vec = runs.csflux.shelf.west.shelf;
        %laby = 'Shelf water flux (m^3/s)';
        %locy = runs.bathy.xsb/1000; locx = [];

        %%% area plot
        %vec = runs.eddy.vor.lmin .* runs.eddy.vor.lmaj;
        %tvec = runs.eddy.t;
        %laby = 'Surface area (m^2)';
        %locx = []; locy = [];
    end

    if ~exist('ntimes', 'var'), ntimes = length(runs.time); end
    if ~exist('t0', 'var'), t0 = 1; end

    % read zeta if required
    if strcmpi(name, 'zeta')
        runs.read_zeta(t0, ntimes);
        varname = 'zeta';
        titlestr = 'SSH (m)';
    end

    % read eddye if required
    if (dyeplot && isempty(runs.eddsurf)) || strcmpi(name, 'eddye')
        runs.read_eddsurf(t0, ntimes);
        if strcmpi(name, 'eddye')
            varname = 'eddsurf';
            titlestr = 'Eddy dye';
        end
    end

    % read rhosurf if required
    if strcmpi(name, 'rho');
        runs.read_rhosurf(t0, ntimes);
        varname = 'rhosurf';
        titlestr = 'Surface \rho';
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

    %%%% actually plot
    figure; insertAnnotation([runs.name '.animate_field']);
    if subplots_flag == 'x'
        ax = subplot(3,1,[1 2]);
    else
        if subplots_flag == 'y'
            ax = subplot(1,3,[1 2]);
        end
    end
    ii=t0;

    % plot field
    hz = runs.plot_surf(varname, 'pcolor', ii);

    hold on;
    colorbar; freezeColors;

    % bathy
    hbathy = runs.plot_bathy('contour','k');

    % plot track
    plot(runs.eddy.mx/1000, runs.eddy.my/1000);
    plot(runs.eddy.mx/1000, runs.eddy.vor.ne/1000);
    plot(runs.eddy.mx/1000, runs.eddy.vor.se/1000);

    % plot eddy contours
    he = runs.plot_eddy_contour('contour',ii);
    if sshplot
        he2 = runs.plot_eddy_sshcontour('contour',ii);
    end

    % telescoping lines
    if runs.params.flags.telescoping && telesplot
        linex([runs.params.grid.ixn runs.params.grid.ixp], 'telescope','w');
        liney([runs.params.grid.iyp],'telescope','w');
    end

    % misc stuff
    ht = runs.set_title(titlestr,ii);
    xlabel('X (km)');ylabel('Y (km)');
    if isempty(subplots_flag)
        axis image;
    else
        %axis equal;
    end

    if strcmpi(name, 'zeta')
        if dyeplot
            [~,hedd2] = contour(runs.rgrid.x_rho/1000, runs.rgrid.y_rho/1000, ...
                                runs.eddsurf(:,:,ii)', [1 1]*0.95, 'LineWidth', ...
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
    maximize(gcf); pause(0.2);
    beautify([16 16 18]);
    ax = gca;

    % second plot
    if subplots_flag == 'x'
        ax2 = subplot(3,1,3);
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
            hflux = plot(runs.rgrid.xr(2:end-1,1)/1000, ...
                         runs.csflux.shelfxt(:, ii));
            ylim([min(runs.csflux.shelfxt(:)) ...
                  max(runs.csflux.shelfxt(:))]);
            hee = linex(runs.eddy.vor.ee(ii)/1000);
            oldpos = get(ax, 'Position');
            newpos = get(ax2, 'Position');
            newpos(3) = oldpos(3);
            set(ax2, 'Position', newpos);
            linkaxes([ax ax2], 'x');
            liney(0);
            xlabel('X (km)');
            title('\int v(x,z,t)dz (m^2/s)');
        end
    end
    if subplots_flag == 'y'
        ax2 = subplot(1,3,3);

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

            linkaxes([ax ax2], 'y');
        end

        beautify;
    end

    if ntimes > 1
        runs.video_update();
        for ii = t0+1:dt:ntimes
            %L = createLine(runs.eddy.vor.cx(ii)/1000, runs.eddy.vor.cy(ii)/1000, ...
            %           1, -1*runs.eddy.vor.angle(ii)*pi/180);
            %delete(hline);
            if ~isempty(runs.csflux) && csfluxplot > 0
                set(hee_zeta, 'XData', [1 1]* runs.eddy.vor.ee(ii)/ ...
                              1000);
                axes(ax);
                %hline = drawLine(L);
            end

            runs.update_surf(varname, hz,ii);
            runs.update_eddy_contour(he,ii);
            % ssh contour
            if sshplot
                runs.update_eddy_sshcontour(he2,ii);
            end
            runs.update_title(ht,titlestr,ii);

            % eddye contours
            if dyeplot
                set(hedd2, 'ZData', runs.eddsurf(:,:,ii)');
            end

            % shelfwater flux plots
            if ~isempty(runs.csflux) && csfluxplot == 1
                axis(ax2);
                if exist('htime', 'var')
                    set(htime, 'XData', [1 1]*runs.csflux.time(ii)/ ...
                               86400);
                else
                    set(hflux, 'YData', runs.csflux.shelfxt(:,ii));
                    set(hee, 'XData', [1 1]*runs.eddy.vor.ee(ii)/1000);
                end
            end

            % AS eddy water flux plots
            if ~isempty(runs.asflux) && asfluxplot == 1
                axis(ax2);
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
            runs.video_update();
            pause(1);
        end
        runs.video_write();
    end
end
