% this is going to be a runArray version of plot_eddye.
% given some non-dimensional time, allow me to compare sections of
% quantities through the center of the eddies
% each figure shows multiple timesteps for a given field _for a
% given run_.
function [hh] = plot_sections(runArray, varname, ndtimes)

    if isempty(runArray.filter)
        runArray.filter = 1:length(runArray.array);
    end

    nruns = length(runArray.filter);
    if exist('ndtimes', 'var')
        nt = length(ndtimes);
    else
        nt = 1;
    end

    cmap = flipud(cbrewer('div','RdYlBu', 32));
    %cmap = cbrewer('seq','Reds',32);

    % multiple runs on same figure? set to 1
    % multiple timesteps on same figure? set to 0
    runs_on_one_fig = 1;

    for tt = 1:nt
        if runs_on_one_fig
            figure;
        end

        for ii = 1:nruns
            if ~runs_on_one_fig
                figure;
            end

            ff = runArray.filter(ii);
            run = runArray.array(ff);
            name = runArray.getname(ff);

            if exist('ndtimes', 'var')
                ndtime = run.ndtime;
                % find required tindices
                tind = find_approx(ndtime, ndtimes(tt));
            else
                [~,~,tind] = run.locate_resistance;
                ndtimes = run.ndtime(tind);
            end

            yscale = run.eddy.vor.dia(1)/2
            zscale = run.eddy.Lgauss(1)
            y0 = run.eddy.my(tind); %

            if strcmpi(varname, 'v')
                ymat = repmat(run.rgrid.y_v(:,1), [1 run.rgrid.N]);
                zmat = run.rgrid.z_v(:,:,1)';
                xloc = num2str(run.eddy.mx(tind) - ...
                               run.eddy.vor.dia(tind)/2);
            else
                ymat = repmat(run.rgrid.y_rho(:,1), [1 run.rgrid.N]);
                zmat = run.rgrid.z_r(:,:,1)';
                xloc = num2str(run.eddy.mx(tind));
            end

            ynorm = (ymat - y0)./yscale;

            axind = sub2ind([nruns nt], ii,tt);
            if runs_on_one_fig
                ax(axind) = subplot(1,nruns,ii);
            else
                ax(axind) = subplot(1,nt,tt);
            end

            vol = {'x' xloc xloc};

            % read in variable
            if strcmpi(varname, 'rhoanom')
                rback = dc_roms_read_data(run.dir, 'rho', 1, ...
                                        {'x' 1 1}, [], run.rgrid, 'his');
                var = dc_roms_read_data(run.dir, 'rho', tind, ...
                                        vol, [], run.rgrid, 'his');
                var = bsxfun(@minus, var, rback);
            else
                var = dc_roms_read_data(run.dir, varname, tind, ...
                                        vol, [], run.rgrid, 'his');
            end

            % plot variable
            [c0,h0] = contourf(ynorm, zmat./zscale, var, 40);
            hh(ii) = h0;
            hh(ii).LevelList = hh(1).LevelList;

            liney(round(-1*run.eddy.Lgauss(tind)/zscale, 2, ...
                        'significant'));
            linex(0);

            if ii == 1
                clim = caxis;
            else
                caxis(clim);
            end

            colormap(cmap); shading interp
            if ~strcmpi(varname, 'rho')
                center_colorbar;
            end

            % overlay profile of initial eddy
            limy = xlim .* yscale; limz = ylim .* zscale;
            [ymat2,zmat2] = ndgrid(linspace(limy(1),limy(2),100), ...
                                     linspace(limz(1),limz(2),1000));
            % initial profile
            prof = exp(-((ymat2-run.eddy.my(tind))./ ...
                         run.eddy.vor.dia(1)*2).^2) .* ...
                   exp(-(zmat2./run.eddy.Lgauss(1)).^2);
            % current profile
            prof2 = exp(-((ymat2-run.eddy.my(tind))./ ...
                         run.eddy.vor.dia(tt)*2).^2) .* ...
                   exp(-(zmat2./run.eddy.Lgauss(tind)).^2);
            hold on

            % patch bathymetry
            patch(([run.rgrid.y_rho(:,1); min(run.rgrid.y_rho(:,1))] ...
                   - y0)./yscale, ...
                  -1*[min(run.rgrid.h(:)); run.rgrid.h(:,1)]./zscale, 'k');

            % initial profile
            %contour((ymat2-run.bathy.xsb)./yscale, ...
            %        zmat2./zscale, prof, [1 1]*0.5, 'color', 'b', 'LineWidth', ...
            %        2);
            % current profile
            %contour((ymat2-run.bathy.xsb)./yscale, ...
            %        zmat2./zscale, prof2, [1 1]*0.5, 'color', 'b', 'LineWidth', ...
            %        2);

            % line at 1 deformation radius from shelfbreak
            %linex(1.2);

            beautify([18 18 20]); box on;

            if runs_on_one_fig
                %title(['S_\alpha = ', num2str(run.bathy.S_sl)]);
                title(runArray.name{ff});
                set(gca, 'TickDir', 'out');
                set(gca, 'box', 'on');
                set(gcf, 'renderer', 'zbuffer');
                if ii ~= 1
                    set(gca, 'YTickLabel', []);
                end
                %ylim([-round(run.params.eddy.depth/100)*100 0]);
            else
                title(['ndtime = ', num2str(ndtimes(tt))]);
            end

            if axind > 1
                ax(axind).XTick = ax(1).XTick;
                ax(axind).YTick = ax(1).YTick;
            end
        end
        linkaxes(ax, 'xy');

        supAxes = [.1 .1 .85 .85];
        if runs_on_one_fig
            %            suplabel('t', ['ndtime = ' num2str(ndtime)], supAxes);
        else
            % suplabel('t', run.name);
        end
        if runs_on_one_fig
            [~, hx] = suplabel(['Distance from shelfbreak / Deformation ' ...
                                'radius'], 'x', supAxes);
            [~, hy] = suplabel('Depth / Initial vertical scale', 'y', supAxes);
            set(hx, 'fontSize', 18);
            set(hy, 'fontSize', 18);
        end
        title(['ndtime = ', num2str(ndtimes(tt))]);
    end
end