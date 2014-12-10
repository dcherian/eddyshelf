% this is going to be a runArray version of plot_eddye.
% given some non-dimensional time, allow me to compare sections of
% quantities through the center of the eddies
% each figure shows multiple timesteps for a given field _for a
% given run_.
function [] = plot_sections(runArray, varname, ndtimes)

    if isempty(runArray.filter)
        runArray.filter = 1:length(runArray.array);
    end

    nruns = length(runArray.filter);
    nt = length(ndtimes);

    cmap = cbrewer('seq','Reds',32);

    % multiple runs on same figure? set to 1
    % multiple timesteps on same figure? set to 0
    runs_on_one_fig = 1;

    for tt = 1:nt
        figure;
        for ii = 1:nruns
            ff = runArray.filter(ii);
            run = runArray.array(ff);
            name = runArray.getname(ff);

            ndtime = run.time ./ run.tscale;

            % find required tindices
            tind = vecfind(ndtime, ndtimes);

            yscale = run.rrdeep;
            zscale = 1;
            ymat = repmat(run.rgrid.y_rho(:,1), [1 run.rgrid.N]);
            zmat = run.rgrid.z_r(:,:,1)';

            axind = sub2ind([nruns nt], ii,tt);
            if runs_on_one_fig
                ax(axind) = subplot(1,nruns,ii);
            else
                ax(axind) = subplot(1,nt,tt);
            end

            % read in variable
            var = dc_roms_read_data(run.dir, varname, tind(tt), ...
                                   {'x' num2str(run.eddy.cx(tind(tt))) ...
                                num2str(run.eddy.cx(tind(tt)))}, [], ...
                                   run.rgrid, 'his');

            % plot variable
            pcolorcen((ymat-run.bathy.xsb)./run.rrdeep, ...
                      zmat./zscale, var);
            liney(-1*run.eddy.Lgauss(tind(tt))/zscale);
            colormap(cmap);

            % overlay profile of initial eddy
            limy = xlim .* yscale; limz = ylim .* zscale;
            [ymat2,zmat2] = ndgrid(linspace(limy(1),limy(2),100), ...
                                     linspace(limz(1),limz(2),1000));
            % initial profile
            prof = exp(-((ymat2-run.eddy.cy(tind(tt)))./ ...
                         run.eddy.vor.dia(1)*2).^2) .* ...
                   exp(-(zmat2./run.eddy.Lgauss(1)).^2);
            % current profile
            prof2 = exp(-((ymat2-run.eddy.cy(tind(tt)))./ ...
                         run.eddy.vor.dia(tt)*2).^2) .* ...
                   exp(-(zmat2./run.eddy.Lgauss(tind(tt))).^2);
            hold on
            % patch bathymetry
            patch(([run.rgrid.y_rho(:,1); min(run.rgrid.y_rho(:,1))] ...
                   - run.bathy.xsb)./yscale, ...
                  -1*[min(run.rgrid.h(:)); run.rgrid.h(:,1)], 'k');

            % initial profile
            contour((ymat2-run.bathy.xsb)./yscale, ...
                    zmat2./zscale, prof, [1 1]*0.5, 'color', 'b', 'LineWidth', ...
                    2);
            % current profile
            contour((ymat2-run.bathy.xsb)./yscale, ...
                    zmat2./zscale, prof2, [1 1]*0.5, 'color', 'b', 'LineWidth', ...
                    2);

            % line at 1 deformation radius from shelfbreak
            linex(1.2);

            beautify([18 18 20]);

            if runs_on_one_fig
                title(['S_\alpha = ', num2str(run.bathy.S_sl)]);
                set(gca, 'TickDir', 'out');
                set(gca, 'box', 'on');
                set(gcf, 'renderer', 'zbuffer');
                if ii ~= 1
                    set(gca, 'YTickLabel', []);
                end
                ylim([-round(run.params.eddy.depth/100)*100 0]);
            else
                title(['ndtime = ', num2str(ndtimes(tt))]);
            end

        end
        linkaxes(ax, 'xy');

        supAxes = [.1 .1 .85 .85];
        if runs_on_one_fig
            %            suplabel('t', ['ndtime = ' num2str(ndtime)], supAxes);
        else
            suplabel('t', run.name);
        end
        [~, hx] = suplabel('Distance from shelfbreak (km)', 'x', supAxes);
        [~, hy] = suplabel('Depth (m)', 'y', supAxes);
        set(hx, 'fontSize', 18);
        set(hy, 'fontSize', 18);
    end
end