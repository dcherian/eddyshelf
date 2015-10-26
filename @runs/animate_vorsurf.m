function [] = animate_vorsurf(runs, hax, t0, ntimes)
    if isempty(runs.vorsurf)
        runs.calc_vorsurf();
    end

    if ~exist('hax', 'var'), hax = []; end
    if ~exist('t0', 'var'), t0 = 1; end
    if ~exist('ntimes', 'var'), ntimes = length(runs.time); end

    runs.video_init('surfvorcsd');

    ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
    iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

    % read in dye values
    if isempty(runs.csdsurf)
        runs.csdsurf = dc_roms_read_data(runs.dir, runs.csdname, [], {'z' ...
                            runs.rgrid.N runs.rgrid.N}, [], runs.rgrid, ...
                                         'his');
    end

    csdlevels = runs.bathy.xsb + [0 -20] * 1000;

    tt = t0;
    if isempty(hax), figure(); else axis(hax); end
    maximize();
    vorsurf = bsxfun(@rdivide, runs.vorsurf, ...
                     avg1(avg1(runs.rgrid.f', 1), 2));
    if ntimes ~= 1
        vormax = max(abs(vorsurf(:)));
    else
        vormax = max(max(abs(vorsurf(ix,iy,t0))));
    end

    levels = linspace(-vormax,vormax,20);
    [hh] = pcolor(runs.rgrid.xvor(ix,iy)/1000, ...
                  runs.rgrid.yvor(ix,iy)/1000, ...
                  vorsurf(ix,iy,tt));
    caxis([-1 1] * vormax); hcbar = colorbar; shading interp;
    hcbar.Label.String = '\zeta/f';
    hold on
    % bathy contour
    runs.plot_bathy('contour', [1 1 1]*0.7);

    % eddy contours
    hedd = runs.plot_eddy_contour('contour', tt);
    % hssh = runs.plot_eddy_sshcontour('contour', tt);

    % dye-contours (current)
    %[~,hcsd] = contour(runs.rgrid.x_rho'/1000, runs.rgrid.y_rho'/1000, ...
    %                   runs.csdsurf(:,:,tt), csdlevels, 'Color', [1 1 1]*0.5, ...
    %                   'LineWidth', 2);
    % dye contours (initial)
    % contour(runs.rgrid.x_rho'/1000, runs.rgrid.y_rho'/1000, ...
    %        runs.csdsurf(:,:,tt), csdlevels(1), 'Color', [1 1 1]*0.9, ...
    %        'LineWidth', 2)
    xlabel('X (km)'); ylabel('Y (km)');
    axis image;
    ht = title(['Surface vorticity @ t = ' ...
                num2str(runs.time(tt)/86400) ' days']);
    ax = gca;
    runs.add_timelabel(tt);
    colormap(flipud(cbrewer('div','RdBu',20)));
    beautify([18 18 20]);

    for tt = 2:3:ntimes
        set(hh,'ZData', double(vorsurf(ix,iy,tt)));
        %set(hcsd, 'ZData', double(runs.csdsurf(ix,iy,tt)));
        %set(h0vor, 'ZData', double(runs.vorsurf(:,:,tt)));
        shading flat;
        runs.update_eddy_contour(hedd, tt);
        % runs.update_eddy_sshcontour(hssh, tt);
        set(ht,'String',['Surface vorticity | t = ' ...
                         num2str(runs.time(tt)/86400) ' ' ...
                         'days']);
        runs.update_timelabel(tt);
        runs.video_update();
        pause(0.05);
    end

    runs.video_write();
end
