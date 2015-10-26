function [] = animate_center(runs)
    runs.video_init('center');
    eddy = runs.eddy;
    xvec = runs.rgrid.xr(:,1);
    yvec = runs.rgrid.yr(1,:)';

    varname = 'rho';

    % stride values
    % if y is cross-isobath, sx = st, sy = sxy & vice versa
    sxy = 10;
    sz = 1;
    st = 2;

    % this does not work yet.
    t0 = 1;

    ix = vecfind(xvec,eddy.mx([t0:st:end]));
    iy = vecfind(yvec,eddy.my([t0:st:end]));

    ixmax = max(ix); ixmin = min(ix);
    iymax = max(iy); iymin = min(iy);

    if runs.bathy.axis == 'x'
        stride = [sxy 1 sz st];
        temper = dc_roms_read_data(runs.dir,varname,[t0 st Inf], ...
                                   {'y' iymin iymax},stride,runs.rgrid, ...
                                   'his', 'single');
        strat = dc_roms_read_data(runs.dir,varname,[1 1], ...
                                  {'y' Inf Inf},stride, runs.rgrid, 'his', ...
                                  'single');

        temper = bsxfun(@minus,temper,permute(strat,[1 3 2]));
        %temper = roms_read_data(runs.dir,varname,[1 iymin 1 t0], ...
        %                  ceil([Inf iymax-iymin+1 Inf Inf]./stride), stride);
        %              toc;
        %strat  = roms_read_data(runs.dir,varname,[Inf 1 1 1], ...
        %                  ceil([1 1 Inf 1]./stride),stride);
        %              toc

    else
        stride = [1 sxy sz st];
        temper = dc_roms_read_data(runs.dir,varname,[t0 st Inf], ...
                                   {'x' ixmin ixmax},stride, runs.rgrid, ...
                                   'his', 'single');
        strat = dc_roms_read_data(runs.dir,varname,[1 1], ...
                                  {'y' Inf Inf},stride, runs.rgrid, 'his', ...
                                  'single');
        temper = bsxfun(@minus,temper,permute(strat,[3 1 2]));
        %temper = roms_read_data(runs.dir,varname,[ixmin 1  1 t0], ...
        %                ceil([ixmax-ixmin+1 Inf Inf Inf]./stride),stride);
        %            toc;
        %strat  = roms_read_data(runs.dir,varname,[1 1 1 1], ...
        %                ceil([1 Inf Inf 1]./stride),stride);
        %            toc;
    end


    % make plot
    tt = 1;
    figure;
    % first plan view of zeta
    subplot(211)
    hz = runs.plot_zeta('pcolor',tt);
    shading interp
    hold on
    colorbar; freezeColors;
    hb = runs.plot_bathy('contour','k');
    he = runs.plot_eddy_contour('contour',tt);
    ht1 = title(['Free surface | ' num2str(runs.rgrid.ocean_time(tt)/86400)  ' days']);
    xlabel('X (km)');ylabel('Y (km)');
    axis image;
    beautify([16 16 18]);

    % temp following eddy center
    levels = linspace(min(temper(:)),max(temper(:)),25);
    subplot(212)
    if runs.bathy.axis == 'x'
        xzr = repmat(xvec(1:stride(1):end,1),[1 size(temper,3)]);
        [~,hh] = contourf(xzr/1000,squeeze(runs.rgrid.zr(1:stride(1):end,iy(1),:)), ...
                          squeeze(temper(:,iy(1)-iymin + 1,:,1)),levels);
    else
        yzr = repmat(yvec(1:stride(2):end),[1 size(temper,3)]);
        [~,hh] = contourf(yzr/1000,squeeze(runs.rgrid.zr(ix(1),1:stride(2):end,:)), ...
                          squeeze(temper(ix(1)-ixmin + 1,:,:,1)),levels);
    end
    %ht = title(['(mx,my) = (', num2str(eddy.mx(stride(4))/1000) ',' ...
    %        num2str(eddy.my(tt*stride(4))/1000) ') km | t = ' num2str(stride(4)) ' days']);
    xlabel('y (km)'); ylabel('z (m)'); colorbar;
    %caxis([-1 1]*max(mat2vec(abs(temper(ix-ixmin+1,:,:,1:end-10)))));
    caxis([-1 1] *max(abs(temper(:))));
    h1 = liney(-eddy.Lz2(stride(4)),[],'b');
    ylim([-1500 0]);
    title('Cross-shore temperature anomaly - slice through eddy center');
    %h2 = liney(-eddy.Lz3(stride(4)),'3','k');
    maximize(gcf); pause(0.2);
    beautify([16 16 18]);
    runs.video_update();
    % update plots
    for tt=2:size(temper,4)
        if runs.bathy.axis == 'y'
            set(hh,'YData',squeeze(runs.rgrid.zr(ix(tt),1:stride(2):end,:)));
            set(hh,'ZData',squeeze(temper(ix(tt)-ixmin + 1,:,:,tt)));
        else
            set(hh,'ZData',squeeze(temper(:,iy(tt)-iymin + 1,:,tt)));
        end
        tstr = [num2str(runs.time(tt*stride(4))/86400) ' days'];
        set(h1,'ydata',[-eddy.Lz2(tt*stride(4)) -eddy.Lz2(tt*stride(4))]);
        runs.update_zeta(hz,tt*stride(4));

        runs.update_eddy_contour(he,tt*stride(4));
        %set(ht,'String', ['(mx,my) = (', num2str(eddy.mx(tt*stride(4))/1000) ',' ...
        %    num2str(eddy.my(tt*stride(4))/1000) ') | t = ' tstr]);
        set(ht1,'String',['Free surface | ' tstr]);
        runs.video_update();
        pause(0.01);
    end

    runs.video_write();
end
