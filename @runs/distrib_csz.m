% distribution of cs-z dyes
function [] = distrib_csz(runs)

% upper y-limit to save memory
    yend = find_approx(runs.rgrid.y_rho(:,1),130*1000);
    t0 = 65;runs.eddy.trevind;
    read_start = [1 1 1 t0-20];
    read_count = [Inf yend Inf 30];
    tindices = [t0 t0+read_count(end)-1];

    % read to calculate depth integrated upwelling/downwelling
    % before time loop
    w = dc_roms_read_data(runs.dir, 'w', tindices, {'y' 1 yend},[],runs.rgrid);

    % co-ordinate axes

    %[grd.xax,grd.yax,grd.zax,~,~,~] = dc_roms_var_grid(runs.rgrid,'temp');
    %grd.xax = grd.xax(:,1:yend,:);
    %grd.yax = grd.yax(:,1:yend,:);
    %grd.zax = grd.zax(:,1:yend,:);

    % grid matrices required for plotting
    xr = runs.rgrid.xr(:,1:yend)/1000; yr = runs.rgrid.yr(:,1:yend)/1000;
    ix = repmat([1:size(xr,1)]',[1 yend]);
    iy = repmat([1:yend],[size(xr,1) 1]);
    yzw = repmat(yr(1,:)', [1 runs.rgrid.N+1]);
    yzr = repmat(yr(1,:)', [1 runs.rgrid.N]);
    zw = permute(runs.rgrid.z_w(:,1:yend,:),[3 2 1]);

    % NEED TO ACCOUNT FOR TILTING IN VERTICAL?
    cx = runs.eddy.cx(t0:t0+read_count(end)-1)/1000;
    cy = runs.eddy.cy(t0:t0+read_count(end)-1)/1000;
    ee = runs.eddy.ee(t0:t0+read_count(end)-1)/1000;
    % hack if eddy center is outside extracted domain
    cy(cy > max(yr(:))) = max(yr(:));
    cxind = vecfind(xr(:,1),cx);
    cyind = vecfind(yr(1,:),cy)';

    % vertically integrated w - plan view - in streamer
    WS = squeeze( nansum( bsxfun(@times, ...
                                 bsxfun(@times,avg1(w,3), permute(streamer2,[1 2 4 3])), ...
                                 diff(zw,1,3) ), 3) );

    hfig = figure;
    maximize();

    for tt = 1:size(streamer2,3)
        % streamer has been identified - now extract data section
        volume = {'x' min(ixstr) max(ixstr);
                  'y' min(iystr) max(iystr)};

        %wstr = avg1(dc_roms_read_data(runs.dir, 'w', t0+tt-1,volume),3);
        % w was read earlier - just extract once
        wstr = w(volume{1,2}:volume{1,3}, volume{2,2}:volume{2,3}, :,tt);
        zdye = dc_roms_read_data(runs.dir, runs.zdname, t0+tt-1,volume,[],runs.rgrid);
        zr = permute(runs.rgrid.z_r(:,volume{2,2}:volume{2,3}, ...
                                    volume{1,2}:volume{1,3}),[3 2 1]);

        sz = [size(wstr,1) size(wstr,2)];
        wstr = reshape(wstr, sz(1) * sz(2), size(wstr,3));
        zdye = reshape(zdye, sz(1) * sz(2), size(zdye,3));
        zr = reshape(zr, sz(1) * sz(2), size(zr,3));

        % extract streamer section - indicated by suffix 'ex'
        inc = sub2ind(sz, ixstr - min(ixstr(:)) + 1, ...
                      iystr - min(iystr(:)) + 1);
        wex = wstr(inc,:);
        zrex = zr(inc,:);
        zdyeex = zdye(inc,:) - zrex;
        xex = repmat(dstr,[1 size(zrex,2)]);

        % index of western & eastern edges
        %wind = vecfind(xr(:,1), runs.eddy.vor.we/1000);
        %eind = vecfind(xr(:,1), runs.eddy.vor.ee/1000);

        % colorbar for vertical vel cross-section
        %wcolor = sort( [-1 1  ] * max(max(abs( ...
        %                    log10(abs(w(sort([eind wind]),:))) ))) )/2;

        %% animate depth integrated w in streamer

        %windex = wind(tindex)-dx; % for cross-section
        %eindex = eind(tindex)-dx; % for cross-section
        tindex = t0+tt-1;
        zlimit = [-1000 0];

        figure(hfig);
        if tt == 1
            subplot(221)
            titlestr = 'Depth integrated w in streamer (blue)';
            hws = pcolorcen(xr,yr,double(WS(:,:,ii))); shading flat;
            hold on;
            [~,hs] = contour(xr,yr,repnan(streamer(:,:,40,ii),0), ...
                             1,'b','LineWidth',2);
            he = runs.plot_eddy_contour('contour',tindex);
            hstr = plot(xstr,ystr,'kx');
            runs.plot_bathy('contour','k');
            colormap(flipud(cbrewer('div','RdBu',32)));
            caxis([-1 1] * max(abs([nanmin(WS(:)) nanmax(WS(:))])));
            colorbar; %cbunits('m^2/s');
            ht = runs.set_title(titlestr,tindex);

            % depth of 'streamer'
            subplot(223)
            hz = pcolorcen(xr,yr,double(max(abs(zs(:,:,:,ii)),[],3)));
            hold on;
            hcb = colorbar;  caxis([0 max(abs(zs(:)))]);cbunits('[m]');
            hzeta = runs.plot_zeta('contour',tindex);
            title('Depth of ''streamer''');

            % zdye - streamer section
            subplot(222)
            [~,hzdye] = contourf(xex,zrex,zdyeex);
            colorbar;
            ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
            title('\Delta z-dye');

            % vertical vel - streamer section
            subplot(224)
            [~,hw] = contourf(xex,zrex,avg1(wex,2));
            colorbar;
            ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
            title('vertical velocity');

            spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
            pause(0.01);

        else
            set(hws ,'CData',double(WS(:,:,tt)));
            set(hs  ,'ZData',repnan(streamer(:,:,40,tt),0));

            set(hz  ,'CData',double(max(abs(zs(:,:,:,tt)),[],3)));

            set(hstr,'XData',xstr,'YData',ystr);

            % streamer sections
            set(hzdye,'XData',xex,'YData',zrex,'ZData',zdyeex);
            set(hw  , 'XData',xex,'YData',zrex, 'ZData',avg1(wex,2));

            runs.update_zeta(hzeta,tindex);
            runs.update_eddy_contour(he, tindex);
            runs.update_title(ht,titlestr,tindex);
            pause(0.01);
        end
    end
end
