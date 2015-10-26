% depth section through streamer
function [] = animate_streamer_section(runs)

    debug_plot = 0;
    try
        if ~isfield(runs.streamer.west,'mask')
            runs.build_streamer_section();
        end
    catch
        runs.build_streamer_section();
    end
    yend = runs.streamer.yend;
    t0 = 65;runs.eddy.trevind;
    tend = t0+30;
    %read_count = [Inf yend Inf 30];
    tindices = [t0 tend];

    sz4dfull = runs.streamer.sz4dfull;
    sz4dsp = runs.streamer.sz4dsp;
    sz3dfull = runs.streamer.sz3dfull;
    sz3dsp = runs.streamer.sz3dsp;

    sz4dfull(4) = tend-t0+1;
    sz4dsp(2) = tend-t0+1;
    sz4d3d= [sz4dfull(1)*sz4dfull(2) sz4dfull(3) sz4dfull(4)];
    sz3d2d = sz4d3d(1:2);

    % read to calculate depth integrated upwelling/downwelling
    % before time loop
    w = avg1(dc_roms_read_data(runs.dir, 'w', tindices, ...
                               {'y' 1 yend},[],runs.rgrid),3);
    wstr = reshape(w,sz4dsp) .* runs.streamer.west.mask(:,t0:tend);
    clear w

    % grid matrices required for plotting
    xr = runs.rgrid.xr(:,1:yend)/1000; yr = runs.rgrid.yr(:,1:yend)/1000;
    zw = permute(runs.rgrid.z_w(:,1:yend,:),[3 2 1]);
    zr = permute(runs.rgrid.z_r(:,1:yend,:),[3 2 1]);

    % vertically integrated w - plan view - in streamer
    WS = squeeze( nansum( bsxfun(@times, ...
                                 reshape(full(wstr),sz4dfull), ...
                                 diff(zw,1,3) ), 3) );

    hfig = figure;
    maximize();

    for tt = 1:sz4dsp(end)
        tind = t0+tt-1;

        % get section locations & make grid matrices
        xstr = runs.streamer.west.xstr{tind};
        ystr = runs.streamer.west.ystr{tind};
        dstr = repmat(runs.streamer.west.dstr{tind},[1 runs.rgrid.N]);

        ixmin = find_approx(xr(:,1),min(xstr));
        ixmax = find_approx(xr(:,1),max(xstr));
        iymin = find_approx(yr(1,:),min(ystr));
        iymax = find_approx(yr(1,:),max(ystr));

        ixmin = max(ixmin-5,1);
        ixmax = min(ixmax+5,size(runs.bathy.h,1));
        iymin = max(iymin-5,1);
        iymax = min(iymax+5,size(runs.bathy.h,2));

        % streamer has been identified - now extract data section
        volume = {'x' ixmin ixmax;
                  'y' iymin iymax};

        tindex = t0+tt-1;
        zlimit = [-1000 0];

        streamer = reshape(full(runs.streamer.west.mask(:,t0+tt-1)) ...
                           ,sz3dfull);

        % read velocities & dyes in block form
        sznew3d = [(ixmax-ixmin+1) (iymax-iymin+1) 40];
        sznew2d = [sznew3d(1)*sznew3d(2) sznew3d(3)];
        [u,xumat,yumat,zumat] = dc_roms_read_data(runs.dir,'u', ...
                                                  tind,volume,[],runs.rgrid);
        [v,xvmat,yvmat,zvmat] = dc_roms_read_data(runs.dir,'v', ...
                                                  tind,volume,[],runs.rgrid);
        [csdye,xrmat,yrmat,zrmat] = dc_roms_read_data(runs.dir, runs.csdname, ...
                                                      tind,volume,[],runs.rgrid);
        zdye = dc_roms_read_data(runs.dir, runs.zdname, ...
                                 tind,volume,[],runs.rgrid);

        xumat = xumat/1000; yumat = yumat/1000;
        xvmat = xvmat/1000; yvmat = yvmat/1000;
        xrmat = xrmat/1000; yrmat = yrmat/1000;

        if runs.streamer.west.fit_circle
            N = runs.rgrid.N;

            % bathymetry along streamer
            bstr = interp2(xr',yr',runs.bathy.h(:,1:yend)',xstr,ystr);

            %               % zgrid along streamer - RHO points
            zstr = squeeze(set_depth(2,4,runs.rgrid.theta_s, ...
                                       runs.rgrid.theta_b,runs.rgrid.hc,N,1,bstr,...
                                       zeros(size(bstr)),0));

            [I.XR,I.YR] = ndgrid(xstr,ystr);
            hin = interp2(xr',yr',runs.bathy.h(:,1:yend)',I.XR,I.YR);
            zetain = interp2(xr',yr',runs.zeta(:,1:yend,tind)',I.XR,I.YR);
            I.ZR = set_depth(2,4,runs.rgrid.theta_s, ...
                               runs.rgrid.theta_b,runs.rgrid.hc,N,1,hin,...
                               zetain,0);

            % structure for interp_field.m
            % doesn't change
            I.Vname = 'does not matter';
            I.nvdims = ndims(u);
            I.Dmask = ones(size(u)); I.Rmask = ones(size(u));
            I.Zsur = max(I.ZR(:));
            I.Zbot = min(I.ZR(:));
            % indices to extract section
            lstr = length(xstr);
            indin = sub2ind([lstr lstr],[1:lstr],[1:lstr]);

            % now interp variables
            I.VD = u;
            I.XD = xumat; I.YD = yumat; I.ZD = zumat;
            ustr = reshape(interp_field(I),[lstr*lstr N]);
            ustr = ustr(indin,:);

            I.VD = v;
            I.XD = xvmat; I.YD = yvmat; I.ZD = zvmat;
            vstr = reshape(interp_field(I),[lstr*lstr N]);
            vstr = vstr(indin,:);

            I.VD = zdye;
            I.XD = xrmat; I.YD = yrmat; I.ZD = zrmat;
            zdstr = reshape(interp_field(I),[lstr*lstr N]);
            zdstr = zdstr(indin,:);

            I.VD = csdye;
            csstr = reshape(interp_field(I),[lstr*lstr N]);
            csstr = csstr(indin,:);

            I.VD = streamer(ixmin:ixmax,iymin:iymax,:);
            strstr = reshape(interp_field(I),[lstr*lstr N]);
            strstr = round(strstr(indin,:));

            % first interpolate in horizontal on original grid levels
            %                 xin = nan([numel(zstr) 1]);
            %                 yin = xin; zin = xin;
            %                 % build grid vectors
            %                 for mmm=1:length(xstr)
            %                     start = N*(mmm-1) + 1;
            %                     stop = start+N-1;
            %
            %                     xin(start:stop) = xstr(mmm);
            %                     yin(start:stop) = ystr(mmm);
            %                     zin(start:stop) = zstr(mmm,:);
            %                 end
            %                 for nn=1:N
            %                     nel = [numel(xumat(:,:,nn)) 1];
            %                     Fu = scatteredInterpolant( ...
            %                         reshape(xumat(:,:,nn), nel), ...
            %                         reshape(yumat(:,:,nn), nel), ...
            %                         reshape(u(:,:,nn), nel));
            %                     ui(:,nn) = Fu(xstr,ystr);
            %                 end

            % now interpolate
            %                 Fu = scatteredInterpolant(xumat(:),yumat(:),zumat(:),u(:));
            %                 ustr = reshape(Fu(xin,yin,zin), [N numel(xin)/N])';
            %
            %                 Fv = scatteredInterpolant(xvmat(:),yvmat(:),zvmat(:),v(:),'nearest');
            %                 vstr = reshape(Fv(xin,yin,zin), [N numel(xin)/N])';
            %
            %                 Fcs = scatteredInterpolant(xrmat(:),yrmat(:),zrmat(:),csdye(:),'nearest');
            %                 csstr = reshape(Fcs(xin,yin,zin), [N numel(xin)/N])';
            %
            %                 Fz = scatteredInterpolant(xrmat(:),yrmat(:),zrmat(:),zdye(:),'nearest');
            %                 zdstr = reshape(Fz(xin,yin,zin), [N numel(xin)/N])';

            %                 % interpolating streamer mask doesn't work
            %                 zmin = min(runs.streamer.zr .* streamer,[],3);
            %                 zminstr = interp2(xr',yr',zmin',xstr,ystr);
            %                 strstr = bsxfun(@gt,zstr,zminstr);
            %strex = streamer(ixmin:ixmax,iymin:iymax,:);
            %F = scatteredInterpolant(xrmat(:),yrmat(:),zrmat(:),strex(:),'linear');
            %strstr = round(reshape(F(xin,yin,zin), [N numel(xin)/N])');

        else
            ixstr = runs.streamer.west.ixstr{tind};
            iystr = runs.streamer.west.iystr{tind};
            indices = sub2ind(sz4dfull(1:2),ixstr,iystr);

            zlin = reshape(zr,sz3d2d);
            zstr = zlin(indices,:);
            clear zlin;

            bstr = runs.bathy.h(indices)';

            % streamer mask vertical section - along-streamer section
            % points
            strlin = reshape(streamer,sz3d2d);
            strstr = strlin(indices,:);

            ixnew = ixstr - min(ixstr(:)) + 1;
            iynew = iystr - min(iystr(:)) + 1;
            % extract variables at streamer points
            u = reshape(u,sznew2d);
            v = reshape(v,sznew2d);
            csdye = reshape(csdye,sznew2d);
            zdye = reshape(csdye,sznew2d);
            indnew = sub2ind(sznew3d(1:2),ixnew,iynew);
            ustr = u(indnew,:);
            vstr = v(indnew,:);
            zdstr = zdye(indnew,:);
            csstr = csdye(indnew,:);
        end

        % bathy-patch
        bpatch = [-bstr' -max(runs.bathy.h(:))-100 ...
                  -max(runs.bathy.h(:))-100];
        dpatch = [dstr(:,1)' dstr(end,1) 0];

        % streamer mask at surface
        streamer = streamer(:,:,40);

        % rotate velocities to along & cross-streamer dirns.
        angle = atan2d(diff(ystr),diff(xstr));
        angle(end+1) = angle(end);
        angle = repmat(angle,[1 size(ustr,2)]);
        if debug_plot
            figure;
            plot(xstr,ystr); hold on;
            dx = 4;
            for ii=1:size(xstr,1)
                text(xstr(ii),ystr(ii),num2str(angle(ii,1)));
            end
        end
        % normal vel
        Unstr = ustr .* cosd(angle) - vstr .* sind(angle);
        % tangential vel
        Utstr = ustr .* sind(angle) + vstr .* cosd(angle);

        % replace values in the vertical that aren't associated with
        % the streamer with NaNs
        %Utstr(strstr == 0) = NaN;
        %Unstr(strstr == 0) = NaN;
        %zdstr(strstr == 0) = NaN;
        %csstr(strstr == 0) = NaN;

        figure(hfig);
        if tt == 1
            limy = [0 nanmax(cat(1,runs.streamer.west.ystr{:}))+ ...
                    10*runs.rgrid.dy/1000];
            limx = [nanmin(cat(1,runs.streamer.west.xstr{:})) ...
                    400]; % CHANGE THIS
            limz = [-1000 0];

            % normalized depth integrated w
            ax(1) = subplot(231);
            titlestr = 'NORMALIZED \int w dz in streamer (blue)';
            hws = pcolorcen(xr,yr,double(WS(:,:,tt))./...
                            nanmax(nanmax(abs(WS(:,:,tt))))); shading flat;
            hold on;
            [~,hs] = contour(xr,yr,repnan(streamer,0), ...
                             1,'b','LineWidth',2);
            he = runs.plot_eddy_contour('contour',tindex);
            hstr = plot(xstr,ystr,'kx');
            runs.plot_bathy('contour','k');
            colormap(flipud(cbrewer('div','RdBu',32)));
            caxis([-1 1]);
            xlim(limx); ylim(limy);
            %caxis([-1 1] * max(abs([nanmin(WS(:)) nanmax(WS(:))])));
            colorbar; %cbunits('m^2/s');
            ht = runs.set_title(titlestr,tindex);

            % un=normalized depth integrated w
            ax(2) = subplot(234);
            hws2 = pcolorcen(xr,yr,double(WS(:,:,tt))); shading flat;
            hold on;
            [~,hs2] = contour(xr,yr,repnan(streamer,0), ...
                              1,'b','LineWidth',2);
            he2 = runs.plot_eddy_contour('contour',tindex);
            runs.plot_bathy('contour','k');
            colormap(flipud(cbrewer('div','RdBu',32)));
            title('\int w dz in streamer (blue)');
            xlim(limx); ylim(limy);
            caxis([-1 1] * max(abs([nanmin(WS(:)) nanmax(WS(:))])));
            colorbar; %cbunits('m^2/s');

            % zdye - streamer section
            ax(3) = subplot(232);
            [~,hzdye] = contourf(dstr,zstr,zdstr - zstr);
            colorbar; ylim(limz); caxis([-50 50]);
            ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
            title('\Delta z-dye');
            hold on;
            [~,hstrz1] = contour(dstr,zstr,strstr,[1 1],'k');
            set(hstrz1,'LineWidth',2);
            hpatch(3) = patch(dpatch,bpatch,'k');

            % cross-shelf dye - streamer section
            ax(4) = subplot(235);
            [~,hcsd] = contourf(dstr,zstr,csstr/1000 - runs.bathy.xsb/1000);
            colorbar; ylim(limz); caxis([-10 40]);
            ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
            title('Cross-shelf dye - X_{shelfbreak} (km)');
            hold on;
            [~,hstrz2] = contour(dstr,zstr,strstr,[1 1],'k');
            set(hstrz2,'LineWidth',2);
            hpatch(4) = patch(dpatch,bpatch,'k');

            % velocities - streamer section
            ax(5) = subplot(233);
            [~,hun] = contourf(dstr,zstr,Unstr);
            colorbar; ylim(limz);caxis([-1 1]*0.1);
            ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
            title('Normal velocity (m/s)');
            hold on;
            [~,hstrz3] = contour(dstr,zstr,strstr,[1 1],'k');
            set(hstrz3,'LineWidth',2);
            hpatch(5) = patch(dpatch,bpatch,'k');

            ax(6) = subplot(236);
            [~,hut] = contourf(dstr,zstr,Utstr);
            colorbar; ylim(limz); caxis([-1 1]*0.1);
            ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
            title('Tangential velocity (m/s)');
            hold on;
            [~,hstrz4] = contour(dstr,zstr,strstr,[1 1],'k');
            set(hstrz4,'LineWidth',2);
            hpatch(6) = patch(dpatch,bpatch,'k');

            %spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
            pause();

        else
            set(hws ,'CData',double(WS(:,:,tt))./...
                     nanmax(nanmax(abs(WS(:,:,tt)))));
            set(hs  ,'ZData',repnan(streamer,0));
            set(hstr,'XData',xstr,'YData',ystr);

            for mmm=3:6
                set(ax(mmm),'XLim',[0 max(dstr(:,1))]);
                set(hpatch(mmm),'XData',dpatch,'YData',bpatch);
            end

            set(hws2, 'CData', double(WS(:,:,tt)));
            set(hs2  ,'ZData',repnan(streamer,0));

            set(hzdye,'XData',dstr,'YData',zstr,'ZData', zdstr - zstr);
            set(hcsd ,'XData',dstr,'YData',zstr,'ZData', ...
                      csstr/1000 - runs.bathy.xsb/1000);

            set(hun ,'XData',dstr,'YData',zstr,'ZData',Unstr);
            set(hut ,'XData',dstr,'YData',zstr,'ZData',Utstr);

            % update streamer depth contour
            set(hstrz1,'XData',dstr,'YData',zstr,'ZData',strstr);
            set(hstrz2,'XData',dstr,'YData',zstr,'ZData',strstr);
            set(hstrz3,'XData',dstr,'YData',zstr,'ZData',strstr);
            set(hstrz4,'XData',dstr,'YData',zstr,'ZData',strstr);

            %runs.update_zeta(hzeta,tindex);
            runs.update_eddy_contour(he2, tindex);
            runs.update_eddy_contour(he, tindex);
            runs.update_title(ht,titlestr,tindex);
            pause();
        end
    end
end
