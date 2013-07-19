classdef runs < handle
    properties
        % dir & file names
        dir; out_file; ltrans_file; flt_file; givenFile
        % data
        zeta; temp; dye; usurf; vsurf
        % grid & bathymetry
        rgrid; bathy
        % float data
        roms; ltrans;
        % eddy track data
        eddy;
        % initial params
        params
        % transport
        eutrans;
        % make video?
        makeVideo; mm_instance;
        %
        comment = ['prox = distance of edge from shelfbreak in m | '...
                   'eutrans = eulerian transport estimate (structure).'];
    end
    methods
        % constructor
        function [runs] = runs(dir)
            if isdir(dir)
                runs.dir = dir;
                runs.out_file = [dir '/ocean_avg.nc'];
            else
                runs.givenFile = 1;
                runs.out_file = dir;
                dir = strrep(dir,'\','/');
                inds = strfind(dir,'/');
                dir = dir(1:inds(end));
                runs.dir = dir;
            end
            runs.flt_file = [dir '/ocean_flt.nc'];
            runs.ltrans_file = [dir '/ltrans.nc'];
            runs.rgrid = roms_get_grid(runs.out_file,runs.out_file,0,1);
            runs.rgrid.xr = runs.rgrid.x_rho';
            runs.rgrid.yr = runs.rgrid.y_rho';
            runs.rgrid.zr = permute(runs.rgrid.z_r,[3 2 1]);
            runs.makeVideo = 0; % no videos by default.
            
            if ~runs.givenFile
                runs.zeta = roms_read_data(dir,'zeta');
                try
                    runs.dye  = roms_read_data(dir,'dye_02');
                catch ME
                end
            else
                runs.zeta = double(ncread(runs.out_file,'zeta'));
                try
                    runs.dye  = squeeze(double(ncread(runs.out_file,'dye_02', ...
                    [1 1 runs.rgrid.N 1],[Inf Inf 1 Inf])));
                catch ME
                end
            end
            % params & bathy
            runs.params = read_params_from_ini(runs.dir);
            runs.bathy = runs.params.bathy;
            
            % fill bathy
            [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = find_shelfbreak(runs.out_file);
            [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = find_shelfbreak(runs.out_file,'slope');
            runs.bathy.h = runs.rgrid.h';
            
            if exist(runs.flt_file,'file')
                runs.roms = floats('roms',runs.flt_file,runs.rgrid);
            end
            
            if ~exist([dir '/eddytrack.mat'],'file')
                try
                    runs.eddy = track_eddy(dir);
                catch ME
                    disp('Couldn''t run track_eddy.m');
                end
            else
                if strfind(runs.out_file,'_004.nc')
                    edd = load([dir '/eddytrack_004.mat'],'eddy');
                else
                    edd = load([dir '/eddytrack.mat'],'eddy');
                end
                runs.eddy = edd.eddy;
                if runs.bathy.axis == 'y'
                    edge = runs.eddy.se;
                else
                    edge = runs.eddy.we;
                end
                runs.eddy.prox = (edge-runs.bathy.xsb);
                h = runs.bathy.h(2:end-1,2:end-1);
            
                ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx);
                iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.my);

                runs.eddy.hcen = h(sub2ind(size(runs.eddy.xr),ix,iy))';
            end
            
            if exist(runs.ltrans_file,'file')
                runs.ltrans = floats('ltrans',runs.ltrans_file,runs.rgrid);
            end
        end
        
        function [] = animate_zeta(runs)
            if runs.makeVideo
%                 runs.mm_instance = mm_setup;
%                 runs.mm_instance.pixelSize = [1600 900];
%                 runs.mm_instance.outputFile = 'mm_output.avi';
%                 runs.mm_instance.ffmpegArgs = '-q:v 1 -g 1';
%                 runs.mm_instance.InputFrameRate = 3;
%                 runs.mm_instance.frameRate = 3;
                aviobj = VideoWriter('output','MPEG-4');
                open(aviobj);
            end
            figure;
            hz = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,1));
            shading flat
            ax = gca;
            hold on
            colorbar; freezeColors;
            [cc,hh] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.rgrid.h',[200 500 1000 1500 2000],'k');
            clabel(cc,hh);
            [~,he] = contour(runs.eddy.xr/1000,runs.eddy.yr/1000,runs.eddy.mask(:,:,1));
            if runs.bathy.axis == 'y'
                liney(runs.bathy.xsb/1000,'shelfbreak','w');
            else
                linex(runs.bathy.xsb/1000,'shelfbreak','w');
            end
            ht = title([num2str(runs.rgrid.ocean_time(1)/86400)  ' days']);
            xlabel('X (km)');ylabel('Y (km)');
            if runs.makeVideo
                %shading(gca,'interp');
                disp('maximize!');
                pause; 
                mm_addFrame(runs.mm_instance,gcf);
            end
            for ii = 2:size(runs.zeta,3)
                set(hz,'CData',runs.zeta(:,:,ii));
                if runs.makeVideo, shading interp; end
                set(he,'ZData',runs.eddy.mask(:,:,ii));
                set(ht,'String',[num2str(runs.rgrid.ocean_time(ii)/86400) ' days']);
                if runs.makeVideo
                   % shading(gca,'interp');
                    %mm_addFrame(runs.mm_instance,gcf);
                    F = getframe(gcf);
                    writeVideo(aviobj,F);
                end
                pause(0.03);
            end
            if runs.makeVideo
               % mm_render(runs.mm_instance);
               close(aviobj);
            end
        end
        
        function [] = animate_center(runs)
           if runs.makeVideo
%                 runs.mm_instance = mm_setup;
%                 runs.mm_instance.pixelSize = [1600 900];
%                 runs.mm_instance.outputFile = 'mm_output.avi';
%                 runs.mm_instance.ffmpegArgs = '-q:v 1 -g 1';
%                 runs.mm_instance.InputFrameRate = 3;
%                 runs.mm_instance.frameRate = 3;
                aviobj = VideoWriter('animate_center','MPEG-4');
                set(aviobj,'FrameRate',8,'Quality',100);
                open(aviobj);
            end
            eddy = runs.eddy;
            xvec = runs.rgrid.xr(:,1);
            yvec = runs.rgrid.yr(1,:)';
            stride = [1 1 1 1];
            
            ix = vecfind(xvec,eddy.mx([1:stride(4):end]));
            iy = vecfind(yvec,eddy.my([1:stride(4):end]));
            
            ixmax = max(ix); ixmin = min(ix);
            iymax = max(iy); iymin = min(iy);            

            if runs.bathy.axis == 'x'
                temper = double(squeeze(ncread(runs.out_file,'temp',[1 iymin 1 1], ...
                                [Inf iymax-iymin+1 Inf Inf],stride)));
            else
                temper = double(squeeze(ncread(runs.out_file,'temp',[ixmin 1  1 1], ...
                                [ixmax-ixmin+1 Inf Inf Inf],stride)));
            end

            
            tt = 1;
            figure;
            % first plan view of zeta
            subplot(211)
            hz = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,1));
            shading interp
            ax = gca;
            hold on
            colorbar; freezeColors;
            [cc,hb] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.rgrid.h',[200 500 1000 1500 2000],'k');
            clabel(cc,hb);
            [~,he] = contour(runs.eddy.xr/1000,runs.eddy.yr/1000,runs.eddy.mask(:,:,1));
            if runs.bathy.axis == 'y'
                liney(runs.bathy.xsb/1000,'shelfbreak','k');
            else
                linex(runs.bathy.xsb/1000,'shelfbreak','k');
            end
            ht1 = title(['Free surface | ' num2str(runs.rgrid.ocean_time(1)/86400)  ' days']);
            xlabel('X (km)');ylabel('Y (km)');
            axis image;
            beautify([14 14 16]);
            
            % temp following eddy center
            subplot(212)
            if runs.bathy.axis == 'x'
                xzr = repmat(xvec(1:stride(1):end,1),[1 size(temper,3)]);
                [~,hh] = contourf(xzr/1000,squeeze(runs.rgrid.zr(1:stride(1):end,iy(1),:)), ...
                         squeeze(temper(:,iy(1)-iymin + 1,:,1)),10);
            else
                yzr = repmat(yvec(1:stride(2):end),[1 size(temper,3)]);
                [~,hh] = contourf(yzr/1000,squeeze(runs.rgrid.zr(ix(1),1:stride(2):end,:)), ...
                         squeeze(temper(ix(1)-ixmin + 1,:,:,1)),20);
            end
            %ht = title(['(mx,my) = (', num2str(eddy.mx(stride(4))/1000) ',' ...
            %        num2str(eddy.my(tt*stride(4))/1000) ') km | t = ' num2str(stride(4)) ' days']);
            xlabel('y (km)'); ylabel('z (m)'); colorbar; caxis([14 20]);
            h1 = liney(-eddy.Lz2(stride(4)),[],'b');
            ylim([-1000 0]);
            title('Cross-shore temperature slice through eddy center');
            %h2 = liney(-eddy.Lz3(stride(4)),'3','k');
            beautify([14 14  16]);
            if runs.makeVideo
                %shading(gca,'interp');
                disp('maximize!');
                pause; 
%               mm_addFrame(runs.mm_instance,gcf);
                F = getframe(gcf);
                writeVideo(aviobj,F);
            end
            % update plots
            for tt=2:size(temper,4)
                if runs.bathy.axis == 'y'
                    set(hh,'YData',squeeze(runs.rgrid.zr(ix(tt),1:stride(2):end,:)));
                    set(hh,'ZData',squeeze(temper(ix(tt)-ixmin + 1,:,:,tt)));
                else
                    set(hh,'ZData',squeeze(temper(:,iy(tt)-iymin + 1,:,tt)));
                end
                tstr = [num2str(runs.rgrid.ocean_time(tt*stride(4))/86400) ' days'];
                set(h1,'ydata',[-eddy.Lz2(tt*stride(4)) -eddy.Lz2(tt*stride(4))]);
                %set(h2,'ydata',-eddy.Lz3(tt*stride(4)));
                set(hz,'CData',runs.zeta(:,:,tt*stride(4)));
                %if runs.makeVideo, shading interp; end
                set(he,'ZData',runs.eddy.mask(:,:,tt*stride(4)));
                %set(ht,'String', ['(mx,my) = (', num2str(eddy.mx(tt*stride(4))/1000) ',' ...
                %    num2str(eddy.my(tt*stride(4))/1000) ') | t = ' tstr]);
                set(ht1,'String',['Free surface | ' tstr]);
                if runs.makeVideo
                   % shading(gca,'interp');
                    %mm_addFrame(runs.mm_instance,gcf);
                    F = getframe(gcf);
                    writeVideo(aviobj,F);
                end
                pause(0.01); 
            end
            
            if runs.makeVideo
               % mm_render(runs.mm_instance);
               close(aviobj);
            end
        end
        
        function [] = animate_pt(runs)
%             if runs.makeVideo
% %                 runs.mm_instance = mm_setup;
% %                 runs.mm_instance.pixelSize = [1600 900];
% %                 runs.mm_instance.outputFile = 'mm_output.avi';
% %                 runs.mm_instance.ffmpegArgs = '-q:v 1 -g 1';
% %                 runs.mm_instance.InputFrameRate = 3;
% %                 runs.mm_instance.frameRate = 3;
%                 aviobj = VideoWriter('animate_pt','MPEG-4');
%                 open(aviobj);
%             end
            if isempty(runs.usurf)
                if runs.givenFile
                    runs.usurf = double(squeeze(ncread(runs.out_file, ....
                        'u',[1 1 runs.rgrid.N 1],[Inf Inf 1 Inf])));
                else
                    error('no code here yet.');
                end
                runs.usurf = avg1(runs.usurf(:,2:end-1,:),1);
            end
            
            if isempty(runs.vsurf)
                if runs.givenFile
                    runs.vsurf = double(squeeze(ncread(runs.out_file, ....
                        'v',[1 1 runs.rgrid.N 1],[Inf Inf 1 Inf])));
                else
                    error('no code here yet.');
                end
                runs.vsurf = avg1(runs.vsurf(2:end-1,:,:),2);
            end
            
            
            rr = sqrt(runs.params.phys.N2)*runs.bathy.hsb/runs.rgrid.f(runs.bathy.isb,1);
            distance = 5*rr; % 5 times rossby radius
            clim = [runs.bathy.xsb/1000 runs.bathy.xsb/1000+distance/1000];
            
            i = 1; clf;
            hf = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.dye(:,:,i)/1000);         
            caxis(clim);
            shading flat;
            colorbar; hold on
            [~,hc] = contour(runs.eddy.xr/1000,runs.eddy.yr/1000,runs.eddy.mask(:,:,i),1,'k');
            dxi = 5; dyi = 5;
            hq = quiver(runs.eddy.xr(1:dxi:end,1:dyi:end)/1000,runs.eddy.yr(1:dxi:end,1:dyi:end)/1000, ...
                        runs.usurf(1:dxi:end,1:dyi:end,i),runs.vsurf(1:dxi:end,1:dyi:end,i));
            set(hc,'LineWidth',2);    
            if runs.bathy.axis == 'y'
                liney(runs.bathy.xsb/1000,'shelfbreak','w');
                liney(mean(runs.eddy.se)/1000,'mean(se)','k');
            else
                linex(runs.bathy.xsb/1000,'shelfbreak','w');
            end
            ht = title(num2str(i));
            xlabel('X (km)');ylabel('Y (km)');
            %axis image;
            beautify;
            pause();
            for i = 2:size(runs.zeta,3)
                set(hf,'CData',runs.dye(:,:,i)/1000);
                set(hc,'ZData',runs.eddy.mask(:,:,i));
                set(ht,'String',num2str(i));
                set(hq,'UData',runs.usurf(1:dxi:end,1:dyi:end,i));
                set(hq,'VData',runs.vsurf(1:dxi:end,1:dyi:end,i));
                pause();
            end
        end
        
        function [] = animate_floats(runs,type)
            if strcmpi(type,'ltrans')
                runs.ltrans.animate(runs.rgrid,runs.zeta,runs.eddy);
            end
            if strcmpi(type,'roms')
                runs.roms.animate(runs.rgrid,runs.zeta,runs.eddy);
            end
        end
        
        % create initial seed file for ltrans
        function [] = ltrans_create(runs)
            ltrans_create(runs.rgrid,runs.zeta,runs.eddy);
        end
        
        % create ltrans init file from roms out
        function [] = ltrans_create_from_roms(runs)
            ltrans_create_from_roms('ltrans_init_compare.txt',runs.flt_file,runs.rgrid);
        end
        
        function [] = compare_floats(runs)
            ltransc = floats('ltrans',[runs.dir '/ltrans-compare.nc'],runs.rgrid);
            runs.roms.plot_stats;
            ltransc.plot_stats;
        end
        
        function [] = transport(runs)
            % need some kind of initial time instant - probably objective
            % criterion based on distance between shelfbreak and eddy or
            % some sort
            runs.eutrans = [];
            t0 = 1;
            h = runs.bathy.h(2:end-1,2:end-1);
            
            ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx(t0:end));
            iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.my(t0:end));
            hcen = h(sub2ind(size(runs.eddy.xr),ix,iy))';
            
            ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx(t0:end));
            iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.se(t0:end));
            hedge = h(sub2ind(size(runs.eddy.xr),ix,iy))';
            % rossby radius
            rr = sqrt(runs.params.phys.N2)*runs.bathy.hsb/runs.rgrid.f(runs.bathy.isb,1);
            distance = 5*rr; % 5 times rossby radius
            
            if runs.params.bathy.axis == 'x'
                error(' not built for north-south isobaths');
            else
                loc = sort([nanmean(runs.eddy.se(t0:end)) nanmean(runs.eddy.cy(t0:end)) ...
                        runs.bathy.xsb  ... 
                        runs.rgrid.y_rho(vecfind(runs.bathy.h(1,:),[250 1000]),1)']);
            end
            
            runs.eutrans.x = loc;
            runs.eutrans.ix = vecfind(runs.rgrid.yr(1,:),loc);%find_approx(runs.rgrid.yr(1,:),loc,1);
            runs.eutrans.h = ceil(runs.bathy.h(1,runs.eutrans.ix));
            
            dx = runs.rgrid.xr(2,1)-runs.rgrid.xr(1,1);
            % cross-shore velocity
            for kk=1:length(loc)
                cs_vel = double(squeeze(ncread(runs.out_file,'v', ...
                    [1 runs.eutrans.ix(kk) 1 t0],[Inf 1 Inf Inf])));
                % dimensions = (x/y , z , t )
                mask = nan(size(cs_vel));
                iwest = vecfind(runs.eddy.xr(:,1),runs.eddy.we);
                for tt=1:size(cs_vel,3)
                    mask(1:iwest(tt),:,tt) = 1;
                end
                
                zmask = (abs(squeeze(runs.rgrid.z_r(:,runs.eutrans.ix(kk),:))   )' ...
                                < runs.bathy.hsb);
                mask = bsxfun(@times,mask,fillnan(zmask,0));

%                 runs.eutrans.trans(:,:,kk) = squeeze( ...
%                      trapz(runs.rgrid.z_r(:,runs.eutrans.ix(kk),1),mask .* cs_vel,2));
                runs.eutrans.trans(:,:,kk) = squeeze( ...
                            nansum(bsxfun(@times,avg1(mask .* cs_vel,2), ...
                            diff(runs.rgrid.z_r(:,runs.eutrans.ix(kk),1))'),2));

                runs.eutrans.Itrans(:,kk) = squeeze(nansum(runs.eutrans.trans(:,:,kk) .* dx,1))';
                
                % if I have passive tracer info I can calculate transport
                % using that
                if ~isempty(runs.dye)
                    mask = nan(size(cs_vel));
                    ieast = vecfind(runs.eddy.xr(:,1),runs.eddy.ee);
                    for tt=1:size(cs_vel,3)
                        mask(1:ieast(tt),:,tt) = 1;
                    end
                    
                    mask = bsxfun(@times,mask,fillnan(zmask,0));

                    dye = double(squeeze(ncread(runs.out_file,'dye_02', ...
                    [1 runs.eutrans.ix(kk) 1 t0],[Inf 1 Inf Inf])));
                    dyemask = (dye >= runs.bathy.xsb) & (dye <=(runs.bathy.xsb + distance));
                    mask = mask .* fillnan(dyemask,0); 
                    runs.eutrans.dye.trans(:,:,kk) = squeeze( ...
                            nansum(bsxfun(@times,avg1(mask .* cs_vel,2), ...
                            diff(runs.rgrid.z_r(:,runs.eutrans.ix(kk),1))'),2));
                        
                    runs.eutrans.dye.Itrans(:,kk) = squeeze(nansum(runs.eutrans.dye.trans(:,:,kk) .* dx,1))';
                end
            end

            %% plot transport
            
            figure;
            subplot(6,1,[1 2])
            plot(runs.rgrid.ocean_time(t0:end)/86400,runs.eutrans.Itrans/1e6);
            hold on
            %plot(runs.rgrid.ocean_time(t0:end)/86400,runs.eutrans.dye.Itrans/1e6,'--');
            limx = xlim;
            legend(num2str(runs.eutrans.h'),'Location','NorthWest');
            ylabel('Eulerian Transport (Sv)');
            title(['Isobaths in legend | Z < ' num2str(ceil(runs.bathy.hsb)) ' m ' ...
                '| mean eddy center isobath = '  num2str(mean(hcen)) ' m ' ...
                '| mean eddy edge isobath = ' num2str(mean(hedge)) 'm']);
            beautify;
            subplot(6,1,[3 4 5])
            plot(runs.rgrid.ocean_time/86400,runs.eutrans.dye.Itrans/1e6,'-');
            limx = xlim;
            legend(num2str(runs.eutrans.h'),'Location','NorthWest');
            ylabel('Dye Transport (Sv)');
            ylim([-0.05 0.3]); liney(0.1,[])
            beautify;
            subplot(6,1,6)
            [ax,~,~] = plotyy(runs.eddy.t,runs.eddy.prox/1000,runs.eddy.t,runs.eddy.hcen);
            set(ax(1),'XLim',limx);set(ax(2),'XLim',limx);
            set(ax(1),'XTickLabel',[]); axes(ax(2));
            set(get(ax(1),'ylabel'),'String','Proximity (km)');
            set(get(ax(2),'ylabel'),'String','h @ center of eddy');
            xlabel('Time (days)');
            
            % throw out locations where dye trans is pretty much zero to
            % make plot cleaner
            arr = [1:length(loc)];
            for kk=1:length(loc)
                if median(runs.eutrans.dye.Itrans(:,kk)) < 1
                    arr(arr == kk) = [];
                end
            end
            figure
            plot(runs.rgrid.ocean_time(t0:end)/86400,(runs.eutrans.Itrans(:,arr) - runs.eutrans.dye.Itrans(:,arr)) ...
                ./ runs.eutrans.dye.Itrans(:,arr) * 100);
            ylim([-100 700]); liney(0);
            legend(num2str(runs.eutrans.h(:,arr)'),'Location','NorthWest');
            title('percentage over-estimation = (eulerian - dye)/ dye');
            beautify;
            
            %% normalized transport plot
            %xmat = bsxfun(@minus,repmat(runs.rgrid.xr(:,1)/1000,[1 length(runs.rgrid.ocean_time)]), ...
            %                     runs.eddy.cx/1000);
            %tmat = repmat(runs.rgrid.ocean_time'/86400,[size(xmat,1) 1]);
            %plot(xmat,ntrans); linex(0)
            %disp_plot(runs.eutrans.dye.trans(:,:,4),xmat,runs.rgrid.ocean_time);
            
            time = runs.rgrid.ocean_time/86400;
            % normalize by max.
            mtrans = max(abs(runs.eutrans.dye.trans),[],1);
            ntrans = bsxfun(@rdivide,runs.eutrans.dye.trans, mtrans);
            mtrans = squeeze(mtrans);
            
            scrsz = get(0, 'ScreenSize'); 
            figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
            for kk = 1:size(runs.eutrans.Itrans,2)
                % figure out eddy edges at latitude of transport calculation
                emask = fillnan((bsxfun(@times, ...
                    squeeze(abs(diff(runs.eddy.mask(:,runs.eutrans.ix(kk),:),1))), ...
                    [1:size(runs.eddy.mask,1)-1]')'),0)';
                left = nanmin(emask); right = nanmax(emask);
                tmask = cut_nan(time' .* fillnan(~isnan(left),0));
                cmask = cut_nan(runs.eddy.cx/1000 .* fillnan(~isnan(left),0));
                
                clf;
                set(gcf,'Renderer','painters')
                subplot(1,5,[1 2 3]);
                hold on
                for ii=1:size(runs.rgrid.ocean_time)
                   plot(runs.rgrid.xr(:,1)/1000 - runs.eddy.cx(ii)/1000, ...
                       ntrans(:,ii,kk) + time(ii));
                end
                xlim([-200 50]);ylim([40 90]);
                %plot(runs.eddy.ee/1000 - runs.eddy.cx/1000,time,'r*');
                %plot(runs.eddy.we/1000 - runs.eddy.cx/1000,time,'r*');
                plot(runs.rgrid.xr(cut_nan(left),1)'/1000 - cmask,tmask,'r*');
                plot(runs.rgrid.xr(cut_nan(right),1)'/1000 - cmask,tmask,'k*');
                
                linex(0,'eddy center'); linex(-75);
                ylabel('Time (days)'); xlabel('X - X_{center}');
                title(['Normalized Transport (m^2/s) across '  ...
                    num2str(runs.eutrans.h(kk)) 'm isobath | red dots = edges']);
                beautify([14 14 16]);
                subplot(154)
                hold on
                if kk ~=1
                    plot(mtrans(:,1:kk-1),time,'Color',0.75*[1 1 1]);
                end
                if kk ~= size(runs.eutrans.Itrans,2)
                    plot(mtrans(:,kk+1:end),time,'Color',0.75*[1 1 1]);
                end
                plot(mtrans(:,kk),time,'b');
                xlabel('Max. Transport (m^2/day)');
                ylim([40 90]);xlim([0 40]);
                beautify([14 14 16]);

                subplot(155)
                hold on
                if kk ~=1
                    plot(runs.eutrans.dye.Itrans(:,1:kk-1)/1e6,time,'Color',0.75*[1 1 1]);
                end
                if kk ~= size(runs.eutrans.Itrans,2)
                    plot(runs.eutrans.dye.Itrans(:,kk+1:end)/1e6,time,'Color',0.75*[1 1 1]);
                end
                plot(runs.eutrans.dye.Itrans(:,kk)/1e6,time,'b'); 
                ylim([40 90]);xlabel('Total Transport (Sv)');
                title(sprintf('Max transport = %.2f Sv',(max(runs.eutrans.dye.Itrans(:,kk)/1e6))));
                xlim([-0.2 0.2]); linex(0);
                
                export_fig(sprintf('images/transport/%04d.png',runs.eutrans.h(kk)));
            end
            
            %% plotting tests
%            figure;
%             clim = [runs.bathy.xsb/1000 runs.bathy.xsb/1000+distance/1000];
%             rrfac = 7;
%             for ind = 1:size(runs.eddy.mask,3)
%                 clf
%                 subplot(211)
%                 pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.dye(:,:,ind)/1000); hold on
%                 dxi = 7; dyi = 7;
%                 if ~isempty(runs.usurf) && ~isempty(runs.vsurf)
%                     hq = quiver(runs.eddy.xr(1:dxi:end,1:dyi:end)/1000,runs.eddy.yr(1:dxi:end,1:dyi:end)/1000, ...
%                         runs.usurf(1:dxi:end,1:dyi:end,ind),runs.vsurf(1:dxi:end,1:dyi:end,ind));
%                 end
%                 caxis(clim);
%                 hold on
%                 title(['t = ' num2str(runs.rgrid.ocean_time(ind)/86400) ' days']);
%                 [~,hh] = contour(runs.eddy.xr/1000, runs.eddy.yr/1000,runs.eddy.mask(:,:,ind),1,'k');
%                 set(hh,'LineWidth',2);
%                 [~,hz] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,ind),5,'k');
%                 plot(runs.rgrid.xr/1000,runs.eutrans.trans(:,ind)*10);
%                 xlim([0 max(runs.rgrid.xr(:))/1000])
%                 linex(runs.eddy.we(ind)/1000); liney(runs.bathy.xsb/1000,'shelfbreak','b');
%                 linex(runs.eddy.we(ind)/1000 - rrfac*rr/1000,[num2str(rrfac) ' * RR']);
%                 linex(runs.eddy.we(ind)/1000 - 50,'center - 50 km');
%                 
%                 subplot(212)
%                 plot(runs.rgrid.xr/1000,runs.eutrans.trans(:,ind)); 
%                 title('Transport (m^2/sec)');
%                 ylim([floor(min(runs.eutrans.trans(:))) ceil(max(runs.eutrans.trans(:)))]);
%                 xlim([0 max(runs.rgrid.xr(:))/1000]);linex(runs.eddy.we(ind)/1000);
%                 
%                 linex(runs.eddy.we(ind)/1000 - rrfac*rr/1000,[num2str(rrfac) ' * RR']);
%                 liney(0);linex(runs.eddy.we(ind)/1000 - 50,'center - 50 km');
%                 pause
%             end
            
        end
        
        function [] = eddyevol(runs)
            eddy = runs.eddy;
            ii = 1; colors(1) = 'b';
            aa = 5; bb = aa*2;
            figure;
            subplot(aa,2,[1:2:bb]); hold on
            pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);
            xlabel('X (km)'); ylabel('Y (km)');
            plot(eddy.cx/1000,eddy.cy/1000,'Color',colors(ii,:),'LineWidth',2);
            colorbar
            if runs.bathy.axis == 'x'
                plot(eddy.we/1000,eddy.cy/1000,'k');
            else
                plot(eddy.cx/1000,eddy.se/1000,'k');
            end
            axis image; axis tight
            subplot(aa,2,2); hold on
            plot(eddy.t,eddy.amp,'Color',colors(ii,:));
            ylabel('amplitude (m)');
            subplot(aa,2,4); hold on
            plot(eddy.t,eddy.dia/1000,'Color',colors(ii,:));
            ylabel('diameter (km)');
            subplot(aa,2,6); hold on
            plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
            ylabel('x - center (km)');
            subplot(aa,2,8); hold on
            plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
            ylabel('y - center (km)');
            subplot(aa,2,10); hold on
            plot(eddy.t,eddy.Lz2,'Color',colors(ii,:));
            ylabel('vertical scale (m)');
            xlabel('time (days)');
        end
        
    end
end