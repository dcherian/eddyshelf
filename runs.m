classdef runs < handle
    properties
        % dir & file names
        name; dir; out_file; ltrans_file; flt_file; givenFile
        % data
        zeta; temp; usurf; vsurf;vorsurf;
        % dimensional and non-dimensional time
        time; ndtime;
        % barotropic vel (geostrophic)
        ubarg; vbarg;
        % dyes
        csdye; asdye; zdye; % cross-shore, along-shore, z dyes
        % grid & bathymetry
        rgrid; bathy
        % float data
        roms; ltrans;
        % eddy track data
        eddy; noeddy;
        % initial params
        params
        % transport
        eutrans;
        % make video?
        makeVideo; mm_instance;
        %
        comment = ['eddy.prox = distance of edge from shelfbreak in m | '...
                   'eutrans = eulerian transport estimate (structure).' ...
                   'eddy.trev = time at which eddy reverses direction (first)'];
    end
    methods
        % constructor
        function [runs] = runs(dir,reset)
            if ~exist('reset','var')
                reset = 0;
            end
            
            if isdir(dir)
                runs.dir = dir;
                runs.out_file = [dir '/ocean_avg.nc'];
                runs.givenFile = 0;
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
            runs.rgrid.z_uw = [];
            runs.rgrid.z_vw = [];
            runs.makeVideo = 0; % no videos by default.
            
            % make run-name
            ind1 = strfind(runs.dir,'/run');
            runs.name = runs.dir(ind1+4:end);
            if runs.name(end) == '/'
                runs.name(end) = [];
            end
            
            % params & bathy
            runs.params = read_params_from_ini(runs.dir);
            runs.bathy = runs.params.bathy;
            
            % fill bathy
            [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = find_shelfbreak(runs.out_file);
            [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = find_shelfbreak(runs.out_file,'slope');
            runs.bathy.h = runs.rgrid.h';
            
            if ~runs.givenFile
                runs.zeta = roms_read_data(dir,'zeta');
                runs.time = roms_read_data(dir,'ocean_time');
                try
                    runs.dye  = roms_read_data(dir,'dye_02');
                catch ME
                end
            else
                runs.zeta = double(ncread(runs.out_file,'zeta'));
                runs.time = double(ncread(runs.out_file,'ocean_time'));
                try
                    runs.csdye  = squeeze(double(ncread(runs.out_file,'dye_02', ...
                    [1 1 runs.rgrid.N 1],[Inf Inf 1 Inf])));
                catch ME
                end
            end
            
            if exist(runs.flt_file,'file')
                runs.roms = floats('roms',runs.flt_file,runs.rgrid);
            end
            
            if ~exist([dir '/eddytrack.mat'],'file') || reset == 1 ...
                    %|| ~exist('runs.eddy.cvx','var')
                try
                    runs.eddy = track_eddy(dir);
                    runs.noeddy = 0;
                catch ME
                    disp(ME.message);
                    disp('Couldn''t run track_eddy.m');
                    runs.noeddy = 1;
                end
            else
                if strfind(runs.out_file,'_004.nc')
                    edd = load([dir '/eddytrack_004.mat'],'eddy');
                else
                    edd = load([dir '/eddytrack.mat'],'eddy');
                end
                runs.eddy = edd.eddy;
                runs.noeddy = 0;
            end
            
            if ~runs.noeddy
               if isfield(runs.eddy,'cvx')
                   if runs.eddy.cvx(1) == 0 || runs.eddy.cvy(1) == 0
                    runs.eddy.cvx(1) = NaN;
                    runs.eddy.cvy(1) = NaN;
                   end
               end
            
                if runs.bathy.axis == 'y'
                    edge = runs.eddy.se;
                else
                    edge = runs.eddy.we;
                end
                % proximity to shelfbreak
                runs.eddy.prox = (edge-runs.bathy.xsb);
                % time of reversal
                runs.eddy.trev = runs.time(find(runs.eddy.cvx < 0,1,'first'));
            
                % water depth at eddy center
                h = runs.bathy.h(2:end-1,2:end-1);
                ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx);
                iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.my);
                runs.eddy.hcen = h(sub2ind(size(runs.eddy.xr),ix,iy))';
                % non-dimensionalized time
                %runs.ndtime = (runs.eddy.cx - runs.eddy.cx(1))./ ...
                %        (runs.params.bg.ubt -  ...
                %        runs.params.phys.beta/2*(runs.eddy.dia/2),^2);
            end
            if exist(runs.ltrans_file,'file')
                runs.ltrans = floats('ltrans',runs.ltrans_file,runs.rgrid);
            end
        end
       
        function [] = info(runs)
            roms_info(runs.dir);
        end
        
       %% floats
        function [] = compare_floats(runs)
            ltransc = floats('ltrans',[runs.dir '/ltrans-compare.nc'],runs.rgrid);
            runs.roms.plot_stats;
            ltransc.plot_stats;
        end
        
        % create initial seed file for ltrans
        function [] = ltrans_create(runs)
            ltrans_create(runs.rgrid,runs.zeta,runs.eddy);
        end
        
        % create ltrans init file from roms out
        function [] = ltrans_create_from_roms(runs)
            ltrans_create_from_roms('ltrans_init_compare.txt',runs.flt_file,runs.rgrid);
        end
        
       %% analysis
       
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
                if ~isempty(runs.csdye)
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
            
            time = runs.time/86400;
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
            
            tind = find(runs.time == runs.eddy.trev);
            
            figure;
            subplot(aa,2,[1:2:bb-2*2]); hold on
            pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);
            xlabel('X (km)'); ylabel('Y (km)');
            plot(eddy.cx/1000,eddy.cy/1000,'Color',colors(ii,:),'LineWidth',2);
            plot(eddy.cx(tind)/1000,eddy.cy(tind)/1000,'*','Color',colors(ii,:), ...
                    'MarkerSize',12);
            colorbar
            if runs.bathy.axis == 'x'
                plot(eddy.we/1000,eddy.cy/1000,'k');
            else
                plot(eddy.cx/1000,eddy.se/1000,'k');
            end
            axis image; axis tight
            title('Bathy + eddy track');
            subplot(aa,2,2); hold on
            plot(eddy.t,eddy.amp,'Color',colors(ii,:));
            ylabel('amplitude (m)');
            linex(tind);
            subplot(aa,2,4); hold on
            plot(eddy.t,eddy.dia/1000,'Color',colors(ii,:));
            ylabel('diameter (km)');
            linex(tind);
            subplot(aa,2,6); hold on
            plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
            ylabel('x - center (km)');
            linex(tind);
            subplot(aa,2,8); hold on
            plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
            ylabel('y - center (km)');
            linex(tind);
            subplot(aa,2,10); hold on
            plot(eddy.t,eddy.prox/1000,'Color',colors(ii,:));
            liney(min(eddy.prox/1000));
            ylabel('Proximity (km)');
            xlabel('time (days)');
            linex(tind);
            
            subplot(aa,2,[7 9] ); hold on
            plot(eddy.t,eddy.hcen/2,'b');
            plot(eddy.t,runs.params.phys.f0 / sqrt(runs.params.phys.N2) * runs.eddy.dia,'r');
            plot(eddy.t,eddy.Lz2,'k');
            legend('H_{center}/2','f*L/N','vertical scale','Location','SouthWest');
            xlabel('Time (days)');
            ylabel('(m)');
            linex(tind);
            suplabel(runs.dir,'t');
        end
        
        function [] = compare_plot(runs,num)
            eddy = runs.eddy;
            ii = num;
            colors = distinguishable_colors(10);
            aa = 6; bb = aa*2;
            tloc = [50 60 70];
            subplot(aa,2,[1:2:bb-2*2]); hold on
            %pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);            colorbar
            xlabel('X (km)'); ylabel('Y (km)');
            plot(eddy.cx/1000,eddy.cy/1000,'Color',colors(ii,:),'LineWidth',2);
            plot(eddy.cx(tloc)/1000,eddy.cy(tloc)/1000,'*','MarkerSize',12,'Color',colors(ii,:),'LineWidth',2);
            %if runs.bathy.axis == 'x'
            %    plot(eddy.we/1000,eddy.cy/1000,'Color',colors(ii,:),'LineStyle','--');
            %else
            %    plot(eddy.cx/1000,eddy.se/1000,'Color',colors(ii,:),'LineStyle','--');
            %end
            axis image; axis tight
            subplot(aa,2,2); hold on
            plot(eddy.t,eddy.amp,'Color',colors(ii,:));
            ylabel('amplitude (m)');
            subplot(aa,2,4); hold on
            plot(eddy.t,eddy.dia/1000,'Color',colors(ii,:));
            ylabel('diameter (km)');
            subplot(aa,2,6); hold on
            %plot(eddy.t,eddy.mvx,'Color',colors(ii,:));
            %ylabel('x - vel (km/day)');
            %ylim([-5 5]);
            plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
            ylabel('x - center (km)');
            subplot(aa,2,8); hold on
            %plot(eddy.t,eddy.mvy,'Color',colors(ii,:));
            %ylabel('y - vel (km/day)');
            %ylim([-5 5]);
            plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
            ylabel('y - center (km)');
            subplot(aa,2,10); hold on
            plot(eddy.t,eddy.Lz2,'Color',colors(ii,:));
            ylabel('vertical scale (m)');
            %xlabel('time (days)');
            subplot(aa,2,12); hold on
            plot(eddy.t,eddy.prox/1000,'Color',colors(ii,:));
            xlabel('time (days)');
            ylabel('Proximity (km)');
            
            subplot(aa,2,[9 11] ); hold on
            plot(eddy.t,eddy.hcen,'Color',colors(ii,:));
%             plot(eddy.t,runs.params.phys.f0 / sqrt(runs.params.phys.N2) * runs.eddy.dia,'Color',colors(ii,:),'LineStyle','--');
%             legend('H_{center}','f/N*dia');
            xlabel('Time (days)');
            ylabel('H_{center}(m)');
        end
        
        function [] = tracer_budget(runs)
            tracer = roms_read_data(runs.out_file,'dye_02');
            s = size(tracer);
            %Itracer = domain_integrate(tracer, ...
            %                runs.rgrid.xr,runs.rgrid.yr,runs.rgrid.zr);
            
            clear N
            
            lim = linspace(min(min(tracer(:,:,end,1))),max(max(tracer(:,:,end,1))),90);
            tracer = reshape(tracer,[s(1)*s(2)*s(3) s(4)]);
            
            for i=1:s(4)
                [N(:,i),bins] = histc(tracer(:,i),lim);
            end
            
            colors = flipud(repmat(linspace(0,0.9,s(4))',[1 3]));
            figure
            set(gca,'ColorOrder',colors); hold all
            plot(lim/1000,N);
            set(gcf,'Colormap',colors); 
            hcbar = colorbar;
            tlab = ceil(runs.rgrid.ocean_time(get(hcbar,'YTick'))/86400);
            set(hcbar,'YTickLabel',num2str(tlab))
            xlabel('Cross-shore axis (km)');
            ylabel('Count');
            cblabel('Time (days)');
            beautify ([14 14 16]);
            
        end
        
        % read surface velocities for animate_pt & surf vorticity plot
        function [] = read_velsurf(runs)
            disp('Reading surface velocity fields...');
            tic;
            start = [1 1 runs.rgrid.N 1];
            count = [Inf Inf 1 Inf];
            stride = [1 1 1 1];
            
            if runs.givenFile
                runs.usurf = double(squeeze(ncread(runs.out_file, ....
                    'u',start,count,stride)));
            else
                runs.usurf = roms_read_data(runs.dir,'u' ...
                    ,start,count,stride);
            end
            runs.usurf = avg1(runs.usurf(:,2:end-1,:),1);
            toc;
            if runs.givenFile
                runs.vsurf = double(squeeze(ncread(runs.out_file, ....
                    'v',start,count,stride)));
            else
                runs.vsurf = roms_read_data(runs.dir,'v',start,count,stride);
            end
            runs.vsurf = avg1(runs.vsurf(2:end-1,:,:),2);
            toc;
        end
        
        % calculate surface vorticity field
        function [] = calc_vorsurf(runs)
            if isempty(runs.usurf) || isempty(runs.vsurf)
                runs.read_velsurf;
            end
            
            if isempty(runs.vorsurf)
                runs.vorsurf = avg1(bsxfun(@rdivide,diff(runs.vsurf,1,1), ...
                                        diff(runs.rgrid.xr(2:end-1,2:end-1),1,1)),2) ...
                                - avg1(bsxfun(@rdivide,diff(runs.usurf,1,2), ...
                                        diff(runs.rgrid.yr(2:end-1,2:end-1),1,2)),1);
                                    
                runs.rgrid.xvor = avg1(avg1(runs.rgrid.xr(2:end-1,2:end-1),1),2);
                runs.rgrid.yvor = avg1(avg1(runs.rgrid.yr(2:end-1,2:end-1),1),2);
            end
        end
        
        % calculate geostrophically balanced barotropic velocities
        function [] = calc_ubarg(runs)
            runs.ubarg = -1 * 9.81 .* bsxfun(@rdivide,diff(runs.zeta,1,2), ...
                                avg1(runs.rgrid.f',2).*diff(runs.rgrid.yr,1,2));
        
            runs.vbarg =      9.81 .* bsxfun(@rdivide,diff(runs.zeta,1,1), ...
                                avg1(runs.rgrid.f',1).*diff(runs.rgrid.xr,1,1));
        end
                
       %% animation functions
        
        function [] = animate_zeta(runs)
            runs.video_init('zeta');
            
            figure;
            ii=1;
            hz = runs.plot_zeta('pcolor',ii);
            ax = gca;
            hold on
            colorbar; freezeColors;
            hbathy = runs.plot_bathy('contour','k');
            he = runs.plot_eddy_contour('contour',ii);
            ht = title([' SSH (m) | ' num2str(runs.time(1)/86400)  ' days']);
            xlabel('X (km)');ylabel('Y (km)');
            axis image;
            maximize(gcf); pause(0.2);  
            beautify([16 16 18]);
            runs.video_update();
            for ii = 2:size(runs.zeta,3)
                runs.update_zeta(hz,ii);
                runs.update_eddy_contour(he,ii);
                set(ht,'String',[' SSH (m) | ' num2str(runs.time(ii)/86400) ' days']);
                runs.video_update();
                pause(0.03);
            end
            runs.video_write();
        end
        
        function [] = animate_vorsurf(runs)         
            if isempty(runs.vorsurf)
                runs.calc_vorsurf();
            end
           
            tt = 1;
            vormax = max(abs(runs.vorsurf(:)))/4;
            levels = linspace(-vormax,vormax,20);
            [~,hh] = contourf(runs.rgrid.xvor/1000,runs.rgrid.yvor/1000, ...
                runs.vorsurf(:,:,tt),levels);
            caxis([-1 1] * vormax); colorbar; 
            xlabel('X (km)'); ylabel('Y (km)');
            axis image;
            ht = title(['Surface vorticity @ t = ' num2str(tt) ' days']);
            for tt = 2:size(runs.vorsurf,3)
                set(hh,'ZData',runs.vorsurf(:,:,tt));
                set(ht,'String',['Surface vorticity @ t = ' num2str(tt) ' days']);
                pause(0.02);
            end
                        
        end
        
        function [] = animate_vor(runs,tind)
            if ~exist('tind','var')
                tind = [];
            end
            if ~exist([runs.dir '/ocean_vor.nc'],'file')
                dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
            end
            
            tt = 1;
            
            for tt=2:read_count(end)
                
            end
            
        end
        
        function [] = animate_center(runs)
            runs.video_init('center');
            eddy = runs.eddy;
            xvec = runs.rgrid.xr(:,1);
            yvec = runs.rgrid.yr(1,:)';
            
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
            
            disp('Reading data.');
            tic;
            if runs.bathy.axis == 'x'
                stride = [sxy 1 sz st];
                temper = roms_read_data(runs.dir,'temp',[1 iymin 1 t0], ...
                                  ceil([Inf iymax-iymin+1 Inf Inf]./stride), stride);
                              toc;
                strat  = roms_read_data(runs.dir,'temp',[Inf 1 1 1], ...
                                  ceil([1 1 Inf 1]./stride),stride);
                              toc
            else
                stride = [1 sxy sz st];
                temper = roms_read_data(runs.dir,'temp',[ixmin 1  1 t0], ...
                                ceil([ixmax-ixmin+1 Inf Inf Inf]./stride),stride);
                            toc;
                strat  = roms_read_data(runs.dir,'temp',[1 1 1 1], ...
                                ceil([1 Inf Inf 1]./stride),stride);               
                            toc;
            end

            temper = bsxfun(@minus,temper,permute(strat,[3 1 2]));
            
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
            caxis([-1 1]*max(mat2vec(abs(temper(ix-ixmin+1,:,:,1:end-10)))));
            h1 = liney(-eddy.Lz2(stride(4)),[],'b');
            ylim([-1500 0]);
            title('Cross-shore temperature anomaly - slice through eddy center');
            %h2 = liney(-eddy.Lz3(stride(4)),'3','k');
            maximize(gcf); pause(0.2);  
            beautify([16 16 18]);
            
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
            if isempty(runs.usurf) || isempty(runs.vsurf)
                runs.read_velsurf;
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
        
        function [] = animate_zslice(runs,varname,depth,tind)
            % process tind
            if ~exist('tind','var'), tind = []; end
            [~,tind,~,nt,stride] = roms_tindices(tind,Inf,length(runs.time));
            
            read_start = [1 1 1 tind(1)];
            read_count = [Inf Inf Inf nt];
            
            if strcmp(varname,'vor');
                grids = [runs.dir '/ocean_vor.nc'];
            else
                grids = runs.rgrid;
            end
            
            [grd.xax,grd.yax,grd.zax,~] = dc_roms_extract(grids,varname,{},1);
            datain= 0;
            if nt < 20
                tic; disp('Reading data...');
                data = roms_read_data(runs.dir,varname, ...
                        read_start,read_count,stride);
                datain = 1;
                var = nan([size(data,1) size(data,2) nt]);
                toc;
            end
            % read data
            for mmm = 1:nt
                
                if ~datain
                    disp(['reading & interpolating timestep ' num2str(mmm) '/' ...
                                num2str(nt)]);
                    data = roms_read_data(runs.dir,varname, ...
                            [read_start(1:3) read_start(4)+mmm-1], ...
                            [read_count(1:3) 1],stride);
                    if mmm == 1
                        var = nan([size(data,1) size(data,2) nt]);
                    end
                    var(:,:,mmm) = dc_roms_zslice_var(data,depth,grd);
                else
                    disp(['interpolating timestep ' num2str(mmm) '/' ...
                                num2str(nt)]);
                    var(:,:,mmm) = dc_roms_zslice_var(data(:,:,:,mmm),depth,grd);
                end
            end
            clear data
            
            % animate
            xax = grd.xax/1000; yax=  grd.yax/1000; clear grd;
            tt = 1;
            [hc] = pcolor(xax,yax,var(:,:,tt)); shading interp
            hold on
            he = runs.plot_eddy_contour('contour',tind(1) + tt-1);
            ht = title([varname ' | z = ' num2str(depth) ' m | t = ' ...
                num2str(runs.time(tind(1)+tt-1)/86400) ' days']);
            axis image;
            xlim([min(xax(:)) max(xax(:))]);
            ylim([min(yax(:)) max(yax(:))]);
            colorbar; caxis([min(var(:)) max(var(:))]);
            xlabel('X (km)'); ylabel('Y (km)');
            runs.plot_bathy('contour','k');
            pause();
            for tt=2:nt
                set(hc,'CData',var(:,:,tt));
                runs.update_eddy_contour(he,tind(1) + tt-1);
                set(ht,'String',[varname ' | z = ' num2str(depth) ' m | t = ' ...
                num2str(runs.time(tind(1)+tt-1)/86400) ' days']);
                pause(0.01);
            end
            
        end
        
       %% generic plotting functions
        function [hplot] = plot_zeta(runs,plottype,tt)
            if ~exist('tt','var'), tt = 1; end
            
            if strcmpi(plottype,'pcolor')
                hplot = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,tt));
                if runs.makeVideo
                    shading interp; 
                else
                    shading flat
                end
            else
                if strcmpi(plottype,'contourf')
                    hplot = contourf(runs.rgrid.xr/1000,runs.rgrid.yr/1000, ...
                        runs.zeta(:,:,tt));
                    shading flat
                end
            end
        end
        function update_zeta(runs,handle,tt)
            try
                set(handle,'CData',runs.zeta(:,:,tt));
            catch ME
                set(handle,'ZData',runs.zeta(:,:,tt));
            end
        end
        
        function [hplot] = plot_eddy_contour(runs,plottype,tt)
            try 
                mask = runs.eddy.vormask(:,:,tt);
            catch
                mask = runs.eddy.mask(:,:,tt);
            end
            [~,hplot] = contour(runs.eddy.xr/1000,runs.eddy.yr/1000, ...
                        mask(:,:,tt),'Color','k');
        end
        function update_eddy_contour(runs,handle,tt)
            try 
                mask = runs.eddy.vormask(:,:,tt);
            catch
                mask = runs.eddy.mask(:,:,tt);
            end
            set(handle,'ZData',mask);
        end
        
        function [hplot] = plot_bathy(runs,plottype,color)
            if ~exist('color','var'), color = 'w'; end
            if strcmpi(plottype,'contour')
                [cc,hplot] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000, ...
                                runs.rgrid.h',[200 500 1000 1500 2000],'k');
                clabel(cc,hplot);
                if runs.bathy.axis == 'y'
                    liney(runs.bathy.xsb/1000,'shelfbreak',color);
                else
                    linex(runs.bathy.xsb/1000,'shelfbreak',color);
                end
            end
        end
        
       %% video functions
        function [] = video_init(runs,filename)
            if runs.makeVideo
                runs.makeVideo
                runs.mm_instance = mm_setup('frameDir',['videos/' runs.name '-' filename]);
                runs.mm_instance.pixelSize = [1600 900];
                runs.mm_instance.outputFile = ['videos/' runs.name '-' filename '.mp4'];
                runs.mm_instance.ffmpegArgs = '-q:v 1 -g 1';
                runs.mm_instance.InputFrameRate = 5;
                runs.mm_instance.frameRate = 5;
%                 aviobj = VideoWriter('output','MPEG-4');
%                 open(aviobj);
            end
        end
        
        function [] = video_update(runs)
            if runs.makeVideo
                mm_addFrame(runs.mm_instance,gcf);
                %F = getframe(gcf);
                %writeVideo(aviobj,F);
            end
        end
        
        function [] = video_write(runs)
            if runs.makeVideo
               mm_render(runs.mm_instance);
               %close(aviobj);
            end
        end
        
        
        function [] = imageEffect(runs)
            dx = runs.rgrid.xr(2,1)-runs.rgrid.xr(1,1);
            dy = runs.rgrid.yr(1,2)-runs.rgrid.yr(1,1);
            % eddy vorticity
            if isempty(runs.vorsurf)
                runs.calc_vorsurf();
            end
            w = avg1(avg1(runs.eddy.mask,1),2).*runs.vorsurf;
            % circulation
            circ = squeeze(dx*dy * sum(sum(w,1),2));
            
            plot(runs.time/86400,circ);
            ylabel('Surface Circulation');
            xlabel('Time (days)');
        end
        
        
        
       
        
    end
end