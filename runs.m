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
            
            if ~runs.givenFile
                runs.zeta = roms_read_data(dir,'zeta');
                runs.dye  = roms_read_data(dir,'dye_02');
            else
                runs.zeta = double(ncread(runs.out_file,'zeta'));
                runs.dye  = squeeze(double(ncread(runs.out_file,'dye_02', ...
                    [1 1 runs.rgrid.N 1],[Inf Inf 1 Inf])));
            end
            
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
                edd = load([dir '/eddytrack.mat'],'eddy');
                runs.eddy = edd.eddy;
            end
            
            if exist(runs.ltrans_file,'file')
                runs.ltrans = floats('ltrans',runs.ltrans_file,runs.rgrid);
            end
            
            runs.params = read_params_from_ini(runs.dir);
            runs.bathy = runs.params.bathy;
            
            % fill bathy
            [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = find_shelfbreak(runs.out_file);
            [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = find_shelfbreak(runs.out_file,'slope');
            runs.bathy.h = runs.rgrid.h';
        end
        
        function [] = animate(runs)
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
            t0 = 1;
            if runs.params.bathy.axis == 'x'
                error(' not built for north-south isobaths');
            else
                loc = nanmean(runs.eddy.se(t0:end));
            end
            
            runs.eutrans.x = loc;
            runs.eutrans.ix = find_approx(runs.rgrid.yr(1,:),loc,1);
            runs.eutrans.h = runs.bathy.h(1,runs.eutrans.ix);
            % cross-shore velocity
            cs_vel = double(squeeze(ncread(runs.out_file,'v', ...
                [1 runs.eutrans.ix 1 t0],[Inf 1 Inf Inf])));
            
            mask = nan(size(cs_vel));
            iwest = vecfind(runs.eddy.xr(:,1),runs.eddy.we);
            for tt=1:size(cs_vel,3)
                mask(1:iwest(tt),:,tt) = 1;
            end
            
            runs.eutrans.trans = squeeze( ...
                 trapz(runs.rgrid.z_r(:,runs.eutrans.ix,1),mask .* cs_vel,2));
            dx = runs.rgrid.xr(2,1)-runs.rgrid.xr(1,1);
             
            runs.eutrans.Itrans = squeeze(nansum(runs.eutrans.trans .* dx,1));
            
            %% plotting tests
            rr = sqrt(runs.params.phys.N2)*runs.bathy.hsb/runs.rgrid.f(runs.bathy.isb,1);
            distance = 5*rr; % 5 times rossby radius
            clim = [runs.bathy.xsb/1000 runs.bathy.xsb/1000+distance/1000];
            rrfac = 7;
            for ind = 1:size(runs.eddy.mask,3)
                clf
                subplot(211)
                pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.dye(:,:,ind)/1000); hold on
                dxi = 7; dyi = 7;
                if ~isempty(runs.usurf) && ~isempty(runs.vsurf)
                    hq = quiver(runs.eddy.xr(1:dxi:end,1:dyi:end)/1000,runs.eddy.yr(1:dxi:end,1:dyi:end)/1000, ...
                        runs.usurf(1:dxi:end,1:dyi:end,ind),runs.vsurf(1:dxi:end,1:dyi:end,ind));
                end
                caxis(clim);
                hold on
                title(['t = ' num2str(runs.rgrid.ocean_time(ind)/86400) ' days']);
                [~,hh] = contour(runs.eddy.xr/1000, runs.eddy.yr/1000,runs.eddy.mask(:,:,ind),1,'k');
                set(hh,'LineWidth',2);
                [~,hz] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,ind),5,'k');
                plot(runs.rgrid.xr/1000,runs.eutrans.trans(:,ind)*10);
                xlim([0 max(runs.rgrid.xr(:))/1000])
                linex(runs.eddy.we(ind)/1000); liney(runs.eutrans.x/1000,'shelfbreak','b');
                linex(runs.eddy.we(ind)/1000 - rrfac*rr/1000,[num2str(rrfac) ' * RR']);
                linex(runs.eddy.we(ind)/1000 - 50,'center - 50 km');
                
                subplot(212)
                plot(runs.rgrid.xr/1000,runs.eutrans.trans(:,ind)); 
                title('Transport (m^2/sec)');
                ylim([floor(min(runs.eutrans.trans(:))) ceil(max(runs.eutrans.trans(:)))]);
                xlim([0 max(runs.rgrid.xr(:))/1000]);linex(runs.eddy.we(ind)/1000);
                
                linex(runs.eddy.we(ind)/1000 - rrfac*rr/1000,[num2str(rrfac) ' * RR']);
                liney(0);linex(runs.eddy.we(ind)/1000 - 50,'center - 50 km');
                pause
            end
            
        end
    end
end