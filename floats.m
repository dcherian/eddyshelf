classdef floats < handle
    properties
        x; y; z; time; age;
        fac; type
        temp; salt
        hitLand; hitBottom
        init; N = 0;
        comment = ['init = (x,y,z,t) = initial location, release time in meters, seconds | ' ...
                   'mom1 = first moment | disp = cloud dispersion | kur = kurtosis | ' ...
                   'N = number of floats at a time instant.'];
        disp; mom1; kur;
    end
    methods
        % reads data
        function [floats] = floats(type,file,rgrid)
            if strcmpi(type,'roms')
                floats.x = ncread(file,'x')';
                floats.y = ncread(file,'y')';
                floats.z = ncread(file,'depth')';
                floats.time = ncread(file,'ocean_time');
                floats.temp = ncread(file,'temp')';
                floats.salt = ncread(file,'salt')';
                floats.type = 'roms';
            end

            if strcmpi(type,'ltrans')
                floats.y = ncread(file,'lat')';
                floats.x = ncread(file,'lon')';
                floats.z = ncread(file,'depth')';
                floats.age = ncread(file,'age')'; 
                floats.time = ncread(file,'model_time');
                floats.temp = ncread(file,'temperature')';
                floats.salt = ncread(file,'salinity')';
                floats.hitLand = ncread(file,'hitLand')';
                floats.hitBottom = ncread(file,'hitBottom')';
                mask = zeros(size(floats.x));
                if max(floats.hitLand(:) ~= zeros(size(floats.hitLand(:))))
                    nf = length(cut_nan(fillnan(floats.hitLand,0)));
                    warning([num2str(nf) ' floats have hit land.']);
                    mask = (mask | cumsum(floats.hitBottom,1));
                end
                if max(floats.hitBottom(:) ~= zeros(size(floats.hitBottom(:))))
                    nf = length(cut_nan(fillnan(floats.hitBottom,0)));
                    warning([num2str(nf) ' floats have hit bottom.']);
                    % remove record after float has hit bottom i.e., fill
                    % with NaNs
                    mask = (mask | cumsum(floats.hitBottom,1));
                end
                floats.type = 'ltrans';
                mask = (mask | bsxfun(@eq,floats.x,floats.x(1,:))); 
                floats.x(mask) = nan; floats.y(mask) = nan; floats.z(mask) = nan;
                floats.temp(mask) = nan; floats.salt(mask) = nan;
            end

            if strcmpi(type,'tracmass')
                % ASSUMES REGULAR GRID
                % from the manual
                % The trajectories are stored in <outDataDir>/<outDataFile>_run.asc, which has to be
                % specified in <project>_run.in
                % The trajectories are stored in colums of
                %
                %              ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens
                %
                % where
                % ntrac is the trajectory number
                % niter is the TRACMASS code iteration (only important for TRACMASS modellers)
                % x1 is the zoonal position of the trajectory particle - INDEX
                % y1 is the meridional position of the trajectory particle - INDEX
                % z1 is the vertical position of the trajectory particle - INDEX
                % tt the time of the trajectory particle (in days)
                % t0 the initial time of the trajectory particle
                % subvol is the the "volume" transport in m3/s of the trajectory
                % temp is the temperature of the trajectory particle
                % salt is the salinity/specific humidity of the trajectory particle
                % dens is the density of the trajectory particle
                               tic 
                   [ntrac,~,ix,iy,iz,tt,t0,subvol,temp,salt,~] = ...
                                textread(fname,'%d%d%f%f%f%f%f%f%f%f%s');
                   toc
                   disp('Finished opening file. Now processing for unique records');
                   tic         
                %   data = sortrows(data,1); % sort according to trac number
                   floats.time = unique(tt);

                   xr = rgrid.x_rho(1,:)';
                   yr = rgrid.y_rho(:,1);

                   dx = xr(2)-xr(1); dy = yr(2)-yr(1);
                   cpb = progressbar();
                   for i = 1:length(unique(ntrac)) % ith drifter
                      %k = 1
                       j = find(ntrac == i);
                      %for j=1:length(ntrac)    
                        %if ntrac(j) == i
                            floats.t0(i) = t0(j(1));
                            fx = floor(ix(j)); fy = floor(iy(j)); fz = floor(iz(j));
                            cz = ceil(iz(j));
                            fy(fy == 0) = 1; fx(fx == 0) = 1;
                            % not all floats start at t=0
                            dt = (floats.t0(i)-floats.time(1))./ (floats.time(2)-floats.time(1));
                            k = 1:length(j);
                            floats.x(k+dt,i) = xr(fx) + (ix(j)-fx) * dx;
                            floats.y(k+dt,i) = yr(fy) + (iy(j)-fy) * dy;
                            for kk = 1:length(fz)
                                try
                                    if fz(kk) == 0, z1 = NaN; end
                                    z1 = rgrid.z_r(fz(kk),fy(kk),fx(kk));
                                catch ME
                                    disp('!!!')
                                end
                                z2 = rgrid.z_r(cz(kk),fy(kk),fx(kk));
                                dz = z2 - z1;
                                floats.z(kk+dt,i) = z1 + (iz(j(kk))-fz(kk)) * dz;
                            end
                            floats.iz(k+dt,i) = iz(j);
                            floats.temp(k+dt,i) = temp(j);
                            floats.salt(k+dt,i) = salt(j);
                            floats.t(k+dt,i)    = tt(j);
                            floats.subvol(k+dt,i) = subvol(j);
                            progressbarupdate(cpb,i/length(unique(ntrac)) * 100);
                        %    k=k+1;
                        %end
                      %end
                   end
                   cpb.stop();

                   for ii = 1:size(floats.x,2)
                        clear ind
                        ind = find(floats.x(:,ii) > 0);
                        ind = ind(1);

                        floats.init(ii,:) = [floats.x(ind,ii) floats.y(ind,ii) floats.z(ind,ii) floats.time(ind)];
                   end

                   % fill with nans
                   names = fieldnames(floats);
                   for ii=1:length(names)
                      floats.(names{ii}) = fillnan(floats.(names{ii}),0); 
                   end

                   floats.fac = 1; % outputs at model output
                   floats.type = 'tracmass';
                   toc
                   disp(' read float data ');
            end
            floats.N = sum(repnan(floats.x,0)~=0,2);
            floats.comment = ['init = (x,y,z,t) = initial location, release time in meters, seconds | ' ...
                'fac = number float timesteps per ROMS output timestep'];
            
            % calculate fac
            dtroms = rgrid.ocean_time(2)-rgrid.ocean_time(1);
            dtltr  = floats.time(2,1)-floats.time(1,1);
            floats.fac = dtroms/dtltr;
            
            % initial locations and seeding time
            for ii = 1:size(floats.x,2)
                clear ind
                ind = find(floats.x(:,ii) > 0);
                if isempty(ind)
                    ind = 1;
                else
                    ind = ind(1);
                end

                floats.init(ii,:) = [floats.x(ind,ii) floats.y(ind,ii) floats.z(ind,ii) floats.time(ind)];
            end
        end % read_floats   
        
        % plots displacements
        function [] = plot_displacements(floats)
            fig; 
            subplot(131);
            plot(floats.time/86400,floats.x/1000); 
            ylim([0 300]); xlim([0 200]);
            title('x')
            subplot(132);
            plot(floats.time/86400,floats.y/1000); 
            ylim([0 300]); xlim([0 200]);
            title('y')
            subplot(133);
            plot(floats.time/86400,floats.z); 
            ylim([-300 0]); xlim([0 200]);
            title('z')
            suplabel(upper(floats.type),'t');
        end
        
        % calculate lagrangian statistics
        function [] = stats(floats)
            T = size(floats.x,1);
            %first moment
            floats.mom1 = bsxfun(@times, ...
                         [nansum(bsxfun(@minus,floats.x,floats.init(:,1)')/1000,2) ...
                          nansum(bsxfun(@minus,floats.y,floats.init(:,2)')/1000,2) ...
                          nansum(bsxfun(@minus,floats.z,floats.init(:,3)')     ,2)] ...
                          ,1./floats.N);
            tic;  
            s = zeros(T,3);
            for i=1:size(floats.x,2)
                for j=1:size(floats.x,2)
                    if i~=j
                        s(:,1) = s(:,1) + repnan((floats.x(:,i) - floats.x(:,j)),0).^2/1e6;
                        s(:,2) = s(:,2) + repnan((floats.y(:,i) - floats.y(:,j)),0).^2/1e6;
                        s(:,3) = s(:,3) + repnan((floats.z(:,i) - floats.z(:,j)),0).^2;
                    end
                end
            end
            toc;
            floats.disp = bsxfun(@times,s,1./(2*floats.N.*(floats.N-1)));
            
            floats.kur = nan(T,3);
            floats.kur(:,1) = nansum(bsxfun(@minus,floats.x,floats.mom1(:,1)).^4,2)./ ...
                             (nansum(bsxfun(@minus,floats.x,floats.mom1(:,1)).^2,2)).^2;
            floats.kur(:,2) = nansum(bsxfun(@minus,floats.y,floats.mom1(:,2)).^4,2)./ ...
                             (nansum(bsxfun(@minus,floats.y,floats.mom1(:,2)).^2,2)).^2;
            floats.kur(:,3) = nansum(bsxfun(@minus,floats.z,floats.mom1(:,3)).^4,2)./ ...
                             (nansum(bsxfun(@minus,floats.z,floats.mom1(:,3)).^2,2)).^2;
        end
        
        % plot stats in subplots
        function [hfig] = plot_stats(floats,hfig)
            if isempty(floats.kur), floats.stats; end
            if ~exist('hfig','var'), hfig = gcf; end
            if strcmp(floats.type,'ltrans')
                fmt = '--';
            else
                fmt = '-';
            end
            figure(hfig);
            subplot(311)
            plot(floats.time/86400,floats.mom1,fmt); hold on;
            xlabel('time (days)'); ylabel('first moment');
            legend('x (km)','y (km)','z (m)','Location','NorthWest');
            title(' dashed = LTRANS, line = ROMS');
            subplot(312)
            plot(floats.time/86400,floats.disp,fmt); hold on;
            xlabel('time (days)'); ylabel(' Cloud Dispersion (mean sq. pair separation)');
            legend('x (km^2)','y (km^2)','z (m^2)','Location','NorthWest');
            subplot(313)
            plot(floats.time/86400,floats.kur,fmt); hold on;
            xlabel('time (days)'); ylabel('Kurtosis');
            legend('x','y','z','Location','NorthWest');
        end
        
        % animates with zeta plot
        function [] = animate(floats,rgrid,zeta)

            cmap = flipud(cbrewer('div', 'BrBG', 32));

            figure;
            tfilt = cut_nan(fillnan(floats.init(:,4),0));
            ind0 = find_approx(rgrid.ocean_time, tfilt(1),1);

            for i=ind0:size(zeta,3)
                cla
                contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,i)');
                shading flat; axis image
                colormap(cmap); 
                caxis([min(zeta(:)) max(zeta(:))]);
                hold on
                [C,hc] = contour(rgrid.x_rho./1000,rgrid.y_rho./1000,rgrid.h,[114 500 750 1100],'k');
                clabel(C,hc);
                plot(floats.init(:,1)/1000,floats.init(:,2)/1000,'x','MarkerSize',12);
                floats.fac = floor(floats.fac);
                nn = find_approx(floats.time,rgrid.ocean_time(i),1);
                plot(floats.x(nn,:)/1000,floats.y(nn,:)/1000,'k.','MarkerSize',10);
                title(['t = ' num2str(floats.time(1)+i) ' days']);
                pause(0.01)
            end
        end
    end
%% tracks at timestep movie
% dir = 'E:\Work\eddyshelf\runs\topoeddy\runteb-04-hires-6\';
% fname = [dir 'ocean_avg.nc'];
% 
% rgrid = roms_get_grid(fname,fname,1,1);
% redo = 1;
% if ~exist('trac','var') || redo == 1
%     trac = tracmass_read([dir '\eddytest_run.asc'],rgrid);
% end
% save([dir '/trac.mat'],'trac');
% 
% %load trac_movie
% load([dir '\eddytrack.mat']);
% figure
% floats = roms;
% for ii=1:length(floats.time)
%     kk = find_approx(eddy.t,floats.time(ii));
%     clf
%     plot(floats.x(ii,:)/1000,floats.y(ii,:)/1000,'.','MarkerSize',12);
%     hold on;
%     plot(eddy.cx(kk)/1000,eddy.cy(kk)/1000,'r.','MarkerSize',16);
%     plot(eddy.cx/1000,eddy.cy/1000,'k');
%     [cc,hh] = contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,kk),1,'k');
%     set(hh,'LineWidth',2);
%     xlim([0 180]); ylim([0 300]);
%     xlabel('X (km)'); ylabel('Y (km)');
%     pause(0.01);
% end

end