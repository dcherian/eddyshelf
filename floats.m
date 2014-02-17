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
            disp('Reading float data.');
            tic;
            if strcmpi(type,'roms')
                try
                    floats.x = ncread(file,'lon')';
                catch ME
                    floats.x = ncread(file,'x');
                end
                try
                    floats.y = ncread(file,'lat')';
                catch ME
                    floats.y = ncread(file,'y');
                end
                floats.z = ncread(file,'depth')';
                floats.time = ncread(file,'ocean_time');
                floats.temp = ncread(file,'temp')';
                try
                    floats.salt = ncread(file,'salt')';
                catch
                    floats.salt = nan(size(floats.temp));
                end
                floats.type = 'roms';
                disp('Read data. Now processing.');
                toc;tic;
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
                disp('Read data. Now processing.');
                toc; tic;
                mask = zeros(size(floats.x));
                if max(floats.hitLand(:) ~= zeros(size(floats.hitLand(:))))
                    % hitLand is 1 at time the float hits land
                    nf = sum(floats.hitLand(:));
                    warning([num2str(nf) ' floats have hit land.']);
                    mask = (mask | cumsum(floats.hitLand,1));
                end
                if max(floats.hitBottom(:) ~= zeros(size(floats.hitBottom(:))))
                    maskbot = (cumsum(floats.hitBottom) >= 1);
                    nf = diff(maskbot);
                    nf = sum(nf(:));
                    warning([num2str(nf) ' floats have hit bottom.']);
                    % remove record after float has hit bottom i.e., fill
                    % with NaNs
                    mask = (mask | maskbot);
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
            
            % remove floats that are all NaN
            tmat = repmat(floats.time,[1 size(floats.x,2)]);
            nanmask = ~all(isnan(floats.x),1);
            floats.x = floats.x(:,nanmask);
            floats.y = floats.y(:,nanmask);
            floats.z = floats.z(:,nanmask);
            if strcmpi(type,'ltrans')
                floats.age = floats.age(:,nanmask);
                floats.hitLand = floats.hitLand(:,nanmask);
                floats.hitBottom = floats.hitBottom(:,nanmask);
            end
            floats.temp = floats.temp(:,nanmask);
            floats.salt = floats.salt(:,nanmask);
            
            % for ROMS floats remove those that have stuck to bottom / land
            % i.e., they go to (x,y,z) = (0,0,0);
            if strcmpi(type,'roms')
                mask2 = (floats.x == 0) & (floats.y == 0) & (floats.z == 0);
                floats.x(mask2) = NaN;
                floats.y(mask2) = NaN;
                floats.z(mask2) = NaN;
            end
            
            %
            floats.N = sum(repnan(floats.x,0)~=0,2);
            floats.comment = ['init = (x,y,z,t) = initial location, release time in meters, seconds | ' ...
                'fac = number float timesteps per ROMS output timestep'];
            
            % calculate fac
            try
                dtroms = rgrid.ocean_time(2)-rgrid.ocean_time(1);
                dtltr  = floats.time(2,1)-floats.time(1,1);
                floats.fac = dtroms/dtltr;
            catch ME
                 warning('havent calculated fac');
            end
            % initial locations and seeding time
            initmask = find([zeros([1 size(floats.x,2)]); abs(diff((cumsum(repnan(floats.x,0)) >= 1)))] == 1);
            floats.init = [floats.x(initmask) floats.y(initmask) floats.z(initmask) tmat(initmask)];
    %             tic;
    %             indices = max(bsxfun(@times,abs(diff(isnan(floats.x),1,1)),[1:size(floats.x,1)-1]'+1));
    %             for ii = 1:size(floats.x,2)
    %                 clear ind
    %                 ind = find(floats.x(:,ii) > 0);
    %                 if isempty(ind)
    %                     ind = 1;
    %                 else
    %                     ind = ind(1);
    %                 end
    % 
    %                 floats.init(ii,:) = [floats.x(ind,ii) floats.y(ind,ii) floats.z(ind,ii) floats.time(ind)];
    %             end
            disp(['Finished processing ' upper(type) ' floats.']);
            toc;
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
                         [nansum(bsxfun(@minus,floats.x,floats.init(:,1)'),2) ...
                          nansum(bsxfun(@minus,floats.y,floats.init(:,2)'),2) ...
                          nansum(bsxfun(@minus,floats.z,floats.init(:,3)'),2)] ...
                          ,1./floats.N);
                      
             % cloud / relative dispersion
%             tic;  
%             s = zeros(T,3);
%             for i=1:size(floats.x,2)
%                 for j=1:size(floats.x,2)
%                     if i~=j
%                         s(:,1) = s(:,1) + repnan((floats.x(:,i) - floats.x(:,j)),0).^2/1e6;
%                         s(:,2) = s(:,2) + repnan((floats.y(:,i) - floats.y(:,j)),0).^2/1e6;
%                         s(:,3) = s(:,3) + repnan((floats.z(:,i) - floats.z(:,j)),0).^2;
%                     end
%                 end
%             end
%             toc;
%             tic;
%             nn = size(floats.x,2);
%             %[p,q] = meshgrid(1:nn,1:nn);
%             p = nchoosek(1:nn,2);
%             q = p(:,2); p = p(:,1);
%             s(:,1) = nansum((floats.x(:,p(:)) - floats.x(:,q(:))).^2,2)/1e6;
%             s(:,2) = nansum((floats.y(:,p(:)) - floats.y(:,q(:))).^2,2)/1e6;
%             s(:,3) = nansum((floats.z(:,p(:)) - floats.z(:,q(:))).^2,2);
%             toc;
            s = nan(T,3);
            tic;
            nn = size(floats.x,2);
            indices = [1:size(floats.x,2)];
            %[p,q] = meshgrid(1:nn,1:nn);
            p = nchoosek(1:nn,2);
            q = p(:,2); p = p(:,1);
            for kk = 1:size(floats.x,1)
                try
                    s(kk,1) = nansum((floats.x(kk,p(:)) - floats.x(kk,q(:))).^2,2)/1e6;
                    s(kk,2) = nansum((floats.y(kk,p(:)) - floats.y(kk,q(:))).^2,2)/1e6;
                    s(kk,3) = nansum((floats.z(kk,p(:)) - floats.z(kk,q(:))).^2,2);
                catch ME
                    s(kk,:) = NaN;
                end
            end
            toc;
            floats.disp = bsxfun(@times,s,1./(floats.N.*(floats.N-1)));
            
            % kurtosis - this could be improved?
            floats.kur = nan(T,3);
            floats.kur = [kurtosis(floats.x,0,2) kurtosis(floats.y,0,2) kurtosis(floats.z,0,2)];
%             dfx = bsxfun(@minus,floats.x,floats.mom1(:,1));
%             dfy = bsxfun(@minus,floats.y,floats.mom1(:,2));
%             dfz = bsxfun(@minus,floats.z,floats.mom1(:,3));
%             
%             floats.kur(:,1) = nansum(dfx.^4,2)./ ...
%                              (nansum(dfx.^2,2)).^2;
%             floats.kur(:,2) = nansum(dfy.^4,2)./ ...
%                              (nansum(dfy.^2,2)).^2;
%             floats.kur(:,3) = nansum(dfz.^4,2)./ ...
%                              (nansum(dfz.^2,2)).^2;
        end
        
        % plot stats in subplots
        function [hfig] = plot_stats(floats,hfig)
            if isempty(floats.kur), floats.stats; end
            if ~exist('hfig','var'), hfig = gcf; end
            if strcmp(floats.type,'ltrans')
                fmt = '--';
                %floats.time = floats.time + 4147200;
            else
                fmt = '-';
            end
            figure(hfig);
            subplot(311)
            plot(floats.time/86400,bsxfun(@times, floats.mom1,[ 1/1000 1/1000 10]),fmt); hold on;
            xlabel('time (days)'); ylabel('first moment');
            legend('x (km)','y (km)','10 * z (m)','Location','NorthWest');
            title(' dashed = LTRANS, line = ROMS');
            beautify([14 14 16]);
            subplot(312)
            plot(floats.time/86400,bsxfun(@times,floats.disp,[1 10 10]),fmt); hold on;
            xlabel('time (days)'); ylabel('Relative Dispersion');
            legend('x (km^2)','10*y (km^2)','10*z (m^2)','Location','NorthWest');
            beautify([14 14 16]);
            subplot(313)
            plot(floats.time/86400,floats.kur,fmt); hold on;
            xlabel('time (days)'); ylabel('Kurtosis');
            legend('x','y','z','Location','NorthWest');
            beautify([14 14 16]);
        end
        
        % animates with zeta plot
        function [] = animate(floats,rgrid,zeta,eddy)

            cmap = flipud(cbrewer('div', 'BrBG', 32));

            figure;
            tfilt = cut_nan(fillnan(floats.init(:,4),0));
            ind0 = find_approx(rgrid.ocean_time, tfilt(1),1);
            
            i = ind0;
            [~,hz] = contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,i)',25);
            shading flat; axis image
            colormap(cmap); 
            caxis([min(zeta(:)) max(zeta(:))]);
            hold on
            [C,hc] = contour(rgrid.x_rho./1000,rgrid.y_rho./1000,rgrid.h,[114 500 750 1100],'k');
            clabel(C,hc);
            floats.fac = floor(floats.fac);
            nn = find_approx(floats.time,rgrid.ocean_time(i),1);
            hplot = plot(floats.x(nn,:)/1000,floats.y(nn,:)/1000,'k.','MarkerSize',10);
            if exist('eddy','var')
               [~,hh] = contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,i),1);
               set(hh,'LineWidth',2);
            end
            ht = title(['t = ' num2str((rgrid.ocean_time(i)+1)/86400) ' days']);

            for i=ind0+1:size(zeta,3)
                set(hz,'ZData',zeta(:,:,i)'); shading flat
                nn = find_approx(floats.time,rgrid.ocean_time(i),1);
                set(hplot,'XData',floats.x(nn,:)/1000);
                set(hplot,'YData',floats.y(nn,:)/1000);
                if exist('eddy','var')
                   set(hh,'ZData',eddy.mask(:,:,i));
                end
                set(ht,'String',['t = ' num2str((rgrid.ocean_time(i)+1)/86400) ' days'])
                pause(0.03)
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