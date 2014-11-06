classdef floats < handle
    properties
        x; y; z; time; age;
        fac; type
        xsb; % shelfbreak location
        winds; % indices floats that crossed shelfbreak WEST of
               % eddy.
        einds; % indices of floats that crossed EAST of the eddy
        reenter; % number of floats that re-enter shelf after
                 % crossing shelfbreak
        reenter_inds; % indices of floats that crossed the
                      % shelfbreak to the WEST of the eddy and the
                      % returned to the shelf
        roms_inds; % indices of floats that started at same place,
                   % time as ROMS deployment
        actual_init; % for ROMS floats only. added if float output
                     % doesn't coincide with history file
                     % output. This is needed because stats are
                     % calculated based on first location in
                     % *flt.nc output.
        temp; salt; rho;
        hitLand; hitBottom
        init; N = 0;
        comment = ['init = (x,y,z,t) = initial location, release time in meters, seconds | ' ...
                   'mom1 = first moment | disp = cloud dispersion | kur = kurtosis | ' ...
                   'N = number of floats at a time instant.'];
        disp; mom1; kur;
    end
    methods
        % reads data
        function [floats] = floats(type, file, rgrid, xsb, fltname)
            disp('Reading float data.');
            tic;
            if strcmpi(type,'roms')
                try
                    floats.x = dc_roms_read_data(file,'lon',[],{}, ...
                                                 [],[],'flt',[])';
                catch ME
                    floats.x = dc_roms_read_data(file,'x',[],{}, ...
                                                 [],[],'flt',[])';
                end
                try
                    floats.y = dc_roms_read_data(file,'lat',[],{}, ...
                                                 [],[],'flt',[])';
                catch ME
                    floats.y = dc_roms_read_data(file,'y',[],{}, ...
                                                 [],[],'flt',[])';
                end
                floats.z = dc_roms_read_data(file,'depth',[],{}, ...
                                             [],[],'flt',[])';
                floats.time = dc_roms_read_data(file,'ocean_time',[],{}, ...
                                                [],[],'flt',[])';
                try
                    floats.temp = dc_roms_read_data(file,'temp',[],{}, ...
                                                 [],[],'flt',[])';
                    floats.salt = dc_roms_read_data(file,'salt',[],{}, ...
                                                 [],[],'flt',[])';
                catch
                    floats.rho = dc_roms_read_data(file,'rho',[],{}, ...
                                                 [],[],'flt',[])';
                    floats.salt = nan(size(floats.x));
                end
                floats.type = 'roms';
                disp('Read data. Now processing.');
                toc;tic;
            end

            if strcmpi(type,'tracpy')
                floats.x = ncread(file,'lonp');
                floats.y = ncread(file,'latp');
                try
                    floats.z = ncread(file,'zp');
                catch ME
                    floats.z = nan(size(floats.x));
                end
                floats.time = ncread(file,'tp');
                floats.temp = nan(size(floats.x));
                floats.salt = nan(size(floats.x));
                floats.rho = nan(size(floats.x));

                floats.type = 'tracpy';
                disp('Read data. Now processing.');
                toc;tic;
            end

            if strcmpi(type,'ltrans')
                floats.y = ncread(file,'lat')';
                floats.x = ncread(file,'lon')';
                floats.z = ncread(file,'depth')';
                floats.age = ncread(file,'age')';
                floats.time = ncread(file,'model_time');
                try
                    floats.temp = ncread(file,'temperature')';
                    floats.salt = ncread(file,'salinity')';
                catch ME
                    warning('No temperature or salt saved');
                end
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

                        floats.init(ii,:) = [floats.x(ind,ii) ...
                                            floats.y(ind,ii) ...
                                            floats.z(ind,ii) floats.time(ind)];
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
            if isvector(floats.time)
                tmat = repmat(floats.time,[1 size(floats.x,2)]);
            else
                tmat = floats.time;
                floats.time = floats.time(:,1);
            end
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
            floats.comment = ['init = (x,y,z,t) = initial location, ' ...
                              'release time in meters, seconds | ', ...
                              'fac = number float timesteps per ROMS output ' ...
                              'timestep'];

            % assign shelfbreak depth
            floats.xsb = xsb;
            % calculate fac
            try
                dtroms = rgrid.ocean_time(2)-rgrid.ocean_time(1);
                dtltr  = floats.time(2,1)-floats.time(1,1);
                floats.fac = dtroms/dtltr;
            catch ME
                 warning('havent calculated fac');
            end

            % initial locations and seeding time
            if strcmpi(type, 'tracpy')
                floats.init = [floats.x(1,:)' floats.y(1,:)' floats.z(1,:)' ...
                               tmat(1,:)'];
            end

            if strcmpi(type, 'roms')
                % The ROMS float output file does not store the
                % location when the float is deployed, so
                % parse the floats.in file to figure it out. This
                % requires fltname to be provided to the
                % constructor.

                floats.actual_init = floats.parse_roms_fltname(rgrid, fltname);
            end

            if strcmpi(type, 'ltrans') || strcmpi(type, 'roms')
                initmask = find([zeros([1 size(floats.x,2)]); ...
                                 abs(diff((cumsum(repnan(floats.x,0)) >= 1)))] == 1);
                floats.init = [floats.x(initmask) floats.y(initmask) ...
                               floats.z(initmask) tmat(initmask)];
            end

            disp(['Finished processing ' upper(type) ' floats.']);
            toc;
        end % read_floats

        function [] = plot_xz(floats, rgrid)
            hmin = min(rgrid.h(:));

            ind = floats.winds;

            figure;
            % draw bathy cross-section
            patch([rgrid.y_rho(:,1); min(rgrid.y_rho(:,1))]/1000, ...
                  -1*[hmin; rgrid.h(:,1)], 'k');
            limx = xlim;
            xlim([0 limx(2)]); ylim([-1*max(rgrid.h(:)) 0]);
            xlabel('Y (km)'); ylabel('Z (m)');

            % add float tracks
            hold on;

            latp = floats.y(:, ind)/1000;
            zp = floats.z(:, ind);

            plot(latp, zp, '-', 'Color', [1 1 1]*0.75);

            figure;
            ax(1) = subplot(2,1,1);
            hist(zp(end,:));
            beautify;
            title('Final z-locations of particles');
            ax(2) = subplot(2,1,2);
            hist(min(zp, [], 1));
            beautify;
            title(['Minimum z-location that a particle has ' ...
                   'visited']);
            linkaxes(ax, 'xy');
        end

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
        function [hfigs] = plot_stats(floats, hfigs)
            if isempty(floats.kur), floats.stats; end
            if ~exist('hfigs','var')
                hfigs(1) = figure;
                hfigs(2) = figure;
            end
            if strcmp(floats.type,'ltrans') || strcmp(floats.type, 'tracpy')
                fmt = '--';
                %floats.time = floats.time + 4147200;
            else
                fmt = '-';
            end

            ind = [];
            if isempty(ind)
                ind = 1:size(floats.mom1, 1);
            end

            figure(hfigs(1));
            subplot(311)
            plot(floats.time/86400,bsxfun(@times, floats.mom1(ind,1:2), ...
                                          [1/1000 1/1000]),fmt);
            hold on;
            xlabel('time (days)'); ylabel('first moment');
            legend('x (km)','y (km)','Location','NorthWest');
            title([' dashed = ' floats.type ', line = ROMS']);
            beautify([14 14 16]);

            subplot(312)
            plot(floats.time/86400, ...
                 bsxfun(@times,floats.disp(ind,1:2), [1 1]),fmt); hold on;
            xlabel('time (days)'); ylabel('Relative Dispersion');
            legend('x (km^2)','y (km^2)','Location','NorthWest');
            beautify([14 14 16]);

            subplot(313)
            plot(floats.time/86400, floats.kur(ind,1:2),fmt); hold on;
            xlabel('time (days)'); ylabel('Kurtosis');
            legend('x','y','Location','NorthWest');
            beautify([14 14 16]);

            figure(hfigs(2));
            subplot(311)
            plot(floats.time/86400,bsxfun(@times, floats.mom1(ind,3), ...
                                          [1]),fmt);
            hold on;
            xlabel('time (days)'); ylabel('first moment');
            legend('x (km)','y (km)','Location','NorthWest');
            title([' dashed = ' floats.type ', line = ROMS']);
            beautify([14 14 16]);

            subplot(312)
            plot(floats.time/86400, ...
                 bsxfun(@times,floats.disp(ind,3), [1]),fmt); hold on;
            xlabel('time (days)'); ylabel('Relative Dispersion');
            legend('x (km^2)','y (km^2)','Location','NorthWest');
            beautify([14 14 16]);

            subplot(313)
            plot(floats.time/86400, floats.kur(ind,3),fmt); hold on;
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
            [C,hc] = contour(rgrid.x_rho./1000,rgrid.y_rho./1000, ...
                             rgrid.h,[114 500 750 1100],'k');
            clabel(C,hc);
            floats.fac = floor(floats.fac);
            nn = find_approx(floats.time,rgrid.ocean_time(i),1);
            hplot = plot(floats.x(nn,:)/1000,floats.y(nn,:)/1000,'k.','MarkerSize',10);
            if exist('eddy','var')
               [~,hh] = contour(eddy.xr/1000,eddy.yr/1000,eddy.vormask(:,:,i),1);
               set(hh,'LineWidth',2);
            end
            ht = title(['t = ' num2str((rgrid.ocean_time(i)+1)/86400) ' days']);

            for i=ind0+1:size(zeta,3)
                set(hz,'ZData', double(zeta(:,:,i)')); shading flat
                nn = find_approx(floats.time,rgrid.ocean_time(i),1);
                set(hplot,'XData',floats.x(nn,:)/1000);
                set(hplot,'YData',floats.y(nn,:)/1000);
                if exist('eddy','var')
                   set(hh,'ZData',eddy.vormask(:,:,i));
                end
                set(ht,'String',['t = ' num2str((rgrid.ocean_time(i)+1)/86400) ' days'])
                pause(0.03)
            end
        end

        function [init] = parse_roms_fltname(floats, rgrid, fpos_file)

            % find line just before float location specification
            [~,p] = grep('POS = ', fpos_file);

            % make sure there's only one match above
            assert(length(p.line) == 1);

            % Grid units  = (i,j)
            % Initial floats locations for all grids:
            %
            %   1 G      Nested grid number
            %   2 C      Initial horizontal coordinate type (0: grid units, 1: spherical)
            %   3 T      Float trajectory type (1: Lagrangian, 2: isobaric, 3: Geopotential)
            %   4 N       Number floats to be released at (Fx0,Fy0,Fz0)
            %   5 Ft0    Float release time (days) after model initialization
            %   6 Fx0    Initial float X-location (grid units or longitude)
            %   7 Fy0    Initial float Y-location (grid units or latitude)
            %   8 Fz0    Initial float Z-location (grid units or depth)
            %   9 Fdt    Float cluster release time interval (days)
            %  10 Fdx    Float cluster X-distribution parameter
            %  11 Fdy    Float cluster Y-distribution parameter
            %  12 Fdz    Float cluster Z-distribution parameter

            fid = fopen(fpos_file);
            readstr = textscan(fid, '%d%d%d%d%f%f%f%f%f%f%f%f', ...
                               'HeaderLines', p.line+1, 'CommentStyle', ...
                               '!');
            fclose(fid);

            % go through each line
            for ii=1:length(readstr{1})
                if readstr{2}(ii) == 0

                    N = readstr{4}(ii);
                    t0 = readstr{5}(ii);
                    x0 = readstr{6}(ii);
                    y0 = readstr{7}(ii);
                    z0 = readstr{8}(ii) + 1; % +1 since ROMS has 0
                                             % as the bottom level
                                             % but MATLAB uses 1
                    dt = readstr{9}(ii);
                    dx = readstr{10}(ii);
                    dy = readstr{11}(ii);
                    dz = readstr{12}(ii);

                    % figure out how many dimensions to distribute
                    % across. ASSUMING DC_FLOATS_DEPLOYMENT IS ENABLED
                    count = 0;
                    if dt > 0, count = count + 1; end
                    if dx > 0, count = count + 1; end
                    if dy > 0, count = count + 1; end
                    if dz > 0, count = count + 1; end

                    nlimit = floor(double(N)^(1/count));
                    % this is the actual number of floats deployed
                    N = nlimit ^ (count);

                    tlimit = 0; xlimit = 0; ylimit = 0; zlimit = 0;

                    % limits for loops
                    if dt > 0, tlimit = nlimit; end
                    if dx > 0, xlimit = nlimit; end
                    if dy > 0, ylimit = nlimit; end
                    if dz > 0, zlimit = nlimit; end

                    % lengths of grid vectors
                    X = size(rgrid.x_rho, 2);
                    Y = size(rgrid.x_rho, 1);
                    Z = size(rgrid.z_r, 1);

                    k=1;
                    init = nan(N,4);
                    for xx=1:xlimit
                        for yy=1:ylimit
                            for zz=1:zlimit
                                for tt=1:tlimit
                                    init(k,1) = interp1(1:X, ...
                                                               rgrid.x_rho(1,:), ...
                                                               x0 + (xx-1)*dx);
                                    init(k,2) = interp1(1:Y, ...
                                                               rgrid.y_rho(:,1), ...
                                                               y0 + (yy-1)*dy);
                                    %init(k,3) = interp1(1:Z, ...
                                    %                           rgrid.z_r(:, ...
                                    %                                     y0+(yy-1)*dy, ...
                                    %                                     x0+(xx-1)*dx), ...
                                    %                           z0 + (zz-1)*dz);
                                    init(k,3) = z0+(zz-1)*dz;
                                    init(k,4) = (t0 + (tt-1)*dt)*86400;
                                    k=k+1;
                                end
                            end
                        end
                    end
                else
                    warning(['Not ready to interpret non-grid-co-ordinate ' ...
                             'floats.in file']);
                end
            end
        end
    end
end