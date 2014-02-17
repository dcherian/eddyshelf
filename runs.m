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
    csdye; asdye; zdye; eddye; % cross-shore, along-shore, z dyes, eddy dye
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
    % streamer properties
    streamer;
    % rossby radii
    rrdeep; rrshelf
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
        
        % get grid
        zeta0 = double(ncread(runs.out_file,'zeta',[1 1 1],[Inf Inf 1]));
        runs.rgrid = roms_get_grid(runs.out_file,runs.out_file, ...
                        zeta0',1);
        runs.rgrid.xr = runs.rgrid.x_rho';
        runs.rgrid.yr = runs.rgrid.y_rho';
        runs.rgrid.zr = permute(runs.rgrid.z_r,[3 2 1]);
        runs.rgrid.z_uw = [];
        runs.rgrid.z_vw = [];
        runs.rgrid.zeta = [];

        % read zeta
        if ~runs.givenFile
            runs.zeta = dc_roms_read_data(dir,'zeta',[],{},[],runs.rgrid);
            runs.time = roms_read_data(dir,'ocean_time');%,[],{},[],runs.rgrid);
            %try
            %    runs.csdye  = roms_read_data(dir,'dye_01', ...
            %        [1 1 runs.rgrid.N 1],[Inf Inf 1 Inf]);
            %catch ME
            %end
        else
            runs.zeta = double(ncread(runs.out_file,'zeta'));
            runs.time = double(ncread(runs.out_file,'ocean_time'));
        end


        warning('Assuming uniform grid for dx,dy');
        runs.rgrid.dx = mean(diff(runs.rgrid.xr(:,1),1,1));
        runs.rgrid.dy = mean(diff(runs.rgrid.yr(1,:),1,2));
        runs.rgrid.dV = runs.rgrid.dx * runs.rgrid.dy * ...
            diff(permute(runs.rgrid.z_w,[ 3 2 1]),1,3); 

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
        [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = ...
                        find_shelfbreak(runs.out_file);
        [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = ...
                        find_shelfbreak(runs.out_file,'slope');
        runs.bathy.h = runs.rgrid.h';

        % rossby radii
        runs.rrdeep = sqrt(runs.params.phys.N2)*max(runs.bathy.h(:)) ...
                    /mean(runs.rgrid.f(:));
        runs.rrshelf = sqrt(runs.params.phys.N2)*max(runs.bathy.hsb) ...
                    /mean(runs.rgrid.f(:));

%             % read in dye surface fields
%             for ii=1:4
%                 % dye name
%                 dname = ['dye_0' num2str(ii)];
%                 try % see if variable exists in ini
%                     vname = [];
%                     % dye description
%                     ddesc = ncreadatt([runs.dir roms_find_file(runs.dir,'ini')], ...
%                                         dname,'long_name');
%                     if strfind(ddesc,'cross shelf'), vname = 'csdye'; end
%                     if strfind(ddesc,'z dye'), vname = 'zdye'; end
%                     if strfind(ddesc,'along shelf'), vname = 'asdye'; end
%                     if strfind(ddesc,'eddy dye'), vname = 'eddye'; end
%                     
%                     % see if variable is in output files
%                     try
%                         runs.(vname) = roms_read_data(filename,dname ...
%                            ,[1 1 runs.rgrid.N 1],[Inf Inf 1 Inf]);
%                     catch ME
%                         warning([dname 'not in output files']);
%                     end
%                 catch ME
%                     warning([dname 'not found in ini file']);
%                 end
%             end

        try
            runs.roms = floats('roms',runs.flt_file,runs.rgrid);
        catch
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

        if ~isfield(runs.eddy,'vor')
            runs.eddy = track_eddy(dir);
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
            try
                runs.eddy.trevind = find(runs.eddy.cvx < 0,1,'first');
                runs.eddy.trev = runs.time(runs.eddy.trevind);
            catch ME
                runs.eddy.trev = nan;
            end
            if isempty(runs.eddy.trev), runs.eddy.trev = NaN; end

            % estimate southward vel.
            % (beta * Lr^2)^2 *1/2 * 1/amp * NH/g
            runs.eddy.Vy = -(runs.params.phys.beta*(runs.params.eddy.dia(1)/2)^2)^2/2 ...
                                    /runs.eddy.amp(1) * ...
                            sqrt(runs.params.phys.N2)/runs.params.phys.g*max(runs.bathy.h(:));

%                            -(run3.params.phys.beta*(run3.params.eddy.dia(1)/2)^2)^2 ...
%                                        *1/2 * 1/run3.eddy.amp(1) * ...
%                                sqrt(run3.params.phys.N2)/run3.params.phys.g*max(run3.bathy.h(:))

          %  % water depth at eddy center
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

    %% conservation checks
    function [] = check_temp(runs)

        visc2 = ncread(runs.out_file,'visc2_r');
        visc2 = visc2 - min(visc2(:));

        figure;
        subplot(121)
        pcolorcen(runs.zeta(:,:,1)');
        hold on
        contour(visc2',[1 1]*3,'k');
        caxis([min(runs.zeta(:)) max(runs.zeta(:))]);

        n = 15;
        [x,y] = ginput(n);
        xi = ceil(x); yi = ceil(y);
        plot(xi,yi,'x','MarkerSize',12);

        for ii=1:n
            text(xi(ii),yi(ii),num2str(ii));
            temp(:,:,ii) = dc_roms_read_data(runs.dir,'temp', [], ...
                    {'x' xi(ii) xi(ii); 'y' yi(ii) yi(ii)});

            dz(:,ii) = diff(runs.rgrid.z_w(:,yi(ii),xi(ii)));
        end

        % depth integrated
        itemp = squeeze(sum( ...
                    bsxfun(@times, temp, permute(dz,[1 3 2])), 1));

        % depth averaged
        atemp = bsxfun(@rdivide,itemp,diag(runs.rgrid.h(yi,xi))');
        subplot(122)
        plot(bsxfun(@minus,atemp, mean(atemp,1)));
        legend(gca,'show');
        xlabel('Time (days)'); ylabel('Depth averaged temperature (without mean)');
    end

   %% analysis

    function [] = transport(runs)
        % need some kind of initial time instant - probably objective
        % criterion based on distance between shelfbreak and eddy or
        % some sort
        runs.eutrans = [];
        t0 = 1;
        revind = runs.eddy.trevind;
        h = runs.bathy.h(2:end-1,2:end-1);

        ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx(t0:end));
        iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.my(t0:end));
        hcen = h(sub2ind(size(runs.eddy.xr),ix,iy))';

%redundant            ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx(t0:end));
        iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.se(t0:end));
        hedge = h(sub2ind(size(runs.eddy.xr),ix,iy))';
        % rossby radius @ SHELFBREAK
        rr =  sqrt(runs.params.phys.N2)*runs.bathy.hsb ...
                    /runs.rgrid.f(runs.bathy.isb,1)
        distance = 5*runs.rr; % 5 times rossby radius

        if runs.params.bathy.axis == 'x'
            csvelid = 'u';
            error(' not built for north-south isobaths');
        else
            csvelid = 'v';
            loc = sort([nanmean(runs.eddy.se(revind:end)) ...
                    nanmean(runs.eddy.cy(revind:end)) ...
                    runs.bathy.xsb  ... 
                    runs.bathy.xsl]); 
                %runs.rgrid.y_rho(vecfind(runs.bathy.h(1,:),[250 1000]),1)']);
        end

        % save locations
        runs.eutrans.x = loc;
        % save indices for locations
        runs.eutrans.ix = vecfind(runs.rgrid.yr(1,:),loc);%find_approx(runs.rgrid.yr(1,:),loc,1);
        % save isobath values
        runs.eutrans.h = ceil(runs.bathy.h(1,runs.eutrans.ix));
        % find west edge indices
        iwest = vecfind(runs.eddy.xr(:,1),runs.eddy.we);
        dx = runs.rgrid.xr(2,1)-runs.rgrid.xr(1,1);

        % loop over all isobaths
        for kk=1:length(loc)
            % read along-shore section of cross-shore vel.
            % dimensions = (x/y , z , t )
            cs_vel = double(squeeze(ncread(runs.out_file,csvelid, ...
                [1 runs.eutrans.ix(kk) 1 t0],[Inf 1 Inf Inf])));
            mask = nan(size(cs_vel)); 
            for tt=1:size(cs_vel,3)
                mask(1:iwest(tt),:,tt) = 1;
            end
            % restrict calculation to region above shelfbreak depth
            zmask = (abs(squeeze(runs.rgrid.z_r(:,runs.eutrans.ix(kk),:))   )' ...
                            < runs.bathy.hsb);
            mask = bsxfun(@times,mask,fillnan(zmask,0));

%                 runs.eutrans.trans(:,:,kk) = squeeze( ...
%                      trapz(runs.rgrid.z_r(:,runs.eutrans.ix(kk),1),mask .* cs_vel,2));
            %runs.eutrans.trans(:,:,kk) = squeeze( ...
            %            nansum(bsxfun(@times,avg1(mask .* cs_vel,2), ...
             %           diff(runs.rgrid.z_r(:,runs.eutrans.ix(kk),1))'),2));
            runs.eutrans.nodye.trans(:,:,kk) = squeeze(trapz( ...
                            runs.rgrid.z_r(:,runs.eutrans.ix(kk),1), ...
                            mask .* cs_vel,2));

            runs.eutrans.nodye.Itrans(:,kk) = squeeze(nansum( ...
                runs.eutrans.nodye.trans(:,:,kk) ...
                                        .* dx,1))';

            % if I have passive tracer info I can calculate transport
            % using that
            if ~isempty(runs.csdye)
                mask = nan(size(cs_vel));
                % mark eastern edge as edge of region I'm interested in
                % removes streamer associated with cyclone running away
                ieast = vecfind(runs.eddy.xr(:,1),runs.eddy.ee);
                for tt=1:size(cs_vel,3)
                    mask(1:ieast(tt),:,tt) = 1;
                end

                mask = bsxfun(@times,mask,fillnan(zmask,0));
                % dye_01 is always cross-shore dye
                dye = double(squeeze(ncread(runs.out_file,'dye_01', ...
                    [1 runs.eutrans.ix(kk) 1 t0],[Inf 1 Inf Inf])));
                dyemask = (dye >= runs.bathy.xsb) & ...
                            (dye <=(runs.bathy.xsb + distance));
                mask = mask .* fillnan(dyemask,0); 
                runs.eutrans.trans(:,:,kk) = squeeze(trapz( ...
                        runs.rgrid.z_r(:,runs.eutrans.ix(kk),1), ...
                        mask .* cs_vel,2));
                %runs.eutrans.trans(:,:,kk) = squeeze( ...
                %        nansum(bsxfun(@times,avg1(mask .* cs_vel,2), ...
                %        diff(runs.rgrid.z_r(:,runs.eutrans.ix(kk),1))'),2));

                runs.eutrans.Itrans(:,kk) = squeeze(nansum( ...
                            runs.eutrans.trans(:,:,kk) .* dx,1))';
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
        plot(runs.rgrid.ocean_time/86400,runs.eutrans.Itrans/1e6,'-');
        limx = xlim;
        legend(num2str(runs.eutrans.h'),'Location','NorthWest');
        ylabel('Dye Transport (Sv)');
        ylim([-0.05 0.3]); liney(0.1,[])
        beautify;
        subplot(6,1,6)
        [ax,~,~] = plotyy(runs.eddy.t,runs.eddy.prox/1000,runs.eddy.t, ...
                runs.eddy.hcen);
        set(ax(1),'XLim',limx);set(ax(2),'XLim',limx);
        set(ax(1),'XTickLabel',[]); axes(ax(2));
        set(get(ax(1),'ylabel'),'String','Proximity (km)');
        set(get(ax(2),'ylabel'),'String','h @ center of eddy');
        xlabel('Time (days)');

        % throw out locations where dye trans is pretty much zero to
        % make plot cleaner
        arr = [1:length(loc)];
        for kk=1:length(loc)
            if median(runs.eutrans.Itrans(:,kk)) < 1
                arr(arr == kk) = [];
            end
        end
        figure
        plot(runs.rgrid.ocean_time(t0:end)/86400,(runs.eutrans.nodye.Itrans(:,arr) - runs.eutrans.Itrans(:,arr)) ...
            ./ runs.eutrans.Itrans(:,arr) * 100);
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

    function [] = plot_simplepv(runs)
       % this function contours the qgpv approximation of the
       % background pv

       if runs.bathy.axis == 'y'
           dhdx = diff(runs.bathy.h,1,2)./diff(runs.rgrid.yr,1,2);
           ax = 2;
       else
           ax = 1;
           dhdx = diff(runs.bathy.h,1,1)./diff(runs.rgrid.xr,1,1);
       end

       beta_t = runs.params.phys.f0 * dhdx/max(runs.rgrid.zr(:));

       q = runs.params.phys.f0 + ...
           (runs.params.phys.beta + beta_t) .* avg1(runs.rgrid.yr,ax);

       clf;
       subplot(211);
       contourf(q');
       subplot(212);
       hold on
       plot(q(2,:));
       plot(-runs.bathy.h(2,:)/max(runs.bathy.h(:)),'k');
       legend('qgpv','bathy');

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
        tloc = [100 120];
        limx = [0 max(runs.time)/86400];
        hold on
        subplot(aa,2,1); hold all
        %pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);            colorbar
        xlabel('X (km)'); ylabel('Y (km)');
        plot((eddy.cx-eddy.cx(1))/1000, ...
            (eddy.cy-eddy.cy(1))/1000,'Color',colors(ii,:),'LineWidth',2);
        hold all;
        try
            plot((eddy.cx(tloc)-eddy.cx(1))/1000,...
                (eddy.cy(tloc)-eddy.cy(1))/1000,'*','MarkerSize',12,'Color',colors(ii,:),'LineWidth',2);
        catch ME
        end
            %if runs.bathy.axis == 'x'
        %    plot(eddy.we/1000,eddy.cy/1000,'Color',colors(ii,:),'LineStyle','--');
        %else
        %    plot(eddy.cx/1000,eddy.se/1000,'Color',colors(ii,:),'LineStyle','--');
        %end
        %axis image; axis tight

        subplot(aa,2,2); hold on
        plot(eddy.t,eddy.vor.amp./eddy.amp(1),'Color',colors(ii,:));
        ylabel('amp/amp(1) ');xlim(limx);

        subplot(aa,2,4); hold on
        plot(eddy.t,eddy.vor.dia/1000,'Color',colors(ii,:));
        ylabel('diameter (km)');xlim(limx);

        subplot(aa,2,6); hold on
        plot(eddy.t,eddy.cvx,'Color',colors(ii,:));
        ylabel('cvx(km/day)');
        ylim([-5 5]);
        liney(0); xlim(limx);
        %plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
        %ylabel('x - center (km)');

        subplot(aa,2,8); hold on
        plot(eddy.t,eddy.cvy,'Color',colors(ii,:));
        ylabel('cvy (km/day)');xlim(limx);
        ylim([-5 5]);
        %plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
        %ylabel('y - center (km)');

        subplot(aa,2,10); hold on
        plot(eddy.t,eddy.Lz2,'Color',colors(ii,:));
        ylabel('vertical scale (m)');xlim(limx);
        %xlabel('time (days)');

        subplot(aa,2,12); hold on
        plot(eddy.t,eddy.prox/1000,'Color',colors(ii,:));
        xlabel('time (days)');
        ylabel('Proximity (km)');xlim(limx);

        subplot(aa,2,11); hold on
        hp = plot(eddy.t,eddy.hcen,'Color',colors(ii,:));
        addlegend(hp,runs.name,'SouthWest');
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

    % detect streamer contour and figure out cross-section points
    function [] = detect_streamer_mask(runs)

        % upper y-limit to save memory
        yend = find_approx(runs.rgrid.y_rho(:,1),130*1000);
        t0 = runs.eddy.trevind;
        %read_start = [1 1 1 t0-20];
        slab = 40;
        
        runs.streamer.yend = yend;
        
        szeta = size(runs.zeta);
        
        % allocate memory
        nanvec = nan(size(runs.time));
        runs.streamer.time = nanvec;
        runs.streamer.west.vol = nanvec;
        runs.streamer.west.zcen = nanvec;
        runs.streamer.west.zdcen = nanvec;
        
        % grid matrices required for plotting     
        xsb = runs.bathy.xsb/1000;
        xr = runs.rgrid.xr(:,1:yend)/1000; 
        yr = runs.rgrid.yr(:,1:yend)/1000;
        zr = permute(runs.rgrid.z_r(:,1:yend,:),[3 2 1]);
        ix = repmat([1:size(xr,1)]',[1 yend]);
        iy = repmat([1:yend],[size(xr,1) 1]);  
       
        runs.streamer.xr = xr;
        runs.streamer.yr = yr;
        runs.streamer.zr = zr;
        
        % size matrices to make processing easier
        runs.streamer.sz4dfull = [size(zr) szeta(3)];
        runs.streamer.sz4dsp = [numel(zr) szeta(3)];
        runs.streamer.sz3dsp = [numel(zr) 1];
        runs.streamer.sz3dfull = size(zr);
        
        % allocate streamer mask variable
        runs.streamer.west.mask = sparse(runs.streamer.sz4dsp(1),szeta(3));
        
        % grid cell volume
        dVs = reshape(runs.rgrid.dV(:,1:runs.streamer.yend,:), ...
                        runs.streamer.sz3dsp); 

        for ii=1:floor(szeta(3)/slab)
            
            tstart = t0+slab*(ii-1);
            tend = tstart+slab-1;
            if tend > szeta(3), tend = szeta(3); end
            
            sz4dfull = [runs.streamer.sz4dfull(1:3) tend-tstart+1];
            sz4dsp = [runs.streamer.sz4dsp(1) tend-tstart+1];
            sz3dsp = runs.streamer.sz3dsp;
            sz3dfull = runs.streamer.sz3dfull;
        
            
            runs.streamer.time(tstart:tend) = runs.time(tstart:tend);
            tindices = [tstart tend];
            
            csdye = dc_roms_read_data(runs.dir, 'dye_01', tindices, ...
                        {'y' 1 yend},[],runs.rgrid)/1000;
            zdye  = dc_roms_read_data(runs.dir, 'dye_02', tindices, ...
                        {'y' 1 yend},[],runs.rgrid);
            %asdye = dc_roms_read_data(runs.dir, 'dye_03', tindices, ...
            %            {'y' 1 yend});

            % identify streamer with 4D data
            % preliminary detection
            % I use cross-shore label to roughly filter first
            % eliminate this step somehow and remove the temporary array?
            streamer1 = (csdye > xsb-10) & (csdye < xsb+30);

            % (xs,ys,zs) are the Eulerian x,y,z values
            %xs = bsxfun(@times, streamer, grd.xax)/1000;
            %ys = bsxfun(@times, streamer1, yr);
            %zs = bsxfun(@times, streamer, grd.zax);

            % (as,cs,z) dyes contain the Lagrangian labels
            % some distance metric between the two will give me an idea of
            % what's happening
            %if runs.bathy.axis == 'y'
            %    das = asdye - xs;
            %    dcs = csdye - ys;
            %else
            %    das = asdye - ys;
            %    dcs = csdye - xs;
            %end
            %dz = zdye - zs;

            % NEED TO ACCOUNT FOR TILTING IN VERTICAL?
            %cx = runs.eddy.cx(tstart:tend)/1000;
            %cy = runs.eddy.cy(tstart:tend)/1000;
            ee = runs.eddy.ee(tstart:tend)/1000;
            % hack if eddy center is outside extracted domain
            %cy(cy > max(yr(:))) = max(yr(:));
            %cxind = vecfind(xr(:,1),cx);
            %cyind = vecfind(yr(1,:),cy)';

            % pick out western streamer by chucking points that are >
            % eastern edge + initial radius. This allows streamer to wrap
            % around eddy and not be chopped off.
            west_mask = bsxfun(@lt, repmat(xr,[1 1 runs.rgrid.N]), ...
                          permute(ee + runs.params.eddy.dia/2000,[3 4 1 2]));
                      
            % mask of points west of eddy center - OLD ATTEMPT
            %west_mask = bsxfun(@lt,repmat(runs.rgrid.x_rho',[1 1 runs.rgrid.N]), ...
            %               permute(runs.eddy.cx(runs.eddy.trevind:runs.eddy.trevind+19), [1 3 4 2]));


            %r = sqrt(bsxfun(@minus,xr,permute(cx,[3 1 2])).^2 ...
            %       + bsxfun(@minus,yr,permute(cy,[3 1 2])).^2);

            % picking only western streamer
            streamer1 = squeeze(streamer1  ... % original streamer
                        ... % parcels have moved more than 5 km 
                        ... %   in the cross-shelf dirn.
                             .* (abs(csdye - bsxfun(@times,streamer1,yr))>5)) ...
                        ... % remove eastern half
                             .* (west_mask);
                        %...     % streamer depth is not total depth
                        %.* squeeze(bsxfun(@lt,max(abs(zs),[],3), runs.rgrid.h(1:yend,:)'));

            % pick out biggest surface piece
            % it looks like the surface is the biggest so we look for
            % only look for everything under it - i.e., hopefully no
            % tilting
            
            stream = streamer1(:,:,40,:);
            for tt=1:size(stream,4)
                % get biggest part - assume it's what i'm interested in
                strtemp = stream(:,:,1,tt);
                strcomps = bwconncomp(strtemp);
                numPixels = cellfun(@numel,strcomps.PixelIdxList);
                [~,bigidx] = max(numPixels);
                strtemp(strcomps.PixelIdxList{bigidx}) = 2;
                strtemp(strtemp < 2) = 0; 
                strtemp(strtemp == 2) = 1;
                stream(:,:,1,tt) = strtemp;
            end
            
            % filter and save
            runs.streamer.west.mask(:,tstart:tend) = sparse(reshape( ...
                            bsxfun(@times,streamer1,stream), ...
                            sz4dsp));
            %clear west_mask streamer1 stream strtemp;
            
            % compress somehow
            %streamnan = fillnan(runs.streamer.west.mask,0);
            % calculate statistics
            %xs = bsxfun(@times, streamnan, xr);
            %ys = bsxfun(@times, streamnan, yr);
            zs = bsxfun(@times, runs.streamer.west.mask(:,tstart:tend), ...
                            reshape(zr,sz3dsp));

            zdyestr = runs.streamer.west.mask(:,tstart:tend) .* ...
                            reshape(zdye,sz4dsp);
            %csdyestr = bsxfun(@times, streamnan, csdye);

            %dcs  = abs(csdyestr - ys);
            %dzd  = abs(zdyestr - zs);

            % calculate volume
            runs.streamer.west.vol(tstart:tend) = runs.domain_integratesp( ...
                runs.streamer.west.mask(:,tstart:tend), dVs);

            % Haven't used temperature yet

            % suffix cen = just centroids
            % suffix dcen = centroid weighted by dye value
            runs.streamer.west.zcen(tstart:tend) = bsxfun(@rdivide, ...
                runs.domain_integratesp(zs,dVs), ...
                runs.streamer.west.vol(tstart:tend));
            runs.streamer.west.zdcen(tstart:tend) = bsxfun(@rdivide,...
                runs.domain_integratesp(zdyestr,dVs), ...
                runs.streamer.west.vol(tstart:tend));
            
            % volume v/s depth plot for streamer
            % VECTORIZE SOMEHOW
            disp('Binning streamer volume...');
            tic;
            dbin = 20;
            bins = -1*0:dbin:1000;
            % required so that 0 bin doesn't get a ton of points
            %zsf = fillnan(full(zs),0);
            %sz = size(runs.streamer.west.mask(:,tstart:tend));
            for kk=1:length(bins)-1
               runs.streamer.west.Vbin(kk,tstart:tend) =  ...
                           sum(bsxfun(@times, (zs < bins(kk) & zs >= bins(kk+1)), ...
                               dVs),1);
            end
            runs.streamer.bins = bins;
            toc;
        end
    end
    
    % extract points for streamer section
    function [] = build_streamer_section(runs)
        
        % make plots to check?
        debug_plot = 1;
        
        runs.streamer.west.fit_circle = 1;
        
        if ~isfield(runs.streamer,'yend')
            runs.detect_streamer_mask();
        end
        yend = runs.streamer.yend;
        xr = runs.rgrid.xr(:,1:yend)/1000; 
        yr = runs.rgrid.yr(:,1:yend)/1000;
        
        cx = runs.eddy.mx/1000;
        cy = runs.eddy.my/1000;
        cy(cy > max(yr(:))) = max(yr(:));
        
        cxind = vecfind(xr(:,1),runs.eddy.mx/1000);
        cyind = vecfind(yr(1,:),cy)';
        
        for tind=1:size(runs.streamer.west.mask,2)
            % now pick ONLY SURFACE
            stream = reshape(full(runs.streamer.west.mask(:,tind)), ...
                runs.streamer.sz3dfull);
            stream = stream(:,:,end); % SURFACE ONLY
            
            % if no streamer or too small, skip
            if isequal(stream,zeros(size(stream))) ...
                    || numel(find(stream(:) ~= 0)) < 150
                continue;
            end
            
            % code from
            % http://blogs.mathworks.com/steve/2014/01/07/automating-data-extraction-2/x
            skeleton = bwmorph(stream,'skel','inf');
            skel = breakapart(skeleton);
            skelcomps = bwconncomp(skel);
            % find distance from eddy center?
            distcen = sqrt( (skel.*xr - runs.eddy.cx(tind)).^2 +  ...
                            (skel.*yr - runs.eddy.cy(tind)).^2 );
            distcen = distcen .* fillnan(skel,0);
            meandist = nan([skelcomps.NumObjects 1]);
            
            icen = nan(skelcomps.NumObjects,2);
            
            % process the branches for mean distance, centroid, and sort
            % clockwise
            for mm = 1:skelcomps.NumObjects
                meandist(mm) = nanmean(distcen(skelcomps.PixelIdxList{mm}));
                
                [ixtemp,iytemp] = ind2sub(size(skel), ...
                    skelcomps.PixelIdxList{mm});
                % don't remap to preserve order of points crossing the
                % horizontal axis
                %[~,sorttang] = angleSort([ixtemp iytemp], ...
                %                [cxind(tind) cyind(tind)],-pi/2);
                %sorttang = flipdim(sorttang,1);
                % works with 0 crossing
                tempang = atan2(iytemp-cyind(tind),ixtemp-cxind(tind));
                
                if max(diff(tempang) > 5.9)
                    tempang = mod(tempang + 2*pi,2*pi);
                end
                %tempang(tempang < 0) = tempang(tempang < 0) + 360;
                %[~,sorttang] = sort(tempang,'descend');
                %skelcomps.PixelIdxList{mm} = skelcomps.PixelIdxList{mm}(sorttang);
                %if ~isclockwise(ixtemp,iytemp)
                                
                % sort points in each branch clockwise. This is imposed by
                % setting the reference angle (w.r.t eddy center) to be the
                % minimum of all point angles in the branch
                refAngle = min(tempang(:));
                [out,~] = angleSort([ixtemp iytemp], ...
                    [cxind(tind) cyind(tind)],refAngle);
                skelcomps.PixelIdxList{mm} = sub2ind(size(skel), ...
                    flipud(out(:,1)),flipud(out(:,2)));
                %testBranch(skelcomps.PixelIdxList{mm},size(skel));
                
                % store centroid and find it's angle w.r.t eddy center
                icen(mm,:) = centroid([ixtemp(:) iytemp(:)]);
            end
            
            % sort by distance
            [~,sortdist] = sort(meandist);
            %, then chuck top 20%
            %indices = cat(1, ...
            %    skelcomps.PixelIdxList{ sortdist(1: floor(0.8*length(sortdist)) ) });
            %indices = skelcomps.PixelIdxList{sortdist(1)};
            
            % measure number of pixels in each branch and 
            % throw out small branches
            numPixels = cellfun(@numel,skelcomps.PixelIdxList);
            numPixels(numPixels < 5) = NaN;
            [~,sortnum] = sort(numPixels);
            nanindices = cut_nan(fillnan(isnan(numPixels) ...
              .* (1:skelcomps.NumObjects),0));
            for mm=1:length(nanindices)
              sortnum(sortnum == nanindices(mm)) = NaN;
            end
            
            % remove farthest away segment for sure
            if skelcomps.NumObjects > 1
                sortnum( sortnum == sortdist(end) ) = NaN;
            end
            % chuck out indices I'm not interested in
            sortnum = cut_nan(sortnum);
            
            % if region is too small, exit
            if isempty(sortnum)
                warning(['skipping @ tt=' num2str(tind)]);
                continue;
            end
            
            if runs.streamer.west.fit_circle
                % first get discrete points
                % old version without joining
                indices = cat(1,skelcomps.PixelIdxList{sortnum});
                [ixstr,iystr] = ind2sub(size(skel),indices);
                xstr = xr(ixstr,1);
                ystr = yr(1,iystr)';
                
                % fit circle
                circ = CircleFitByPratt([xstr ystr]);
                Cx = circ(1); Cy = circ(2); R = circ(3);
                theta0 = unwrap(atan2(ystr-Cy,xstr-Cx));
                % i want 2 km resolution i.e., R * dtheta = 2 km
                dtheta = 2/R;
                theta = min(theta0(:)):dtheta:max(theta0(:));
                xstr = Cx + R .* cos(theta);
                ystr = Cy + R .* sin(theta);
                
                strmask = round(interp2(xr',yr',stream',xstr,ystr));
                xstr(strmask == 0) = [];
                ystr(strmask == 0) = [];
                
                if ~isclockwise(xstr,ystr)
                    xstr = fliplr(xstr);
                    ystr = fliplr(ystr);
                end
                
                ixstr = []; iystr = [];
            else
                % use angleSort on branch centroids to order regions appropriately
                [~,sortcen] = angleSort(icen(sortnum,:), ...
                                [cxind(tind) cyind(tind)],-pi/2);
                sortcen = flipdim(sortcen,1);
                % alternative to above - sortcen code
                %[~,sortang] = sort(meanangle(sortnum),'descend');

                % now actually select the remaining regions and figure out
                % (x,y) co-ordinates
                sortnum = sortnum(sortcen);

                % sortnum should be final sorted order here
                [ixstr, iystr] = ind2sub(size(skel), skelcomps.PixelIdxList{sortnum(1)});
                for mm = 1:length(sortnum)-1
                    ix1 = ixstr(end);
                    iy1 = iystr(end);

                    [ix2,iy2] = ind2sub(size(skel), ...
                                        skelcomps.PixelIdxList{sortnum(mm+1)});

                    % use Bresenham's algorithm to join
                    [jx,jy] = bresenham(ix1,iy1,ix2(1),iy2(1));

                    ixstr = [ixstr; jx; ix2];
                    iystr = [iystr; jy; iy2];
                end

                xstr = xr(ixstr,1)';
                ystr = yr(1,iystr);
            end
            
            % fix the starting!!!
            dstr = [0 cumsum(hypot(diff(xstr),diff(ystr)))];
            
            % distance from perimeter - NOT QUITE AS GOOD
            %{
            distper = bwdist(~stream);
            [~,index1] = max(distper(1:cxind(tt),:));
            [~,index2] = max(distper(cxind(tt):end,:));
            index1(index1 == 1) = NaN;
            index2(index2 == 1) = NaN;
            index2 = index2+cxind(tt);
            idxx = [index1(:); fliplr(index2(:))]';
            idxy = [1:size(stream,2) fliplr(1:size(stream,2))];
            %}
            %polyline = [cut_nan(idxx)' (cut_nan(idxy .* idxx)./cut_nan(idxx))'];
            
            % testing streamer cross-section detection
            if debug_plot
                clf
                subplot(211)
                pcolorcen(xr,yr,double(stream));
                hold on;
                plot(cx(tind),cy(tind),'ko','MarkerSize',16);
                %plot(idxx,idxy,'bx','markersize',8);
                plot(xstr,ystr,'k*');
                
                drawCircle(circ(1),circ(2),circ(3));
                %ell = EllipseDirectFit([ixstr iystr]);
                %a = ell(1); b = ell(2); c = ell(3); 
                %d = ell(4); e = ell(5); f = ell(6);
                %x0 = (c*d - b*f)/(b^2-a*c);
                %y0 = (a*f - b*d)/(b^2-a*c);
                
                title(num2str(tind));
                subplot(212)
                for mm=1:length(sortnum)
                    testBranch(skelcomps.PixelIdxList{sortnum(mm)}, ...
                        size(skel),sortnum(mm));
                end
                set(gca,'ydir','normal');
                hold on;
                plot(cxind(tind),cyind(tind),'ko','MarkerSize',16);
                plot(icen(:,1),icen(:,2),'k*','MarkerSize',8);
                pause();
            end
            
            % save locations in runs object
            runs.streamer.west.xstr{tind}  = xstr;
            runs.streamer.west.ystr{tind}  = ystr;
            runs.streamer.west.ixstr{tind} = ixstr;
            runs.streamer.west.iystr{tind} = iystr;
            runs.streamer.west.dstr{tind}  = dstr;
            runs.streamer.comment   = ...
                [' contour = 1 in streamer, 0 outside | ' ...
                ' (xstr,ystr) = cross-section through streamer (cell array) | ' ...
                ' (ixstr,iystr) = indices corresponding to (xstr,ystr) ' ...
                ' - (cell array) | dstr = along-streamer distance (cell array)'];
        end
    end
    
    % plot streamer profiles
    function [] = plot_streamerstats(runs)
        bins = runs.streamer.bins;
        figure
        subplot(121)
        cmap = brighten(cbrewer('seq','YlOrRd',runs.streamer.sz4dsp(2)),0);
        cmap = cmap(3:end,:,:); % chuck out lightest colors
        set(gca,'ColorOrder',cmap); colormap(cmap);
        line(runs.streamer.west.Vbin, repmat(avg1(bins'),[1 runs.streamer.sz4dsp(2)]));
        hold on
        zcenbin = vecfind(bins,runs.streamer.west.zcen);
        Vcenbin = diag(runs.streamer.west.Vbin(zcenbin,1:slab));
         colorbar; cblabel('day');
        scatter(gca,zeros(20,1),runs.streamer.west.zcen, ...
                    96,runs.streamer.time/86400,'filled');
        caxis([min(runs.streamer.time) max(runs.streamer.time)]/86400);
        xlabel('Volume (m^3)'); 
        ylabel(['Depth (' num2str(dbin) ' m bins)']);

        subplot(122); hold on
        plot(runs.streamer.time/86400,runs.streamer.west.zcen,'r');
        plot(runs.streamer.time/86400,runs.streamer.west.zdcen,'b');
        legend('z centroid','z-dye centroid');
        ylabel(' Depth (m) '); xlabel('day');
    end
    
    % domain integration for sparse matrix input
    function [out] = domain_integratesp(runs,in, dV)

        if ~exist('dV','var') % not good idea
            dV = reshape(runs.rgrid.dV, runs.streamer.sz3dsp);
        end

        out = full(sum( bsxfun(@times, in, dV)));
    end

    % domain integration for full matrix input
    function [out] = domain_integrate(runs,in, dV)

        if ~exist('dV','var'), dV = runs.rgrid.dV; end

        sz = size(in);
        if length(sz) == 3, sz(4) = 1; end
        out = nansum( reshape( bsxfun(@times, in, dV), ...
                [prod(sz(1:end-1)) sz(end)]), 1);
    end

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
            zdye = dc_roms_read_data(runs.dir, 'dye_02', t0+tt-1,volume,[],runs.rgrid);
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

    function [] = distrib_csz_old(runs)
        % lets subtract out mean at each z-level to account for near
        % surface and near bottom upwelling.
        % This has to be done after interpolating to constant z-level
        % because you can't take a constant z-level mean otherwise
        yend = find_approx(runs.rgrid.y_rho(:,1),100*1000);
        t0 = runs.eddy.trevind;
        read_start = [1 1 1 t0];
        read_count = [Inf yend Inf 35];

        zdye = ncread(runs.out_file,'dye_02', ...
                        read_start,read_count);
        csdye = ncread(runs.out_file,'dye_01', ...
                        read_start, read_count)/1000;
        w = ncread(runs.out_file,'w', ...
                        read_start, read_count);
        %asdye = double(ncread(runs.out_file,'dye_03', ...
        %                [1 1 1 runs.eddy.trevind],[Inf Inf Inf 20]))/1000;

        % depth to interpolate to
        depth = 100;
        xsb = runs.bathy.xsb/1000;
        [grd.xax,grd.yax,grd.zax,~,~,~] = dc_roms_var_grid(runs.rgrid,'temp');
        grd.xax = grd.xax(:,1:yend,:);
        grd.yax = grd.yax(:,1:yend,:);
        grd.zax = grd.zax(:,1:yend,:);

        % grid matrices required for plotting        
        xr = runs.rgrid.xr(:,1:yend)/1000; yr = runs.rgrid.yr(:,1:yend)/1000;
        yzw = repmat(yr(1,:)', [1 runs.rgrid.N+1]);
        yzr = repmat(yr(1,:)', [1 runs.rgrid.N]);
        zw = permute(runs.rgrid.z_w(:,1:yend,:),[3 2 1]);

    %             figure;
    %             for tt = 1:size(zdye,4)
    %                 clf;
    %                 tind = runs.eddy.trevind + tt;
    %                 % interpolate to a given depth
    %                 zdyein = dc_roms_zslice_var(zdye(:,:,:,tt),depth,grd);
    %                 csdyein = dc_roms_zslice_var(csdye(:,:,:,tt),depth,grd);
    %                 
    %                 % define streamer
    %                 streamer = fillnan((csdyein > xsb-10) & (csdyein < xsb+30) ...
    %                             & (runs.rgrid.x_rho' < runs.eddy.cx(tind)),0);
    %                 %streamer = fillnan( csdyein < xsb, 0);
    %                 
    %                 % remove mean to show up/down-welling
    %                 zdyein_demean = zdyein - nanmean(zdyein(:));
    %                 
    %                 % visualize
    %                 pcolorcen((zdyein_demean .* streamer)');
    %                 hold on
    %                 contour(runs.eddy.mask(:,:,tind)','k','LineWidth',2);
    %                 pause();
    %             end
    %             
        % mask of points west of eddy center
        %west_mask = bsxfun(@lt,repmat(runs.rgrid.x_rho',[1 1 runs.rgrid.N]), ...
        %               permute(runs.eddy.cx(runs.eddy.trevind:runs.eddy.trevind+19), [1 3 4 2]));


        % identify streamer again, but now with 4D data
        % this is more general compared to streamer2
        streamer = fillnan( (csdye > xsb-10) & (csdye < xsb+30) ...
                       , 0);

        % number of west of eddy's west edge for streamer cross section
        dx = 4;

        % (xs,ys,zs) are the Eulerian x,y,z values
        %xs = bsxfun(@times, streamer, grd.xax)/1000;
        ys = bsxfun(@times, streamer, grd.yax)/1000;
        zs = bsxfun(@times, streamer, grd.zax);

        % (as,cs,z) dyes contain the Lagrangian labels
        % some distance metric between the two will give me an idea of
        % what's happening
        if runs.bathy.axis == 'y'
        %    das = asdye - xs;
            dcs = csdye - ys;
        else
        %    das = asdye - ys;
            dcs = csdye - xs;
        end
        %dz = zdye - zs;

        % make streamer section - with more processing
        % NEED TO ACCOUNT FOR TILTING IN VERTICAL?
        cx = runs.eddy.cx(t0:t0+read_count(end)-1)/1000;
        cy = runs.eddy.cy(t0:t0+read_count(end)-1)/1000;
        ee = runs.eddy.ee(t0:t0+read_count(end)-1)/1000;
        % hack if eddy center is outside extracted domain
        cy(cy > max(yr(:))) = max(yr(:));
        cxind = vecfind(xr(:,1),cx);
        cyind = vecfind(yr(1,:),cy)';

        %r = sqrt(bsxfun(@minus,xr,permute(cx,[3 1 2])).^2 ...
        %       + bsxfun(@minus,yr,permute(cy,[3 1 2])).^2);
        % picking only western streamer
        streamer2 = squeeze(streamer(:,:,end,:)  ... % streamer
                    ... % parcels have moved more than 5 km in the cross-shelf dirn.
                         .* (abs(dcs(:,:,end,:))>5)) ...
                    ... % remove eastern half
                         .* ( bsxfun(@lt, xr, ...
                              permute(ee + runs.params.eddy.dia/2000,[3 1 2])));

         stream = repnan(streamer2(:,:,end),0);

    %             streamer2 = fillnan(streamer2,0); 
    %             xs2 = bsxfun(@times, streamer2, xr);
    %             ys2 = bsxfun(@times, streamer2, yr);
    %             rstreamer = r .* streamer2;
    %             find mean r in along-stream direction.
    %             rs = squeeze(nanmean(rstreamer,1));
    %             
    %             divide streamer into E-W & N-S halves to account for
    %             multiple valued contour
    %             for tt = 1:size(streamer2,3)
    %                 xsect = [squeeze(nanmean(xs2(1:cxind,1:end,tt),1)) ...
    %                          ... %cut_nan(squeeze(nanmean(xs2(1:cxind,cyind+1:end,tt),1))) ...
    %                          squeeze(nanmean(xs2(cxind+1:end,1:end,tt),1))];% ...
    %                          ...%cut_nan(squeeze(nanmean(xs2(cxind+1:end,cyind+1:end,tt),1)))]; 
    %                 ysect = fillnan(~isnan(xsect),0) .* [yr(1,:) yr(1,:)];     
    %                 xsect = cut_nan(xsect);
    %                 ysect = cut_nan(ysect);
    %                 ysect = [cut_nan(squeeze(nanmean(ys2(1:cxind,1:cyind,tt),2)))' ...
    %                        cut_nan(squeeze(nanmean(ys2(1:cxind,cyind+1:end,tt),2)))' ...
    %                         cut_nan(squeeze(nanmean(ys2(cxind+1:end,1:cyind,tt),2)))' ...
    %                         cut_nan(squeeze(nanmean(ys2(cxind+1:end,cyind+1:end,tt),2)))']; 
    %             end


        % vertically integrated w in streamer
        WS = squeeze( nansum( bsxfun(@times,avg1(w,3).*streamer, diff(zw,1,3) ) ...
                                , 3) );

        % index of western & eastern edges
        wind = vecfind(xr(:,1), runs.eddy.vor.we/1000);
        eind = vecfind(xr(:,1), runs.eddy.vor.ee/1000);

        % colorbar for vertical vel cross-section
        wcolor = sort( [-1 1] * max(max(abs( ...
                            log10(abs(w(sort([eind wind]),:))) ))) )/2;

        figure; 
       %% animate depth integrated w in streamer
        figure;clf; ii=1; maximize();
        %subplot(231); subplot(232); subplot(233);
        %subplot(234); subplot(235); subplot(236);
        %spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))
        tindex = t0+ii-1;
        windex = wind(tindex)-dx; % for cross-section
        eindex = eind(tindex)-dx; % for cross-section
        zlimit = [-1000 0];

        subplot(231)
        titlestr = 'Depth integrated w in streamer (blue)';
        hws = pcolorcen(xr,yr,double(WS(:,:,ii))); shading flat;
        hxw = linex(xr(windex,1));
        hxe = linex(xr(eindex,1));
        hold on;
        [~,hs] = contour(xr,yr,repnan(streamer(:,:,40,ii),0), ...
                        1,'b','LineWidth',2);
        he = runs.plot_eddy_contour('contour',tindex);
        runs.plot_bathy('contour','k');
        colormap(flipud(cbrewer('div','RdBu',32)));
        caxis([-1 1] * max(abs([nanmin(WS(:)) nanmax(WS(:))])));
        colorbar; %cbunits('m^2/s'); 
        ht = runs.set_title(titlestr,tindex);


        % depth of 'streamer'
        subplot(234)
        hz = pcolorcen(xr,yr,double(max(abs(zs(:,:,:,ii)),[],3)));
        hold on;
        hcb = colorbar;  caxis([0 max(abs(zs(:)))]);cbunits('[m]');
        hzeta = runs.plot_zeta('contour',tindex);
        title('Depth of ''streamer''');
        %colormap(cbrewer('seq','Blues',32));
        %xlim([100 400]); 
        %ylim([0 140]);

        % vertical vel - west cross-section
        subplot(232)
        [hwcs] = pcolorcen(yzw,squeeze(zw(1,:,:)), ...
                        double(squeeze( ...
                            sign(w(windex,:,:,ii)) .* log10(abs(w(windex,:,:,ii))) ...
                                )));
        colorbar; caxis(wcolor);
        ylim(zlimit); linex(runs.bathy.xsl/1000,'slopebreak','w');
        hcw = linex(runs.eddy.vor.cy(tindex)/1000);
        ylabel('Z (m)'); xlabel('Y (km)');
        title('cross-section of log_{10}(w)');

        % vertical vel - east cross-section
        subplot(235)
        [hecs] = pcolorcen(yzw,squeeze(zw(1,:,:)), ...
                        double(squeeze( ...
                        sign(w(eindex,:,:,ii)) .* log10(abs( w(eindex,:,:,ii) )) ...
                        )));
        colorbar; caxis(wcolor);
        ylim(zlimit); 
        linex(runs.bathy.xsl/1000,'slopebreak','w');
        hce = linex(runs.eddy.vor.cy(tindex)/1000);
        ylabel('Z (m)'); xlabel('Y (km)');
        title('cross-section of log_{10}(w)');


        % z-dye - west cross-section
        subplot(233)
        [~,hwz] = contourf(yzr,squeeze(grd.zax(1,:,:)), ...
                        double(squeeze(zdye(windex,:,:,ii))), ...
                        linspace(zlimit(1),zlimit(2),20));
        colorbar; caxis(zlimit);
        ylim(zlimit); linex(runs.bathy.xsl/1000,'slopebreak','w');
        hcw = linex(runs.eddy.vor.cy(tindex)/1000);
        ylabel('Z (m)'); xlabel('Y (km)');
        title('cross-section of z-dye | need to adjust for BC');

        % zdye - east cross-section
        subplot(236)
        [~,hez] = contourf(yzr,squeeze(grd.zax(1,:,:)), ...
                        double(squeeze( zdye(eindex,:,:,ii) )), ....
                        linspace(zlimit(1),zlimit(2),20));
        colorbar; caxis(zlimit);
        ylim(zlimit); 
        linex(runs.bathy.xsl/1000,'slopebreak','w');
        hce = linex(runs.eddy.vor.cy(tindex)/1000);
        ylabel('Z (m)'); xlabel('Y (km)');
        title('cross-section of z-dye | need to adjust for BC');

        spaceplots(0.05*ones([1 4]),0.04*ones([1 2]));
        pause();         

        for ii=2:size(WS,3)
            tindex = t0+ii-1;
            % for cross-section
            windex = wind(tindex)-dx; 
            eindex = eind(tindex)+dx;

            set(hws ,'CData',double(WS(:,:,ii)));
            set(hs  ,'ZData',repnan(streamer(:,:,40,ii),0));

            set(hz  ,'CData',double(max(abs(zs(:,:,:,ii)),[],3)));

            set(hwcs,'CData',double(squeeze( ...
                sign(w(eindex,:,:,ii)) .* log10(abs( w(windex,:,:,ii) )) )));
            set(hecs,'CData',double(squeeze( ...
                sign(w(eindex,:,:,ii)) .* log10(abs( w(eindex,:,:,ii))) )));

            set(hwz ,'ZData',double(squeeze( zdye(windex,:,:,ii) )));
            set(hez ,'ZData',double(squeeze( zdye(eindex,:,:,ii) )));

            set(hxw ,'XData',[1 1]*xr(windex,1));
            set(hxe ,'XData',[1 1]*xr(eindex,1));
            set(hce ,'XData',[1 1]*runs.eddy.vor.cy(tindex)/1000);
            set(hcw ,'XData',[1 1]*runs.eddy.vor.cy(tindex)/1000);

            runs.update_zeta(hzeta,tindex);
            runs.update_eddy_contour(he, tindex);
            runs.update_title(ht,titlestr,tindex);
            pause();
        end

        %% old stuff
        %             figure;
%             for tt = 1:size(zdye,4)
%                 clf;
%                 tind = runs.eddy.trevind + tt;
%                 % interpolate to a given depth
%                 zdyein = dc_roms_zslice_var(zdye(:,:,:,tt),depth,grd);
%                 csdyein = dc_roms_zslice_var(csdye(:,:,:,tt),depth,grd);
%                 
%                 % define streamer
%                 streamer = fillnan((csdyein > xsb-10) & (csdyein < xsb+30) ...
%                             & (runs.rgrid.x_rho' < runs.eddy.cx(tind)),0);
%                 %streamer = fillnan( csdyein < xsb, 0);
%                 
%                 % remove mean to show up/down-welling
%                 zdyein_demean = zdyein - nanmean(zdyein(:));
%                 
%                 % visualize
%                 pcolorcen((zdyein_demean .* streamer)');
%                 hold on
%                 contour(runs.eddy.mask(:,:,tind)','k','LineWidth',2);
%                 pause();
%             end
%             

        %% another streamer attempt

%             streamer2 = fillnan(streamer2,0); 
%             xs2 = bsxfun(@times, streamer2, xr);
%             ys2 = bsxfun(@times, streamer2, yr);
%             rstreamer = r .* streamer2;
%             find mean r in along-stream direction.
%             rs = squeeze(nanmean(rstreamer,1));
%             
%             divide streamer into E-W & N-S halves to account for
%             multiple valued contour
%             for tt = 1:size(streamer2,3)
%                 xsect = [squeeze(nanmean(xs2(1:cxind,1:end,tt),1)) ...
%                          ... %cut_nan(squeeze(nanmean(xs2(1:cxind,cyind+1:end,tt),1))) ...
%                          squeeze(nanmean(xs2(cxind+1:end,1:end,tt),1))];% ...
%                          ...%cut_nan(squeeze(nanmean(xs2(cxind+1:end,cyind+1:end,tt),1)))]; 
%                 ysect = fillnan(~isnan(xsect),0) .* [yr(1,:) yr(1,:)];     
%                 xsect = cut_nan(xsect);
%                 ysect = cut_nan(ysect);
%                 ysect = [cut_nan(squeeze(nanmean(ys2(1:cxind,1:cyind,tt),2)))' ...
%                        cut_nan(squeeze(nanmean(ys2(1:cxind,cyind+1:end,tt),2)))' ...
%                         cut_nan(squeeze(nanmean(ys2(cxind+1:end,1:cyind,tt),2)))' ...
%                         cut_nan(squeeze(nanmean(ys2(cxind+1:end,cyind+1:end,tt),2)))']; 
%             end
    end

   %% animation functions
    function [] = animate_zeta(runs)
        runs.video_init('zeta');

        titlestr = 'SSH (m)';

        figure;
        ii=1;
        hz = runs.plot_zeta('pcolor',ii);
        ax = gca;
        hold on
        colorbar; freezeColors;
        hbathy = runs.plot_bathy('contour','k');
        he = runs.plot_eddy_contour('contour',ii);
        ht = runs.set_title(titlestr,ii);

        xlabel('X (km)');ylabel('Y (km)');
        axis image;
        maximize(gcf); pause(0.2);  
        beautify([16 16 18]);
        runs.video_update();
        for ii = 2:size(runs.zeta,3)
            runs.update_zeta(hz,ii);
            runs.update_eddy_contour(he,ii);
            runs.update_title(ht,titlestr,ii);
            runs.video_update();
            pause(0.03);
        end
        runs.video_write();
    end
    
    function [] = animate_streamer_section(runs)
        
        debug_plot = 0;
        if ~isfield(runs.streamer.west,'mask')
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
        % co-ordinate axes
        
        %[grd.xax,grd.yax,grd.zax,~,~,~] = dc_roms_var_grid(runs.rgrid,'temp');
        %grd.xax = grd.xax(:,1:yend,:);
        %grd.yax = grd.yax(:,1:yend,:);
        %grd.zax = grd.zax(:,1:yend,:);
        
        % grid matrices required for plotting
        xr = runs.rgrid.xr(:,1:yend)/1000; yr = runs.rgrid.yr(:,1:yend)/1000;
        %ix = repmat([1:size(xr,1)]',[1 yend]);
        %iy = repmat([1:yend],[size(xr,1) 1]);
        %yzw = repmat(yr(1,:)', [1 runs.rgrid.N+1]);
        %yzr = repmat(yr(1,:)', [1 runs.rgrid.N]);
        zw = permute(runs.rgrid.z_w(:,1:yend,:),[3 2 1]);
        zr = permute(runs.rgrid.z_r(:,1:yend,:),[3 2 1]);
        
        % NEED TO ACCOUNT FOR TILTING IN VERTICAL?
        %cx = runs.eddy.cx(t0:tend)/1000;
        %cy = runs.eddy.cy(t0:tend)/1000;
        %ee = runs.eddy.ee(t0:tend)/1000;
        % hack if eddy center is outside extracted domain
        %cy(cy > max(yr(:))) = max(yr(:));
        %cxind = vecfind(xr(:,1),cx);
        %cyind = vecfind(yr(1,:),cy)';
        
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
            dstr = repmat(runs.streamer.west.dstr{tind}',[1 runs.rgrid.N]);
            
            if runs.streamer.west.fit_circle
                bstr = interp2(xr',yr',runs.bathy.h(:,1:yend)',xstr,ystr);
            else
                ixstr = runs.streamer.west.ixstr{tind};
                iystr = runs.streamer.west.iystr{tind};
                indices = sub2ind(sz4dfull(1:2),ixstr,iystr);

                zlin = reshape(zr,sz3d2d);
                zstr = zlin(indices,:);
                clear zlin;
                
                bstr = runs.bathy.h(indices)';
            end
            
            ixmin = find_approx(xr(:,1),min(xstr));
            ixmax = find_approx(xr(:,1),max(xstr));
            iymin = find_approx(yr(1,:),min(ystr));
            iymax = find_approx(yr(1,:),max(ystr));
            
            % bathy-patch
            bpatch = [bstr -max(runs.bathy.h(:))-100 ...
                                    -max(runs.bathy.h(:))-100];
            dpatch = [dstr(:,1)' dstr(end,1) 0];
                        
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
            u = dc_roms_read_data(runs.dir,'u', ...
                tind,volume,[],runs.rgrid);
            v = dc_roms_read_data(runs.dir,'v', ...
                tind,volume,[],runs.rgrid);
            csdye = dc_roms_read_data(runs.dir, 'dye_01', ...
                tind,volume,[],runs.rgrid);
            zdye = dc_roms_read_data(runs.dir, 'dye_02', ...
                tind,volume,[],runs.rgrid);
                  
            if runs.streamer.west.fit_circle
                strstr = 
            else
                % streamer mask vertical section - along-streamer section
                % points
                strlin = reshape(streamer,sz3d2d);
                strstr = strlin(indices,:);

                ixnew = ixstr - ixmin + 1;
                iynew = iystr - iymin + 1;
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
            Utstr(strstr == 0) = NaN;
            Unstr(strstr == 0) = NaN;
            zdstr(strstr == 0) = NaN;
            csstr(strstr == 0) = NaN;
            
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
                hpatch(3) = patch(dpatch,bpatch,'k');
                
                % cross-shelf dye - streamer section
                ax(4) = subplot(235);
                [~,hcsd] = contourf(dstr,zstr,csstr/1000 - runs.bathy.xsb/1000);
                colorbar; ylim(limz); caxis([-10 40]);
                ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
                title('Cross-shelf dye - X_{shelfbreak} (km)');
                hold on;
                hpatch(4) = patch(dpatch,bpatch,'k');
                
                % velocities - streamer section
                ax(5) = subplot(233);
                [~,hun] = contourf(dstr,zstr,Unstr);
                colorbar; ylim(limz);caxis([-1 1]*0.1);
                ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
                title('Normal velocity (m/s)');
                hold on;
                hpatch(5) = patch(dpatch,bpatch,'k');
                
                ax(6) = subplot(236);
                [~,hut] = contourf(dstr,zstr,Utstr);
                colorbar; ylim(limz); caxis([-1 1]*0.1);
                ylabel('Z (m)'); xlabel('Along-streamer dist (km)');
                title('Tangential velocity (m/s)');
                hold on;
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
                set(hs2 ,'ZData',repnan(streamer,0));
                
                set(hzdye,'XData',dstr,'YData',zstr,'ZData', zdstr - zstr);
                set(hcsd,'XData',dstr,'YData',zstr,'ZData', ...
                    csstr/1000 - runs.bathy.xsb/1000);
                
                set(hun ,'XData',dstr,'YData',zstr,'ZData',Unstr);
                set(hut ,'XData',dstr,'YData',zstr,'ZData',Utstr);
                
                %runs.update_zeta(hzeta,tindex);
                runs.update_eddy_contour(he2, tindex);
                runs.update_eddy_contour(he, tindex);
                runs.update_title(ht,titlestr,tindex);
                pause();
            end
        end
    end
    
    function [] = animate_3d(runs)
        stride = [1 1 1 1];

        xrmat = repmat(runs.rgrid.xr(1:stride(1):end,1:stride(2):end)', ...
                        [1 1 runs.rgrid.N]);
        yrmat = repmat(runs.rgrid.yr(1:stride(1):end,1:stride(2):end)', ...
                        [1 1 runs.rgrid.N]);
        zrmat = permute(runs.rgrid.zr(1:stride(1):end,1:stride(2):end,:),[2 1 3]);

        %eddye = roms_read_data(runs.dir,'dye_04',[1 1 1 1],[Inf Inf Inf Inf], ...
        %            stride);

        tic;csdye = ncread(runs.out_file,'dye_01');toc;
        tic;eddye = ncread(runs.out_file,'dye_04');toc;
        csdye = permute(csdye,[2 1 3 4]);
        eddye = permute(eddye,[2 1 3 4]);
        mask = zeros(size(eddye,2),size(eddye,1),size(eddye,4));
        mask(2:end-1,2:end-1,:)=runs.eddy.vormask;
        mask = ones(size(mask));

        %% make isosurface plot

        eddlevel = zrmat; %0.8;
        thresh = 0.8;
        xsb = runs.bathy.xsb/1000;
        cslevel = [xsb-10 xsb]*1000;
        sbcolors = distinguishable_colors(length(cslevel));

        clf; clear pcsd pedd;
        hold on
        hbathy = surf(runs.rgrid.xr/1000,runs.rgrid.yr/1000,-runs.bathy.h);
        colormap(copper); freezeColors;
        set(hbathy,'FaceColor','Flat','EdgeColor','None');

        ii=1;
        [faces,verts,colors] = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                bsxfun(@times,eddye(:,:,:,1) > thresh,mask(:,:,1)'),eddlevel);
        pedd = patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors, ...
                'FaceColor','interp','EdgeColor','none');
        colormap(flipud(cbrewer('div', 'RdYlGn', 32))); freezeColors;
        colorbar; cbfreeze;
        %set(pedd,'EdgeColor','none','FaceAlpha',0.5);
        view(3)

        for kk=1:length(cslevel)
            pcsd(kk) = patch(isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                            csdye(:,:,:,ii),cslevel(kk)));
            set(pcsd(kk),'FaceColor',sbcolors(kk,:));
            set(pcsd(kk),'EdgeColor','none');
            %reducepatch(pcsd(kk),0.5,'verbose');
        end
        [~,hedd] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.zeta(:,:,1));

        titlestr = 'dyes';
        ht = runs.set_title(titlestr,ii);
        %view(-104,30);
        view(-150,66);
        xlim([min(xrmat(:)) max(xrmat(:))]/1000)
        xlabel('X'); ylabel('Y'); zlabel('Z');
        pause();
        for ii=2:4:size(eddye,4)
            heddye = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                bsxfun(@times,eddye(:,:,:,ii) > thresh,mask(:,:,ii)'),eddlevel);
            set(pedd,'Vertices',heddye.vertices,'Faces',heddye.faces, ...
                    'FaceVertexCData',heddye.facevertexcdata);
            set(hedd,'ZData',runs.zeta(:,:,ii));
            for kk=1:length(cslevel)
                hcsdye = isosurface(xrmat/1000,yrmat/1000,zrmat, ...
                            csdye(:,:,:,ii),cslevel(kk));
                set(pcsd(kk),'Vertices',hcsdye.vertices,'Faces',hcsdye.faces);
            end
            runs.update_title(ht,titlestr,ii);
            pause();
        end

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
%             if ~exist('tind','var')
%                 tind = [];
%             end
%             if ~exist([runs.dir '/ocean_vor.nc'],'file')
%                 dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
%             end

        tt = 1;
        rvor = double(ncread(runs.out_file,'rvorticity',[1 1 1 tt], ...
            [Inf Inf Inf 1]));
        [rint,ravg] = roms_depthIntegrate(rvor, ...
            runs.rgrid.Cs_r,runs.rgrid.Cs_w, ...
            avg1(avg1(runs.bathy.h,1),2),avg1(avg1(runs.zeta(:,:,tt),1),2), ...
            [0 -max(runs.bathy.h(:))]);

        rplot = ravg;

        figure;
        titlestr = 'Depth avg rvor';
        xvor = avg1(avg1(runs.rgrid.xr,1),2);
        yvor = avg1(avg1(runs.rgrid.yr,1),2);
        hvor = pcolor(xvor/1000,yvor/1000,rplot); hold on; shading flat;
        ht = runs.set_title(titlestr,tt);
        he = runs.plot_eddy_contour('contour',tt);
        hbathy = runs.plot_bathy('contour','k');
        shading flat
        caxis([-1 1] * max(abs(rplot(:)))); colorbar;

        for tt=2:2:size(runs.zeta,3)
            rvor = double(ncread(runs.out_file,'rvorticity',[1 1 1 tt], ...
                [Inf Inf Inf 1]));
            tic;
            [rint,ravg] = roms_depthIntegrate(rvor, ...
                runs.rgrid.Cs_r,runs.rgrid.Cs_w, ...
                avg1(avg1(runs.bathy.h,1),2),avg1(avg1(runs.zeta(:,:,tt),1),2), ...
                [0 -max(runs.bathy.h(:))]);
            rplot = ravg;
            set(hvor,'cdata',rplot);
            runs.update_eddy_contour(he,tt);

            runs.update_title(ht,titelstr,tt);
            toc;
            pause(0.1);
        end

    end

    function [] = animate_vorbudget(runs,tind)
%             if ~exist('tind','var')
%                 tind = [];
%             end
%             if ~exist([runs.dir '/ocean_vor.nc'],'file')
%                 dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
%             end
        % prepare grids for differentiation
        xvor = avg1(avg1(runs.rgrid.xr,1),2);
        yvor = avg1(avg1(runs.rgrid.yr,1),2);
        N = runs.rgrid.N;

        gridu.xmat = repmat(runs.rgrid.x_u',[1 1 N]);
        gridu.ymat = repmat(runs.rgrid.y_u',[1 1 N]);
        gridu.zmat = permute(runs.rgrid.z_u,[3 2 1]);
        gridu.s = runs.rgrid.s_rho;
        gridu.zw = runs.rgrid.z_w;
        gridu.s_w = runs.rgrid.s_w;

        gridv.xmat = repmat(runs.rgrid.x_v',[1 1 N]);
        gridv.ymat = repmat(runs.rgrid.y_v',[1 1 N]);
        gridv.zmat = permute(runs.rgrid.z_v,[3 2 1]);
        gridv.s = runs.rgrid.s_rho;
        gridv.zw = runs.rgrid.z_w;
        gridv.s_w = runs.rgrid.s_w;

        gridw.xmat = repmat(runs.rgrid.x_rho',[1 1 N+1]);
        gridw.ymat = repmat(runs.rgrid.y_rho',[1 1 N+1]);
        gridw.zmat = permute(runs.rgrid.z_w,[3 2 1]);
        gridw.s = runs.rgrid.s_w;
        gridw.zw = runs.rgrid.z_w;
        gridw.s_w = runs.rgrid.s_w;

        gridrv.xmat = repmat(xvor,[1 1 N-1]);
        gridrv.ymat = repmat(yvor,[1 1 N-1]);
        gridrv.zmat = avg1(avg1(avg1(permute(runs.rgrid.z_r,[3 2 1]),1),2),3);
        gridrv.s = avg1(runs.rgrid.s_rho);
        gridrv.zw = avg1(avg1(avg1(permute(runs.rgrid.z_w,[3 2 1]),1),2),3);
        gridrv.s_w = avg1(runs.rgrid.s_w);

        beta = runs.params.phys.beta;

        % for depth integration
        h = runs.bathy.h(2:end-1,2:end-1);
        csr = runs.rgrid.Cs_r(2:end-1); 
        csw = runs.rgrid.Cs_w(2:end-1);

        xavg = avg1(avg1(xvor,1),2)/1000; yavg = avg1(avg1(yvor,1),2)/1000;

        depthRange = [100 -max(runs.bathy.h(:))];
        trange = 1:2:size(runs.zeta,3);

        disp(['starting from t instant = ' num2str(trange(1))]);
        for kk=1:length(trange)
            tt = trange(kk);
            zeta = runs.zeta(2:end-1,2:end-1,tt);

            % read data
            u1 = double(ncread(runs.out_file,'u',[1 1 1 tt],[Inf Inf Inf 2]));
            v1 = double(ncread(runs.out_file,'v',[1 1 1 tt],[Inf Inf Inf 2]));
            w  = double(ncread(runs.out_file,'w',[1 1 1 tt],[Inf Inf Inf 1]));
            %rvor1 = double(ncread(runs.out_file,'rvorticity',[1 1 1 tt], ...
            %    [Inf Inf Inf 2]));
            %rvor = rvor1(:,:,:,1);
            u = u1(:,:,:,1); v = v1(:,:,:,1);
            u1(:,:,:,1) = []; v1(:,:,:,1) = [];

            % differentiate
            ux = diff_cgrid(gridu,u,1); uy = diff_cgrid(gridu,u,2);
                uz = diff_cgrid(gridu,u,3);
            vx = diff_cgrid(gridv,v,1); vy = diff_cgrid(gridv,v,2); 
                vz = diff_cgrid(gridv,v,3);
            wx = avg1(diff_cgrid(gridw,w,1),3); wy = avg1(diff_cgrid(gridw,w,2),3);
                wz = avg1(diff_cgrid(gridw,w,3),3);

            % calculate relative vorticity
            v1x = diff_cgrid(gridv,v1,1); u1y = diff_cgrid(gridu,u1,2);
            rvor = vx-uy; rv1 = v1x-u1y;
            rvx = diff_cgrid(gridrv,rvor,1); rvy = diff_cgrid(gridrv,rvor,2);
                rvz = diff_cgrid(gridrv,rvor,3);


            % check continuity
            cont = ux(:,2:end-1,:) + vy(2:end-1,:,:) + wz(2:end-1,2:end-1,:);
            contfrac = cont./wz(2:end-1,2:end-1,:);
            [~,CONT] = roms_depthIntegrate(cont,csr,csw,h,zeta,depthRange);

            % form terms - avg to interior RHO points
            adv = avg1( avg1(avg1(u(:,:,2:end-1),1),2) .* rvx,2) + ...
                    avg1( avg1(avg1(v(:,:,2:end-1),1),2) .* rvy,1) + ...
                        avg1(avg1( avg1(avg1(avg1(w(:,:,2:end-1),1),2),3) ...
                               .* rvz    ,1),2);
            str = avg1(avg1( ...
                    avg1(   bsxfun(@plus,rvor,avg1(avg1(runs.rgrid.f',1),2)) ...
                            .* avg1(avg1(wz,1),2)   ,1) ...
                                ,2),3);

            bet = avg1(beta * v(2:end-1,:,2:end-1),2);

            tilt = avg1( avg1(avg1(avg1(wy,1).*avg1(uz,2) - ...
                        avg1(wx,2).*avg1(vz,1),1),2),3);

            tend = avg1(avg1( ...
                    avg1(rv1-rvor,3)./diff(runs.rgrid.ocean_time(1:2)) ,1),2);

            % depth integrate
            [ubar,vbar] = uv_barotropic(u,v,runs.rgrid.Hz);
            [rint,ravg] = roms_depthIntegrate(avg1(avg1(rvor,1),2), ...
                            csr,csw,h,zeta,depthRange);
            [ADV ,~] = roms_depthIntegrate(adv ,csr,csw,h,zeta,depthRange);
            [STR ,~] = roms_depthIntegrate(str ,csr,csw,h,zeta,depthRange);
            [BET ,~] = roms_depthIntegrate(bet ,csr,csw,h,zeta,depthRange);
            [TILT,~] = roms_depthIntegrate(tilt,csr,csw,h,zeta,depthRange);
            [TEND,~] = roms_depthIntegrate(tend,csr,csw,h,zeta,depthRange);
            rplot = rint; 

            budget = tend+adv+bet-tilt-str;
            BUD = TEND+ADV+BET-TILT-STR;

            ubar = avg1(ubar(:,2:end-1),1);
            vbar = avg1(vbar(2:end-1,:),2);

            limy = [0 150];
            titlestr = 'Depth avg rvor';
            % plot
            if kk == 1
                figure;
                ax(1) = subplot(2,4,[1:2]);
                hvor = pcolor(xavg,yavg,rplot); hold on; shading flat;
                axis image;
                ht = runs.set_title('Depth avg rvor',tt);
                he(1) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                shading flat
                caxis([-1 1] * max(abs(rplot(:)))); colorbar;
                ylim(limy);

                ax(2) = subplot(2,4,3);
                hbet = pcolor(xavg,yavg,-BET); colorbar; shading flat;
                he(2) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                caxis([-1 1] * max(abs(BET(:))));
                title('- \beta V');

                ax(3) = subplot(2,4,4); cla
                xran = 1:6:size(xavg,1); yran = 1:4:size(yavg,2);
                hquiv = quiver(xavg(xran,yran),yavg(xran,yran), ...
                        ubar(xran,yran), vbar(xran,yran),1.5);
                title('(ubar,vbar)');
                he(3) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');

                ax(4) = subplot(2,4,5);
                htend = pcolor(xavg,yavg,TEND); colorbar; shading flat;
                he(4) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                caxis([-1 1] * max(abs(TEND(:))));
                title('d\xi/dt');

                ax(5) = subplot(2,4,6);
                hadv = pcolor(xavg,yavg,-ADV); colorbar; shading flat;
                he(5) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                caxis([-1 1] * max(abs(ADV(:))));
                title('-Advection');

                ax(6) = subplot(2,4,7);
                htilt = pcolor(xavg,yavg,TILT); colorbar; shading flat;
                he(6) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                caxis([-1 1] * max(abs(TILT(:))));
                title('Tilting');

                ax(7) = subplot(2,4,8);
                hstr = pcolor(xavg,yavg,STR); colorbar; hold on; shading flat;
                %hquiv = quiverclr(xavg(xran,yran),yavg(xran,yran), ...
                %    ubar(xran,yran),vbar(xran,yran),0.3,STR(xran,yran), ...
                %    [-1 1]*1e-11);
                %set(gca,'color',[0 0 0]);
                he(7) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                caxis([-1 1] * max(abs(STR(:))));
                title('Stretching = (f+\xi)w_z')
                pause();
                spaceplots(0.06*ones([1 4]),0.05*ones([1 2]))
                linkaxes(ax,'xy');
            else
                set(hvor ,'cdata',rplot); 
                set(hadv ,'cdata',-ADV);
                set(hbet ,'cdata',-BET);
                set(hstr ,'cdata',STR);
                set(htilt,'cdata',TILT);
                set(htend,'cdata',TEND);
                try
                    set(hquiv,'udata',ubar(xran,yran),'vdata',vbar(xran,yran));
                catch ME
                end

                runs.update_eddy_contour(he,tt);
                runs.update_title(ht,titlestr,tt);
                pause();
            end
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

        if runs.bathy.axis == 'x'
            stride = [sxy 1 sz st];
            temper = dc_roms_read_data(runs.dir,'temp',[t0 st Inf], ...
                            {'y' iymin iymax},stride);
            strat = dc_roms_read_data(runs.dir,'temp',[1 1], ...
                            {'y' Inf Inf},stride);

            temper = bsxfun(@minus,temper,permute(strat,[1 3 2]));
            %temper = roms_read_data(runs.dir,'temp',[1 iymin 1 t0], ...
            %                  ceil([Inf iymax-iymin+1 Inf Inf]./stride), stride);
            %              toc;
            %strat  = roms_read_data(runs.dir,'temp',[Inf 1 1 1], ...
            %                  ceil([1 1 Inf 1]./stride),stride);
            %              toc

        else
            stride = [1 sxy sz st];
            temper = dc_roms_read_data(runs.dir,'temp',[t0 st Inf], ...
                            {'x' ixmin ixmax},stride);
            strat = dc_roms_read_data(runs.dir,'temp',[1 1], ...
                            {'y' Inf Inf},stride);
            temper = bsxfun(@minus,temper,permute(strat,[3 1 2]));
                        %temper = roms_read_data(runs.dir,'temp',[ixmin 1  1 t0], ...
            %                ceil([ixmax-ixmin+1 Inf Inf Inf]./stride),stride);
            %            toc;
            %strat  = roms_read_data(runs.dir,'temp',[1 1 1 1], ...
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

    function [] = animate_pt(runs,depth,t0)

        if ~exist('depth','var'), depth = 0; end
        if ~exist('t0','var'), t0 = 1; end

        runs.video_init(['pt-z-' num2str(abs(depth))]);

        %dye = csdye/1000;
        rr = sqrt(runs.params.phys.N2)*runs.bathy.hsb/runs.rgrid.f(runs.bathy.isb,1);
        distance = 5*rr; % 5 times rossby radius

        cmedd = cbrewer('seq','Greys',32);%flipud(cbrewer('div', 'RdYlGn', 32));
        cmcsd = haxby;
        cmcsd = cmcsd(1:end-3,:,:);
        clim_edd = [0 1];
        clim_csd = [0 runs.bathy.xsb/1000 + 50];

        % stride for quiver
        dxi = 5; dyi = 3;

        figure;
        i = t0;
        if depth == 0
            if isempty(runs.usurf) || isempty(runs.vsurf)
                runs.read_velsurf;
            end
            dye = runs.eddye(:,:,i);
            csdye = runs.csdye(:,:,i);
            u = runs.usurf(1:dxi:end,1:dyi:end,i);
            v = runs.vsurf(1:dxi:end,1:dyi:end,i);
        else
            grdr.xax = repmat(runs.rgrid.xr,[1 1 runs.rgrid.N]);
            grdr.yax = repmat(runs.rgrid.yr,[1 1 runs.rgrid.N]);
            grdr.zax = permute(runs.rgrid.z_r,[3 2 1]);

            grdu.xax = repmat(runs.rgrid.x_u',[1 1 runs.rgrid.N]);
            grdu.yax = repmat(runs.rgrid.y_u',[1 1 runs.rgrid.N]);
            grdu.zax = permute(runs.rgrid.z_u,[3 2 1]);

            grdv.xax = repmat(runs.rgrid.x_v',[1 1 runs.rgrid.N]);
            grdv.yax = repmat(runs.rgrid.y_v',[1 1 runs.rgrid.N]);
            grdv.zax = permute(runs.rgrid.z_v,[3 2 1]);

            % read and interpolate
            disp(['Reading and interpolating ' num2str(i)]);
            dye = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,'dye_04',i), depth,grdr);
            csdye = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,'dye_01',i), depth,grdr);
            u = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,'u',i), depth,grdu);
            v = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,'v',i), depth,grdv);
            % get on interior RHO points
            u = avg1(u(:,2:end-1),1);
            v = avg1(v(2:end-1,:),2);
            % decimate for quiver
            u = u(1:dxi:end,1:dyi:end);
            v = v(1:dxi:end,1:dyi:end);
        end
        % first get z-slice out
        heddye = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000, ...
                    -addnan(-dye,-0.1));

        ylim([0 130]);
        shading flat;
        caxis(clim_edd);colormap(cmedd);freezeColors;
        hcb1 = colorbar; cbunits(hcb1,'Eddy');cbfreeze(hcb1);
        hold on
        he = runs.plot_eddy_contour('contour',i);

        hcsdye = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000, ...
                    fillnan((csdye/1000 < clim_csd(2)) .* csdye/1000,0));
        shading flat;
        caxis(clim_csd);colormap(cmcsd); freezeColors; 
        hcb2 = colorbar; cbunits(hcb2,'Cross-shore dye');cbfreeze(hcb2);

        hq = quiver(runs.eddy.xr(1:dxi:end,1:dyi:end)/1000,runs.eddy.yr(1:dxi:end,1:dyi:end)/1000, ...
                    u,v);
        set(he,'LineWidth',2);    
        hbathy = runs.plot_bathy('Contour','k');
        titlestr = ['CS dye | z = ' num2str(depth) 'm'];
        ht = runs.set_title(titlestr,i);
        xlabel('X (km)');ylabel('Y (km)');
        %axis image;
        beautify;
        pause();
        runs.video_update();
        for i = t0+1:size(runs.zeta,3)
            if depth == 0
                dye = runs.eddye(:,:,i);
                csdye = runs.csdye(:,:,i);
                u = runs.usurf(1:dxi:end,1:dyi:end,i);
                v = runs.vsurf(1:dxi:end,1:dyi:end,i);
            else
                % read and interpolate
                disp(['Reading and interpolating ' num2str(i)]);
                dye = dc_roms_zslice_var( ...
                    dc_roms_read_data(runs.dir,'dye_04',i), depth,grdr);
                csdye = dc_roms_zslice_var( ...
                    dc_roms_read_data(runs.dir,'dye_01',i), depth,grdr);
                u = dc_roms_zslice_var( ...
                    dc_roms_read_data(runs.dir,'u',i), depth,grdu);
                v = dc_roms_zslice_var( ...
                    dc_roms_read_data(runs.dir,'v',i), depth,grdv);
                % get on interior RHO points
                u = avg1(u(:,2:end-1),1);
                v = avg1(v(2:end-1,:),2);
                % decimate for quiver
                u = u(1:dxi:end,1:dyi:end);
                v = v(1:dxi:end,1:dyi:end);
            end
            set(hcsdye,'CData',fillnan( ...
                    (csdye/1000 < clim_csd(2)) ...
                    .* csdye/1000,0));
            caxis(clim_csd);colormap(cmcsd); freezeColors; 
            set(heddye,'CData',-addnan(-dye,-0.1));
            caxis(clim_edd);colormap(cmedd);freezeColors;
            runs.update_eddy_contour(he,i);
            runs.update_title(ht,titlestr,i);
            set(hq,'UData',u);
            set(hq,'VData',v);
            runs.video_update();
            pause(0.01);
        end

        runs.video_write();
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
        xax = grd.xax(:,:,1)/1000; yax=  grd.yax(:,:,1)/1000; clear grd;
        tt = 1;
        [~,hc] = contourf(xax,yax,var(:,:,tt));
        hold on
        he = runs.plot_eddy_contour('contour',tind(1) + tt-1);
        shading flat;
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
            set(hc,'ZData',var(:,:,tt));
            shading flat
            runs.update_eddy_contour(he,tind(1) + tt-1);
            set(ht,'String',[varname ' | z = ' num2str(depth) ' m | t = ' ...
            num2str(runs.time(tind(1)+tt-1)/86400) ' days']);
            pause();
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
            if strcmpi(plottype,'contourf') || strcmpi(plottype,'contour')
                eval(['[cc,hplot] = ' plottype '(runs.rgrid.xr/1000,runs.rgrid.yr/1000,'...
                    'runs.zeta(:,:,tt));']);
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
        hold on;
        [~,hplot] = contour(runs.eddy.xr/1000,runs.eddy.yr/1000, ...
                    mask,'Color','k','LineWidth',1);
    end
    function update_eddy_contour(runs,handle,tt)
        try 
            mask = runs.eddy.vormask(:,:,tt);
        catch
            mask = runs.eddy.mask(:,:,tt);
        end
        for ii=1:length(handle)
            set(handle(ii),'ZData',mask);
        end
    end

    function [ht] = set_title(runs,titlestr,tt)
        ht = title([titlestr ' | ' runs.name ' | ' num2str(runs.time(tt)/86400)  ' days']);
    end
    function update_title(runs,ht,titlestr,tt)
        set(ht,'String',[titlestr ' | ' runs.name ' | ' num2str(runs.time(tt)/86400)  ' days']);
    end

    function [hplot] = plot_bathy(runs,plottype,color)
        if ~exist('color','var'), color = 'w'; end
        if strcmpi(plottype,'contour')
            [cc,hplot] = contour(runs.rgrid.xr/1000,runs.rgrid.yr/1000, ...
                            runs.rgrid.h',[200 500 1000 1500 2000],'k');
            clabel(cc,hplot,'LabelSpacing',108*3);
            if runs.bathy.axis == 'y'
                liney(runs.bathy.xsb/1000,'shelfbreak',color);
                liney(runs.bathy.xsl/1000,'slopebreak',color);
            else
                linex(runs.bathy.xsb/1000,'shelfbreak',color);
                liney(runs.bathy.xsl/1000,'slopebreak',color);
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