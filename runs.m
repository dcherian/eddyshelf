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
    % dye names
    csdname; asdname; zdname; eddname;
    % vorticity budget
    vorbudget;
    % grid & bathymetry
    rgrid; bathy
    % float data
    roms; ltrans;
    % eddy track data
    eddy; noeddy;
    % threshold values
    eddy_thresh = 0.7;
    % initial params
    params
    % fluxes - cross-shore & along-shore
    csflux; asflux;
    % transport - TO BE DEPRECATED
    eutrans;
    % streamer properties
    streamer;
    % water mass statistics
    water;
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
            files = roms_find_file(dir,'avg');
            runs.out_file = [runs.dir '/' files{1}];
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
        runs.rgrid.dx = mean(1./runs.rgrid.pm(:));
        runs.rgrid.dy = mean(1./runs.rgrid.pn(:));

        % read zeta
        if ~runs.givenFile
            runs.zeta = dc_roms_read_data(dir,'zeta',[],{},[],runs.rgrid);
            runs.time = dc_roms_read_data(dir,'ocean_time');%,[],{},[],runs.rgrid);
            %try
            %    runs.csdye  = roms_read_data(dir,runs.csdname, ...
            %        [1 1 runs.rgrid.N 1],[Inf Inf 1 Inf]);
            %catch ME
            %end
        else
            runs.zeta = double(ncread(runs.out_file,'zeta'));
            runs.time = double(ncread(runs.out_file,'ocean_time'));
        end

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
                    /mean(runs.rgrid.f(:))/pi;
        runs.rrshelf = sqrt(runs.params.phys.N2)*max(runs.bathy.hsb) ...
                    /mean(runs.rgrid.f(:))/pi;

             % figure out dye names
             for ii=1:4
                 % dye name
                 dname = ['dye_0' num2str(ii)];
                 try % see if variable exists in ini
                     vname = [];
                     % dye description
                     ddesc = ncreadatt([runs.dir roms_find_file(runs.dir,'ini')], ...
                                         dname,'long_name');
                     if strfind(ddesc,'cross shelf'), runs.csdname = dname; end
                     if strfind(ddesc,'z dye'), runs.zdname = dname; end
                     if strfind(ddesc,'along shelf'), runs.asdname = dname; end
                     if strfind(ddesc,'eddy dye'), runs.eddname = dname; end

%                     % see if variable is in output files
                     %try
                     %    runs.(vname) = roms_read_data(filename,dname ...
                     %       ,[1 1 runs.rgrid.N 1],[Inf Inf 1 Inf]);
                     %catch ME
                     %    warning([dname 'not in output files']);
                     %end
                 catch ME
                     warning([dname 'not found in ini file']);
                 end
             end
        try
            runs.roms = floats('roms',runs.flt_file,runs.rgrid);
        catch
        end

        % load eddy track
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

        % scale time by eddy translation
        runs.eddy.tscaleind = find_approx(runs.eddy.my, runs.bathy.xsl, 1);
        runs.eddy.tscale = runs.time(runs.eddy.tscaleind);

        % rerun track_eddy if not new enough
        if ~isfield(runs.eddy,'vor')
            runs.eddy = track_eddy(dir);
        end

        % load streamer data if it exists.
        if exist([dir '/streamer.mat'], 'file') && reset ~= 1
            disp('Loading streamer data');
            streamer = load([dir '/streamer.mat'],'streamer');
            runs.streamer = streamer.streamer;
            clear streamer;
        end

        % load water mass data
        if exist([dir '/watermass.mat'],'file') && reset ~= 1
            disp('Loading water mass data');
            water = load([dir '/watermass.mat'], 'water');
            runs.water = water.water;
            clear water
        end

        % load fluxes if the file exists
        if exist([dir '/fluxes.mat'],'file') && reset ~= 1
            disp('Loading fluxes');
            data = load([dir '/fluxes.mat']);
            runs.csflux = data.csflux;
            runs.asflux = data.asflux;
            clear data
        end

        % extra processing of eddy track
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
            ix = vecfind(runs.eddy.xr(:,1), runs.eddy.mx);
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
        start = [1 1 runs.rgrid.N 1];
        count = [Inf Inf 1 Inf];
        stride = [1 1 1 1];

        if runs.givenFile
            runs.usurf = double(squeeze(ncread(runs.out_file, ....
                'u',start,count,stride)));
        else
            runs.usurf = dc_roms_read_data(runs.dir,'u', ...
                [],{'z' runs.rgrid.N runs.rgrid.N},[],runs.rgrid);
        end
        if runs.givenFile
            runs.vsurf = double(squeeze(ncread(runs.out_file, ....
                'v',start,count,stride)));
        else
            runs.vsurf = dc_roms_read_data(runs.dir,'v', ...
                [],{'z' runs.rgrid.N runs.rgrid.N},[],runs.rgrid);
        end
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
        % need some kind of initial time instant - decided by streamer mask
        % now
        runs.eutrans = [];
        t0 = find(repnan(runs.streamer.time,0) == 0,1,'last') + 1;
        tinf = length(runs.time);
        revind = runs.eddy.trevind;
        h = runs.bathy.h(2:end-1,2:end-1);

        ix = vecfind(runs.eddy.xr(:,1),runs.eddy.mx(t0:end));
        iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.my(t0:end));
        hcen = h(sub2ind(size(runs.eddy.xr),ix,iy))';

        iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.se(t0:end));
        hedge = h(sub2ind(size(runs.eddy.xr),ix,iy))';
        distance = 5*runs.rrshelf; % 5 times rossby radius

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

        % initialize
        runs.eutrans.Itrans = nan([tinf length(loc)]);
        runs.eutrans.nodye.Itrans = nan([tinf length(loc)]);

        % extract streamer mask
        strmask = reshape(full(runs.streamer.west.mask), runs.streamer.sz4dfull);

        % loop over all isobaths
        for kk=1:length(loc)
            % read along-shore section of cross-shore vel.
            % dimensions = (x/y , z , t )
            %cs_vel = double(squeeze(ncread(runs.out_file,csvelid, ...
            %    [1 runs.eutrans.ix(kk) 1 t0],[Inf 1 Inf Inf])));
            cs_vel = dc_roms_read_data(runs.dir, csvelid, ...
                [t0 Inf],{runs.bathy.axis runs.eutrans.ix(kk) runs.eutrans.ix(kk)}, ...
                [],runs.rgrid);
            mask = nan(size(cs_vel));
            for tt=1:size(cs_vel,3)
                mask(1:iwest(tt),:,tt) = 1;
            end
            % restrict calculation to region above shelfbreak depth
            zmask = (abs(squeeze(runs.rgrid.z_r(:,runs.eutrans.ix(kk),:))   )' ...
                            < runs.bathy.hsb);
            mask = bsxfun(@times,mask,fillnan(zmask,0));

            runs.eutrans.nodye.trans(:,:,kk) = squeeze(trapz( ...
                            runs.rgrid.z_r(:,runs.eutrans.ix(kk),1), ...
                            mask .* cs_vel,2));

            runs.eutrans.nodye.Itrans(t0:tinf,kk) = squeeze(nansum( ...
                runs.eutrans.nodye.trans(:,:,kk) ...
                                        .* runs.rgrid.dx,1))';

            % if I have passive tracer info I can calculate transport
            % using that
            mask = nan(size(cs_vel));
            % mark eastern edge as edge of region I'm interested in
            % removes streamer associated with cyclone running away
            ieast = vecfind(runs.eddy.xr(:,1),runs.eddy.ee);
            for tt=1:size(cs_vel,3)
                mask(1:ieast(tt),:,tt) = 1;
            end

            mask = bsxfun(@times,mask,fillnan(zmask,0));
            % dye_01 is always cross-shore dye
            dye = dc_roms_read_data(runs.dir,runs.csdname, ...
                [t0 Inf],{runs.bathy.axis runs.eutrans.ix(kk) runs.eutrans.ix(kk)}, ...
                [],runs.rgrid);
            dyemask = (dye >= runs.bathy.xsb) & ...
                        (dye <=(runs.bathy.xsb + distance));
            mask = mask .* fillnan(dyemask,0);
            runs.eutrans.trans(:,:,kk) = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.eutrans.ix(kk),1), ...
                    mask .* cs_vel,2));
            runs.eutrans.Itrans(t0:tinf,kk) = squeeze(nansum( ...
                        runs.eutrans.trans(:,:,kk) .* dx,1))';

            % all runs now have passive tracer. I use streamer mask to
            % calculate transport
            mask = squeeze(strmask(:,runs.eutrans.ix(kk),:,t0:tinf));

            % (x,t,location)
            runs.streamer.trans(:,:,kk) = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.eutrans.ix(kk),1), ...
                    mask .* cs_vel,2));
            % integrate in x get (t, location)
            runs.streamer.Itrans(t0:tinf,kk) = squeeze(nansum( ...
                        runs.streamer.trans(:,:,kk) .* dx,1))';
        end

        %% plot transport

        figure;
        subplot(6,1,[1 2])
        plot(runs.time/86400,runs.eutrans.Itrans/1e6);
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
        plot(runs.time/86400,runs.eutrans.Itrans/1e6,'-');
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
        plot(runs.time(t0:end)/86400,(runs.eutrans.nodye.Itrans(:,arr) - runs.eutrans.Itrans(:,arr)) ...
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

    function [] = eddyvordiag(runs)

         if isempty(runs.usurf) || isempty(runs.vsurf)
             runs.read_velsurf;
         end

         ux = bsxfun(@rdivide,diff(runs.usurf,1,1),diff(runs.rgrid.x_u',1,1));
         uy = bsxfun(@rdivide,diff(runs.usurf,1,2),diff(runs.rgrid.y_u',1,2));

         vx = bsxfun(@rdivide,diff(runs.vsurf,1,1),diff(runs.rgrid.x_v',1,1));
         vy = bsxfun(@rdivide,diff(runs.vsurf,1,2),diff(runs.rgrid.y_v',1,2));

         ux = ux(:,2:end-1,:);
         vy = vy(2:end-1,:,:);
         vx = avg1(avg1(vx,1),2);
         uy = avg1(avg1(uy,1),2);

         ow = (ux-vy).^2 + (vx+uy).^2 - (vx-uy).^2;


    end

    % plot eddy parameters with time - good for comparison
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

        %% water mass plots
        time = runs.time/86400;
        figure;
        % by regions
        subplot(3,1,1)
        semilogy(time, runs.water.off.deep, ...
                time, runs.water.sl.deep, ...
                time, runs.water.sh.deep, ...
                time, runs.water.edd.deep, ...
                time, runs.water.mix.deep);
        legend('offshore','slope','shelf','eddy','mix');
        ylabel('Deep region (m^3)');
        subplot(312)
        semilogy(time, runs.water.off.slope, ...
                time, runs.water.sl.slope, ...
                time, runs.water.sh.slope, ...
                time, runs.water.edd.slope, ...
                time, runs.water.mix.slope);
        ylabel('Slope region (m^3)');
        subplot(313)
        semilogy(time, runs.water.off.shelf, ...
                time, runs.water.sl.shelf, ...
                time, runs.water.sh.shelf, ...
                time, runs.water.edd.shelf, ...
                time, runs.water.mix.shelf);
        ylabel('Shelf region (m^3)');
        xlabel('Time (days)');

        %% study plumes
        % offshore plume on shelf &
        % shelf plume on slope
        figure;
        subplot(211)
        ax = plotyy( ...%time,runs.water.eddmix.xshelf/1000 - runs.eddy.vor.ee/1000, ...
               time,runs.water.eddmix.zshelf./runs.bathy.hsb, ...
               time,runs.water.eddmix.yshelf/1000 - runs.bathy.xsb/1000);
        set(get(ax(1),'YLabel'),'String','Zcentroid / Depth at shelfbreak');
        set(get(ax(2),'YLabel'),'String','Distance from shelfbreak (km)');
        set(ax(1),'Ylim',[-1 0],'YTick',[-1 -0.5 0]);
        title('Offshore water plume on shelf');
        subplot(212)
        ax = plotyy( ...%time,runs.water.eddmix.xshelf/1000 - runs.eddy.vor.ee/1000, ...
               time,runs.water.sh.zslope./runs.bathy.hsb, ...
               time,runs.water.sh.yslope/1000 - runs.bathy.xsb/1000);
        set(get(ax(1),'YLabel'),'String','Zcentroid / Depth at shelfbreak');
        set(get(ax(2),'YLabel'),'String','Distance from shelfbreak (km)');
        set(ax(1),'Ylim',[-1 0],'YTick',[-1 -0.5 0]);
        title('Shelf water plume on slope');

        % by water masses
%         if ~isempty(runs.water.off.deep)
%             subplot(622)
%             semilogy(time, runs.water.off.deep, ...
%                     time, runs.water.off.slope, ...
%                     time, runs.water.off.shelf);
%             ylabel('offshore water (m^3)');
%             legend('deep region','slope','shelf');
%
%             subplot(624)
%             semilogy(time, runs.water.sl.deep, ...
%                     time, runs.water.sl.slope, ...
%                     time, runs.water.sl.shelf);
%             ylabel('slope');
%
%             subplot(626)
%             semilogy(time, runs.water.sh.deep, ...
%                     time, runs.water.sh.slope, ...
%                     time, runs.water.sh.shelf);
%             ylabel('shelf (m^3)');
%
%             subplot(628)
%             semilogy(time, runs.water.edd.deep, ...
%                     time, runs.water.edd.slope, ...
%                     time, runs.water.edd.shelf);
%             ylabel('eddy(m^3)');
%
%             subplot(6,2,10)
%             semilogy(time, runs.water.mix.deep, ...
%                     time, runs.water.mix.slope, ...
%                     time, runs.water.mix.shelf);
%             ylabel('eddy mix');
%         end
%
        %% eddy upwelling + vertical scale
        if isfield(runs.eddy.vor,'vol')
            figure;
            subplot(211)
            plot(runs.time/86400,runs.eddy.vor.vol);
            ylabel('Volume (m^3)');
            subplot(212)
            plot(runs.time/86400,abs(runs.eddy.vor.zdcen));
            hold on
            plot(runs.time/86400,abs(runs.eddy.vor.zcen),'g');
            plot(eddy.t,eddy.hcen/2,'b');
            plot(eddy.t, ...
                runs.params.phys.f0 / sqrt(runs.params.phys.N2) * runs.eddy.dia,'r');
            plot(eddy.t,eddy.Lz2,'k');
            xlabel('Time (days)');
            ylabel('Z-scale (m)');
            linex(tind);
            suplabel(runs.dir,'t');
            packrows(2,1);
            legend('z-centroid','zdye-centroid', ...
                'H_{center}/2','f*L/N','vertical scale','Location','SouthWest');
        end
    end

    % study along-shore jet
    function [] = jetdetect(runs)
        % rossby radius
        rr = runs.rrshelf;
        % number of rossby radii east of eddy to plot section
        nrr = 8;

        t0 = 55;
        ix0 = vecfind(runs.rgrid.x_u(1,:),runs.eddy.vor.cx(t0:end));
        % along-shore velocity
        if runs.bathy.axis == 'y'

            uas = dc_roms_read_data(runs.dir,'u',[t0 Inf], ...
                {'y' 1 runs.bathy.isl},[],runs.rgrid);
            zas = permute(runs.rgrid.z_u(:,1:runs.bathy.isl,:),[3 2 1]);
        end

        yz = repmat(runs.rgrid.y_u(1:runs.bathy.isl,1),[1 runs.rgrid.N]);

        % first one moves with eddy
        xind = ix0(1) + [nan 10 20] *ceil(rr/runs.rgrid.dx);
        %%
        tt = 1;
        xind(1) = ix0(tt) + nrr * ceil(rr/runs.rgrid.dx);
        subplot(2,3,[1 2 3])
        hzeta = runs.plot_zeta('contourf',t0+tt-1);
        hlines = linex(xind*runs.rgrid.dx/1000);
        colorbar

        subplot(234)
        [~,huas1] = contourf(yz/1000,squeeze(zas(xind(1),:,:)), ...
                        squeeze(uas(xind(1),:,:,tt)));
        colorbar; caxis([-1 1]*0.1); ylim([-300 0]);


        subplot(235)
        [~,huas2] = contourf(yz/1000,squeeze(zas(xind(2),:,:)), ...
                        squeeze(uas(xind(2),:,:,tt)));
        colorbar; caxis([-1 1]*0.1); ylim([-1000 0]);


        subplot(236)
        [~,huas3] = contourf(yz/1000,squeeze(zas(xind(3),:,:)), ...
                        squeeze(uas(xind(3),:,:,tt)));
        colorbar; caxis([-1 1]*0.1); ylim([-1000 0]);

        for tt=2:size(uas,4)
            runs.update_zeta(hzeta,t0+tt-1);

            xind(1) = ix0(tt) + nrr * ceil(rr/runs.rgrid.dx);
            set(hlines(1),'XData',[1 1]*xind(1)*runs.rgrid.dx/1000);

            set(huas1,'ZData',squeeze(uas(xind(1),:,:,tt)));
            set(huas2,'ZData',squeeze(uas(xind(2),:,:,tt)));
            set(huas3,'ZData',squeeze(uas(xind(3),:,:,tt)));
            pause(0.2);
        end
    end

    % study vorticity export onto shelf
    function [] = vorexport(runs)
        vorname = [runs.dir '/ocean_vor.nc'];
        t0 = 1;

        start = [1 1 1 t0];
        % isb-1 since they are on interior RHO points
        count = [Inf runs.bathy.isb-1 Inf Inf];

        Las = max(runs.rgrid.x_rho(:));
        Lcs = runs.bathy.xsb;

        dx = runs.rgrid.xr(2,1)-runs.rgrid.xr(1,1);

        % both at interior RHO points
        pv = ncread(vorname,'pv',start,count);
        rv = avg1(avg1(ncread(vorname,'rv',start,count+[0 1 0 0]),1),2);

        if runs.bathy.axis == 'y'
            csvel = avg1(avg1(dc_roms_read_data(runs.dir,'v',[t0 Inf], ...
                {'y' runs.bathy.isb-1 runs.bathy.isb},[],runs.rgrid),2),3);
            csvel = csvel(2:end-1,:,:,:);

            % location for calculation along-shore flux of rv/pv
            asloc = runs.eddy.vor.ee(t0:end) + dx*3;

            % full RHO points
            iasmin = find(runs.rgrid.xr(:,1) == min(asloc));
            iasmax = find(runs.rgrid.xr(:,1) == max(asloc));

            % faster to read whole thing and then discard
            asvel = avg1(avg1(dc_roms_read_data(runs.dir,'u',[t0 Inf], ...
                {'y' 2 runs.bathy.isb}, [], runs.rgrid),1),3);
            % average to RHO points
            asvel = avg1(asvel(iasmin-1:iasmax,:,:,:),1);
            % shift indices since I'm now on interior RHO points for pv,rv
            iasmin = iasmin - 1;
            iasmax = iasmax - 1;

            % use center because export occurs west of the eastern edge
            csmask = bsxfun(@gt, runs.eddy.xr(:,1),  ...
                permute(runs.eddy.vor.cx(t0:end),[3 1 2]));
        else
            error('Not implemented for NS isobaths');
            csvel = avg1(avg1(dc_roms_read_data(runs.dir,'u',[t0 Inf], ...
                {'x' runs.bathy.isb runs.bathy.isb+1},[],runs.rgrid),1),3);
            csvel = csvel(2:end-1,:,:,:);
            % use center because export occurs north of the southern edge
            csmask = bsxfun(@gt, runs.eddy.yr(1,:)',  ...
                permute(runs.eddy.vor.cy(t0:end),[3 1 2]));
        end
        %%
        % csvel = csvel(:,:,:,1:38)
        % These are fluxes across the shelfbreak  - technically represent
        % eddy only. I need to quantify vorticity exported permanently onto
        % the shelf , so it might better to do an along shore flux
        % downstream (east) of the eddy, over the shelf. The idea is that
        % vorticity is dumped on the shelf and then moves downstream. So
        % along-shore flux when calculated sufficiently far downstream of
        % the eddy, should represent permanent export of vorticity (and
        % mass) on the shelf

        % isb - 1 to account for being on interior RHO points - isb
        % includes the boundary points too
        pvcsflux = squeeze(bsxfun(@times,pv(:,runs.bathy.isb-1,:,:), csvel));
        rvcsflux = squeeze(bsxfun(@times,rv(:,runs.bathy.isb-1,:,:), csvel));

        asmask = bsxfun(@eq, runs.rgrid.xr(iasmin+1:iasmax+1,2:runs.bathy.isb), ...
            permute(asloc,[1 3 4 2]));

        pvasflux = squeeze(sum( ...
            bsxfun(@times, pv(iasmin:iasmax,:,:,:).*asvel, asmask),1));
        rvasflux = squeeze(sum( ...
            bsxfun(@times, rv(iasmin:iasmax,:,:,:).*asvel, asmask),1));

        dV = avg1(runs.rgrid.dV,3);
        dVas = squeeze(sum( ...
             bsxfun(@times, dV(iasmin:iasmax,2:runs.bathy.isb,:), asmask)));

        runs.csflux.pv = squeeze(sum(sum(bsxfun(@times, ...
                     bsxfun(@times,pvcsflux, squeeze( ...
                     dV(2:end-1,runs.bathy.isb-1,:)))...
                     ,csmask),1),2)) ...
                    /runs.bathy.hsb/Las;

        runs.csflux.rv = squeeze(sum(sum(bsxfun(@times, ...
                     bsxfun(@times,rvcsflux, squeeze( ...
                     dV(2:end-1,runs.bathy.isb-1,:)))...
                     ,csmask),1),2)) ...
                    /runs.bathy.hsb/Las;

        PVASFLUX = squeeze(sum(sum(pvasflux .* dVas,1),2)) ...
                    /runs.bathy.hsb/Lcs;
        RVASFLUX = squeeze(sum(sum(rvasflux .* dVas,1),2)) ...
                    /runs.bathy.hsb/Lcs;

        oPVCSFLUX = orderofmagn(runs.csflux.pv);
        oRVCSFLUX = orderofmagn(runs.csflux.rv);
        oPVASFLUX = orderofmagn(runs.asflux.pv);
        oRVASFLUX = orderofmagn(runs.asflux.rv);

        oPV = min(oPVCSFLUX, oPVASFLUX);
        oRV = min(oRVCSFLUX, oRVASFLUX);

        %%
%        figure;
%        tt0 = 40;
%        isb = runs.bathy.isb;
%        [~,hc] = contourf(runs.rgrid.xr(2:end-1,2:isb)/1000, ...
%                          runs.rgrid.yr(2:end-1,2:isb)/1000, ...
%                          rv(:,:,end,tt0));
%        colorbar;
%         caxis([-1 1]*4*10^(orderofmagn(rv(:,:,end,:))));
%
%         hlines = linex([asloc(1) asloc(1)]/1000);
%         for tt=tt0+1:size(rv,4)
%             set(hc,'ZData',rv(:,:,end,tt));
%
%             set(hlines(1),'XData',[1 1]*runs.eddy.cx(tt)/1000);
%             set(hlines(2),'XData',[1 1]*asloc(tt-tt0+1)/1000);
%             pause();
%         end
%         % quantify loss of vorticity in eddy
%         xmin = min(runs.eddy.vor.we); xmax = max(runs.eddy.vor.ee);
%         ymin = min(runs.eddy.vor.se); ymax = max(runs.eddy.vor.ne);
%
%         ixmin = find_approx(runs.eddy.xr(:,1),xmin);
%         ixmax = find_approx(runs.eddy.xr(:,1),xmax);
%         iymin = find_approx(runs.eddy.yr(1,:),ymin);
%         iymax = find_approx(runs.eddy.yr(1,:),ymax);
%
%         tic;
%         disp('Reading pv and rv for eddy');
%         pveddy = ncread(vorname,'pv',[ixmin+1 iymin+1 1 t0],[ixmax+1 iymax+1 Inf Inf]);
%         rveddy = ncread(vorname,'pv',[ixmin+1 iymin+1 1 t0],[ixmax+1 iymax+1 Inf Inf]);
%         toc;
%
%         % TODO: add dye mask here
%         bsxfun(@times, bsxfun(@times,pveddy, ...
%             permute(runs.eddy.vormask(ixmin:ixmax,iymin:iymax,:),[1 2 4 3])), ...
%             runs.rgrid.dV(ixmin:ixmax,iymin:iymax,:));
        %%
        figure;
        subplot(211)
        plot(runs.time(t0:end)/86400, runs.csflux.pv./10^(oPVCSFLUX), ...
              runs.time(t0:end)/86400, runs.asflux.pv./10^(oPVASFLUX));
         legend(['CS flux x 10^{' num2str(oPVCSFLUX) '}'], ...
                ['AS flux x 10^{' num2str(oPVASFLUX) '}']);
%        legend('cross-shore','along-shore');
        ylabel(['PV flux']);
        set(gca,'YTick',[-5:5]);
        beautify([18 16 16]);
        liney(0);
        title([runs.name ' | CS =  across shelfbreak | AS = east edge + 3dx']);
        subplot(212)
        plot(runs.time(t0:end)/86400, runs.csflux.rv./10^(oRVCSFLUX), ...
             runs.time(t0:end)/86400, runs.asflux.rv./10^(oRVASFLUX));

            set(gca,'YTick',[-5:5]);
         legend(['CS flux x 10^{' num2str(oRVCSFLUX) '}'], ...
                ['AS flux x 10^{' num2str(oRVASFLUX) '}']);
%        legend('cross-shore','along-shore');
        liney(0); ylim([-3 3]);
        beautify([18 16 16]);
        ylabel(['Rel. Vor. Flux']);
        xlabel('Time (days)');

        %% save fluxes
        csflux = runs.csflux;
        asflux = runs.asflux;
        save([runs.dir '/fluxes.mat'], 'csflux', 'asflux');

    end

    % quantify cross-shelfbreak and along-shelfbreak fluxes of a whole
    % bunch of stuff:
    function [] = fluxes(runs)

        % Things I want to calculate fluxes of:
        % 1. RV & PV
        % 2. Shelf water dye & eddy dye
        ticstart = tic;

        vorname = [runs.dir '/ocean_vor.nc'];
        % need some kind of initial time instant - decided by streamer mask
        % now
        runs.csflux = [];
        runs.asflux = [];
        t0 = 1;%find(repnan(runs.streamer.time,0) == 0,1,'last') + 1;
        tinf = length(runs.time);
        revind = runs.eddy.trevind;

        % quantities for area-averaging if needed
        Las = max(runs.rgrid.x_rho(:));
        %Lcs = runs.bathy.xsb;
        h = runs.bathy.h(2:end-1,2:end-1);

        % sort out isobaths across which to calculate transports
        if runs.params.bathy.axis == 'x'
            csvelid = 'u';
            asvelid = 'v';
            bathyax = 1;
            error(' not built for north-south isobaths');
            indices = vecfind(runs.eddy.xr(:,1),loc);
        else
            csvelid = 'v';
            asvelid = 'u';
            bathyax = 2;
            loc = sort([runs.bathy.xsb  runs.bathy.xsl...
                    nanmean(runs.eddy.se(revind:end)) ...
                    nanmean(runs.eddy.cy(revind:end))]);
            indices = vecfind(runs.eddy.yr(1,:),loc);
                %runs.rgrid.y_rho(vecfind(runs.bathy.h(1,:),[250 1000]),1)']);
        end

        % save locations
        runs.csflux.x = loc;
        % save indices for locations - w.r.t INTERIOR RHO POINTS
        runs.csflux.ix = indices;
        % save isobath values
        runs.csflux.h = ceil(runs.bathy.h(1,runs.csflux.ix));

        runs.csflux.comment = ['shelf / slope / eddy= which water mass am I ', ...
                            ' targeting? |\n (x,ix,h) = (location, index, depth) ' ...
                            'at which I''m calculating transport |\n '];

        % initialize
        runs.csflux.west.shelf = nan([tinf length(loc)]);
        runs.csflux.west.slope = nan([tinf length(loc)]);
        runs.csflux.west.eddy = nan([tinf length(loc)]);
        runs.csflux.west.pv = nan([tinf length(loc)]);
        runs.csflux.west.rv = nan([tinf length(loc)]);

        runs.csflux.east.shelf = nan([tinf length(loc)]);
        runs.csflux.east.slope = nan([tinf length(loc)]);
        runs.csflux.east.eddy = nan([tinf length(loc)]);
        runs.csflux.east.pv = nan([tinf length(loc)]);
        runs.csflux.east.rv = nan([tinf length(loc)]);

        % east and west (w.r.t eddy center) masks
        % use center because export occurs west of the eastern edge
        westmask = bsxfun(@lt, runs.eddy.xr(:,1),  ...
                    runs.eddy.vor.cx(t0:end));
        eastmask = bsxfun(@gt, runs.eddy.xr(:,1), ...
                    runs.eddy.vor.cx(t0:end));

        % loop over all isobaths
        for kk=1:length(loc)
            % read along-shore section of cross-shore vel.
            % dimensions = (x/y , z , t )
            % average carefully to get values at RHO point
            csvel = avg1(dc_roms_read_data(runs.dir, csvelid, ...
                [t0 Inf],{runs.bathy.axis runs.csflux.ix(kk)-1 runs.csflux.ix(kk)}, ...
                [],runs.rgrid),bathyax);
            csvel = csvel(2:end-1,:,:,:);
            % process cross-shelf dye
            csdye = dc_roms_read_data(runs.dir,runs.csdname, ...
                [t0 Inf],{runs.bathy.axis runs.csflux.ix(kk)+1 runs.csflux.ix(kk)+1}, ...
                [],runs.rgrid);
            csdye = permute(csdye(2:end-1,:,:), [1 4 2 3]);

            % read eddye
            eddye = dc_roms_read_data(runs.dir, runs.eddname, ...
                [t0 Inf],{runs.bathy.axis runs.csflux.ix(kk)+1 runs.csflux.ix(kk)+1}, ...
                [],runs.rgrid);
            eddye = permute(eddye(2:end-1,:,:), [1 4 2 3]);

            % define water masses
            shelfmask = (csdye <= runs.bathy.xsb);
            slopemask = (csdye > runs.bathy.xsb) & ...
                        (csdye <=runs.bathy.xsl);
            eddymask = eddye > runs.eddy_thresh;

            % transports
            runs.csflux.shelf(:,:,kk) = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1), ...
                    shelfmask .* csvel,3));
            runs.csflux.slope(:,:,kk) = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1), ...
                    slopemask .* csvel,3));
            runs.csflux.eddy(:,:,kk) = squeeze(trapz( ...
                    runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1), ...
                    eddymask .* csvel,3));

            % west of center
            runs.csflux.west.shelf(t0:tinf,kk) = squeeze(nansum( ...
                        bsxfun(@times, runs.csflux.shelf(:,:,kk) .* westmask, ...
                        1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.west.slope(t0:tinf,kk) = squeeze(nansum( ...
                        bsxfun(@times, runs.csflux.slope(:,:,kk) .* westmask, ...
                        1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.west.eddy(t0:tinf,kk) = squeeze(nansum( ...
                        bsxfun(@times, runs.csflux.eddy(:,:,kk) .* westmask, ...
                        1./runs.rgrid.pm(1,2:end-1)'),1))';

            % east of center
            runs.csflux.east.shelf(t0:tinf,kk) = squeeze(nansum( ...
                        bsxfun(@times, runs.csflux.shelf(:,:,kk) .* eastmask, ...
                        1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.east.slope(t0:tinf,kk) = squeeze(nansum( ...
                        bsxfun(@times, runs.csflux.slope(:,:,kk) .* eastmask, ...
                        1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.east.eddy(t0:tinf,kk) = squeeze(nansum( ...
                        bsxfun(@times, runs.csflux.eddy(:,:,kk) .* eastmask, ...
                        1./runs.rgrid.pm(1,2:end-1)'),1))';

            % process pv
            % both at interior RHO points
            start = [1 runs.csflux.ix(kk) 1 t0];
            count = [Inf 1 Inf Inf];

            pv = ncread(vorname,'pv',start,count);
            rv = avg1(avg1(ncread(vorname,'rv',start,count+[0 1 0 0]),1),2);

            % get vorticity fluxes
            % first, depth integrated
            pvcsflux = squeeze(trapz(avg1(runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1),1), ...
                        pv .* avg1(csvel,3),3));
            rvcsflux = squeeze(trapz(avg1(runs.rgrid.z_r(:,runs.csflux.ix(kk)+1,1),1), ...
                        rv .* avg1(csvel,3),3));

            % now east & west of eddy center
            runs.csflux.west.pv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, pvcsflux .* westmask, ...
                1./runs.rgrid.pm(1,2:end-1)'),1))';
            runs.csflux.east.pv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, pvcsflux .* eastmask, ...
                1./runs.rgrid.pm(1,2:end-1)'),1))';

            runs.csflux.west.rv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, rvcsflux .* westmask, ...
                1./runs.rgrid.pm(1,2:end-1)'),1))';
            runs.csflux.east.rv(t0:tinf,kk) = squeeze(nansum( ...
                bsxfun(@times, rvcsflux .* eastmask, ...
                1./runs.rgrid.pm(1,2:end-1)'),1))';
        end
        toc(ticstart);
        % save fluxes
        csflux = runs.csflux;
        asflux = runs.asflux;
        save([runs.dir '/fluxes.mat'], 'csflux', 'asflux');
    end

    function [] = compare_plot(runs,num)
        eddy = runs.eddy;
        % 86400 since eddy.t is already in days
        eddy.t = eddy.t./ (eddy.tscale/86400);
        ii = num;
        %       colors = distinguishable_colors(10);
        colors = cbrewer('qual', 'Dark2', 8);
        aa = 6; bb = aa*2;
        tloc = eddy.tscaleind;;

        % plot eddy tracks
        figure(1)
        hold on
        %subplot(aa,2,bb);
        %subplot(aa,2,1); hold all
        %pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);            colorbar
        xlabel('X (km)'); ylabel('Y (km)');
        he = plot((eddy.mx-eddy.mx(1))/1000, ...
             (eddy.my-eddy.my(1))/1000,'Color',colors(ii,:),'LineWidth',2);
        hold all;
        addlegend(he, runs.name, 'NorthEast');
        try
            plot((eddy.mx(tloc)-eddy.mx(1))/1000,...
                (eddy.my(tloc)-eddy.my(1))/1000,'*', ...
                 'MarkerSize',12,'Color',colors(ii,:),'LineWidth',2);
        catch ME
        end
            %if runs.bathy.axis == 'x'
        %    plot(eddy.we/1000,eddy.cy/1000,'Color',colors(ii,:),'LineStyle','--');
        %else
        %    plot(eddy.cx/1000,eddy.se/1000,'Color',colors(ii,:),'LineStyle','--');
        %end
        %axis image; axis tight

        % plot eddy properties
        figure(2)
        hold on
        subplot(aa,2,2); hold on
        limx = [0 max([xlim eddy.t])];
        plot(eddy.t,eddy.vor.amp./eddy.amp(1),'Color',colors(ii,:));
        ylabel('amp/amp(t=0) ');xlim(limx);

        if isfield(eddy, 'vol')
            subplot(aa,2,1); hold on
            plot(eddy.t, eddy.vol./eddy.vol(1),'Color', colors(ii,:));
            ylabel('Volume');xlim(limx);

            subplot(aa,2,3); hold on
            plot(eddy.t, eddy.PV./abs(eddy.PV(1)), 'Color', colors(ii,:));
            ylabel('PV/|PV0|');xlim(limx);

            subplot(aa,2,5); hold on
            plot(eddy.t, eddy.RV./abs(eddy.RV(1)), 'Color', colors(ii,:));
            ylabel('RV/|RV0|');xlim(limx);

            subplot(aa,2,7); hold on;
            plot(eddy.t, eddy.KE./eddy.KE(1), 'Color', colors(ii,:));
            ylabel('KE');xlim(limx);

            subplot(aa,2,9); hold on;
            plot(eddy.t, eddy.PE./eddy.PE(1), 'Color', colors(ii,:));
            ylabel('PE');xlim(limx);
        end

        subplot(aa,2,4); hold on
        plot(eddy.t, eddy.Ls/runs.rrdeep,'Color',colors(ii,:));
        ylabel('Ls/RRdeep');xlim(limx);

        subplot(aa,2,6); hold on
        plot(eddy.t,eddy.mvx,'Color',colors(ii,:));
        ylabel('cvx(km/day)');
        ylim([-5 5]);
        liney(0); xlim(limx);
        %plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
        %ylabel('x - center (km)');

        subplot(aa,2,8); hold on
        plot(eddy.t,eddy.mvy,'Color',colors(ii,:));
        ylabel('cvy (km/day)');xlim(limx);
        ylim([-5 5]);
        %plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
        %ylabel('y - center (km)');

        subplot(aa,2,10); hold on
        plot(eddy.t,eddy.Lgauss./max(eddy.Lgauss(1)),'Color',colors(ii,:));
        ylabel('H_{eddy}/H_{eddy0}');xlim(limx);
        %xlabel('time (days)');

        subplot(aa,2,12); hold on
        plot(eddy.t,eddy.prox/1000,'Color',colors(ii,:));
        xlabel('time (days)');
        ylabel('Proximity (km)');xlim(limx);

        subplot(aa,2,11); hold on
        hp = plot(eddy.t,eddy.hcen./max(runs.bathy.h(:)),'Color',colors(ii,:));
        addlegend(hp,runs.name,'SouthWest');
        %        plot(eddy.t,runs.params.phys.f0 / sqrt(runs.params.phys.N2) ...
        %             *  runs.eddy.dia,'Color',colors(ii,:),'LineStyle','--');
        %        legend('H_{center}','f/N*dia');
        xlabel('Time / Time at which center reaches slopebreak');
        ylabel('H_{center}(m)/H_{max}');
        xlim(limx);

        %% plot fluxes

        if isfield(runs.csflux,'west')
            figure(4);
            subplot(4,1,1);
            hold on;
            plot(eddy.t, runs.csflux.west.shelf(:,1), 'Color', colors(ii,:));
            ylabel('Shelf water flux - sb');
            title('West');
            xlim(limx);

            subplot(4,1,2);
            hold on;
            plot(eddy.t, runs.csflux.west.slope(:,1), 'Color', colors(ii,:));
            ylabel('Slope water flux - sb');
            xlim(limx);

            subplot(4,1,3);
            hold on;
            plot(eddy.t, runs.csflux.west.pv(:,1), 'Color', colors(ii,:));
            ylabel('PV flux');
            xlim(limx);

            subplot(4,1,4);
            hold on;
            plot(eddy.t, runs.csflux.west.rv(:,1), 'Color', colors(ii,:));
            ylabel('RV flux');
            xlim(limx);

            figure(5);
            subplot(4,1,1);
            hold on;
            plot(eddy.t, runs.csflux.east.shelf(:,1), 'Color', colors(ii,:));
            ylabel('Shelf water flux - sb');
            title('East');
            xlim(limx);

            subplot(4,1,2);
            hold on;
            plot(eddy.t, runs.csflux.east.slope(:,1), 'Color', colors(ii,:));
            ylabel('Slope water flux - sb');
            xlim(limx);

            subplot(4,1,3);
            hold on;
            plot(eddy.t, runs.csflux.east.pv(:,1), 'Color', colors(ii,:));
            ylabel('PV flux');
            xlim(limx);

            subplot(4,1,4);
            hold on;
            plot(eddy.t, runs.csflux.east.rv(:,1), 'Color', colors(ii,:));
            ylabel('RV flux');
            xlim(limx);
        end

        %% plot water masses
        if isfield(runs.water.off, 'deep')
            linestyle = {'-','--','-.','-'};
            markers = {'none','none','none','.'};
            time = eddy.t;
            % normalize volumes by initial eddy volume
            evol0 = runs.eddy.vol(runs.eddy.tscaleind);

            figure(3);
            set(gcf, 'Renderer', 'painters');
            % by regions
            % colors: off = r, sl = g, sh = b , edd = k, mix = m
            subplot(3,1,1)
            hold on;
            hw = plot(time, runs.water.sl.deep/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{1}, 'Marker', ...
                 markers{1});
            plot(time, runs.water.sh.deep/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{2}, 'Marker', ...
                 markers{2});
            ylabel('Deep region ');
            title('All volumes normalized by eddy volume at t=tscale');
            xlim(limx);
            addlegend(hw, runs.name, 'NorthWest');

            subplot(312)
            hold on;
            plot(time, runs.water.off.slope/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{1}, 'Marker', ...
                 markers{1});
            plot(time, runs.water.sh.slope/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{3}, 'Marker', ...
                 markers{3});
            plot(time, runs.water.edd.slope/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{4}, 'Marker', ...
                 markers{4});
            %plot(time, runs.water.mix.slope/evol0, ['m' linestyle{num}]);
            legend('Offshore','Shelf','Eddy');
            ylabel('Slope region');
            xlim(limx);

            subplot(313)
            hold on;
            plot(time, runs.water.off.shelf/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{1}, 'Marker', ...
                 markers{1});
            plot(time, runs.water.sl.shelf/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{2}, 'Marker', ...
                 markers{2});
            plot(time, runs.water.edd.shelf/evol0, 'Color', ...
                 colors(ii,:), 'LineStyle', linestyle{4}, 'Marker', ...
                 markers{4});
            %plot(time, runs.water.mix.shelf/evol0,['m' linestyle{num}]);
            legend('Offshore','Slope','Eddy');
            ylabel('Shelf region');
            xlabel('Time (days)');
            xlim(limx);
        end
    end

    % this is incomplete
    function [] = tracer_budget(runs)
        tracer = roms_read_data(runs.out_file,runs.zdname);
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

    % calculate upwelling in eddy
    function [] = eddy_upwelling(runs)

        % two possibilities here - use eddye as a
        % 1 - use mask
        % 2 - use vormask
        xr = runs.eddy.xr;
        yr = runs.eddy.yr;

        use_sshmask = 0;

        ixmin = vecfind(runs.rgrid.xr(:,1),runs.eddy.vor.we);
        ixmax = vecfind(runs.rgrid.xr(:,1),runs.eddy.vor.ee);
        iymin = vecfind(runs.rgrid.yr(1,:),runs.eddy.vor.se);
        iymax = vecfind(runs.rgrid.yr(1,:),runs.eddy.vor.ne);

        ixm = min(ixmin); ixM = max(ixmax);
        iym = min(iymin); iyM = max(iymax);

        volume = {'x' ixm ixM; 'y' iym iyM};
        zdye = dc_roms_read_data(runs.dir, runs.zdname, [], volume, [], runs.rgrid);

        % from animate_pt output it looks like vormask tracks the edge of
        % eddye contour pretty well. I don't get the stuff that spreads
        % along shelf but get the eddy pretty well.
        try
            eddye = dc_roms_read_data(runs.dir,runs.eddname, [], volume, [], runs.rgrid);
        catch ME
            warning('no eddy dye (dye_04) found');
            return;
        end

        sz4dfull = size(zdye);
        sz4dsp = [prod(sz4dfull(1:3)) sz4dfull(end)];
        sz3dsp = [sz4dsp(1) 1];

        % make my mask matrices 4d and sparse
        eddye = sparse(reshape(eddye > runs.eddy_thresh, sz4dsp));
        vormask = sparse(reshape(repmat( ...
            permute(runs.eddy.vormask(ixm-1:ixM-1,iym-1:iyM-1,:),[1 2 4 3]) ...
            , [1 1 runs.rgrid.N 1]), sz4dsp));

        dV = reshape(runs.rgrid.dV(ixm:ixM,iym:iyM,:),sz3dsp);

        zdye = reshape(zdye,sz4dsp);

        % first with vormask
        zdyevor = zdye .* vormask;
        zedd = bsxfun(@times, vormask, reshape(permute( ...
            runs.rgrid.z_r(:,iym:iyM,ixm:ixM),[3 2 1]), sz3dsp));
        runs.eddy.vor.vol = full(squeeze(sum(bsxfun(@times,vormask,dV),1))');
        runs.eddy.vor.zdcen = runs.domain_integratesp(zdyevor, dV)' ...
                        ./ runs.eddy.vor.vol;
        runs.eddy.vor.zcen = runs.domain_integratesp(zedd,dV)' ...
                        ./ runs.eddy.vor.vol;

        % then with ssh mask - though really vormask is what i'm looking
        % for
        if use_sshmask
            mask = sparse(reshape(repmat( ...
               permute(runs.eddy.mask(ixm-1:ixM-1,iym-1:iyM-1,:),[1 2 4 3]) ...
               , [1 1 runs.rgrid.N 1]), sz4dsp));
            zdyessh = zdye .* mask;
            runs.eddy.vol = full(squeeze(sum(bsxfun(@times,mask,dV) ,1))');
            runs.eddy.zdcen = runs.domain_integratesp(zdyessh, reshape(dV,sz3dsp))' ...
                           ./ runs.eddy.vol;
        end
    end

    % detect streamer contours
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
        %ix = repmat([1:size(xr,1)]',[1 yend]);
        %iy = repmat([1:yend],[size(xr,1) 1]);

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

            csdye = dc_roms_read_data(runs.dir, runs.csdname, tindices, ...
                        {'y' 1 yend},[],runs.rgrid)/1000;
            zdye  = dc_roms_read_data(runs.dir, runs.zdname, tindices, ...
                        {'y' 1 yend},[],runs.rgrid);
            eddye = dc_roms_read_data(runs.dir,runs.eddname, tindices, ...
                        {'y' 1 yend},[],runs.rgrid);
            %asdye = dc_roms_read_data(runs.dir, runs.asdname, tindices, ...
            %            {'y' 1 yend});

            % identify streamer with 4D data
            % preliminary detection
            % I use cross-shore label to roughly filter first

            % eliminate this step somehow and remove the temporary array?
            streamer1 = (csdye > xsb-10) & (csdye < xsb+30) & (eddye < 0.2);

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

            warning('DO I NEED TO ACCOUNT FOR TILTING IN VERTICAL?');
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

            stream = streamer1(:,:,runs.rgrid.N,:);
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
            bins = -1*[0:dbin:1000];
            % required so that 0 bin doesn't get a ton of points
            %zsf = fillnan(full(zs),0);
            %sz = size(runs.streamer.west.mask(:,tstart:tend));
            parfor kk=1:length(bins)-1
                temparray(kk,:) = sum(bsxfun(@times, (zs < bins(kk) & zs >= bins(kk+1)), ...
                               dVs),1);
            end
            runs.streamer.west.Vbin(:,tstart:tend) = temparray;
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
                    xstr = fliplr(xstr)';
                    ystr = fliplr(ystr)';
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

                xstr = xr(ixstr,1);
                ystr = yr(1,iystr)';
            end

            % fix the starting!!!
            dstr = [0; cumsum(hypot(diff(xstr),diff(ystr)))];

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
                [' contour = 1 in streamer, 0 outside |\n ' ...
                ' (xstr,ystr) = cross-section through streamer (cell array) |\n ' ...
                ' (ixstr,iystr) = indices corresponding to (xstr,ystr) ' ...
                ' - (cell array) |\n dstr = along-streamer distance (cell array)'];
        end

        % save to file
        disp('Writing to file');tic;
        streamer = runs.streamer;
        save([runs.dir '/streamer.mat'],'streamer');
        toc;
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
        zcenbin = vecfind(bins, cut_nan(runs.streamer.west.zcen));
        Vcenbin = diag(runs.streamer.west.Vbin(zcenbin,:));
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

    % offshore water census
    function [] = water_census(runs)

        ticstart = tic;
        % if dye_04 > thresh then it is "eddy water"
        % else i classify it as "mixed"
        eddye_thresh = runs.eddy_thresh;

        % check classified vol against total volume
        debug = 1;

        % my region boundaries are based on location of shelfbreak and
        % slopebreak. Let's make it easy.
        xsl = runs.bathy.xsl;
        isl = runs.bathy.isl;
        xsb = runs.bathy.xsb;
        isb = runs.bathy.isb;

        slab = 20; % read 10 at a time

        sz4dfull = [fliplr(size(runs.rgrid.z_r)) slab];
        sz4dsp = [prod(sz4dfull(1:3)) slab];
        sz3dfull = sz4dfull(1:3);
        sz3dsp = [sz4dsp(1) 1];

        % define cross-shore grid co-ordinate
        if runs.bathy.axis == 'y'
            cs = repmat(runs.rgrid.yr,[1 1 runs.rgrid.N]);
        else
            cs = repmat(runs.rgrid.xr,[1 1 runs.rgrid.N]);
        end

        % define regions
        % deep region
        regdp = sparse(reshape(cs > xsl, sz3dsp));
        regsl = sparse(reshape(cs <= xsl & cs > xsb, sz3dsp));
        regsh = sparse(reshape(cs <=xsb, sz3dsp));

        dV = reshape(runs.rgrid.dV, sz3dsp);

        dz = dV ./ runs.rgrid.dx ./ runs.rgrid.dy;
        xsp = reshape(repmat(runs.rgrid.xr,[1 1 runs.rgrid.N]), sz3dsp);
        ysp = reshape(repmat(runs.rgrid.yr,[1 1 runs.rgrid.N]), sz3dsp);
        zsp = reshape(permute(runs.rgrid.z_r,[3 2 1]), sz3dsp);

        ntime = length(runs.time);

        for tt=1:slab:length(runs.time)
            tend = tt + slab -1;
            if tend > length(runs.time)
                tend = length(runs.time);
                sz4dfull(end) = tend-tt+1;
                sz4dsp(end) = tend-tt+1;
            end
            csdye = dc_roms_read_data(runs.dir,runs.csdname,[tt tend],{},[],runs.rgrid);
            eddye = dc_roms_read_data(runs.dir,runs.eddname,[tt tend],{},[],runs.rgrid);
            % define water masses
            % offshore water
            maskoff = sparse(reshape(csdye > xsl, sz4dsp));
            % slope water
            masksl  = sparse(reshape(csdye <= xsl & csdye > xsb, sz4dsp));
            % shelf water
            masksh  = sparse(reshape(csdye <= xsb, sz4dsp));
            % eddy water
            masked  = sparse(reshape(eddye > eddye_thresh, sz4dsp));
            % "mixed water"
            maskmix = sparse(reshape(eddye <= eddye_thresh & eddye > 0.01,sz4dsp));

            % the eddy's velocity field mixes up csdye and makes it look
            % like slope water?
            % in any case i want all 5 to add up to total volume, so let's
            % remove the volume that's in the eddy from csdye.
            maskoff = maskoff & ~(masked | maskmix);
            masksl  = masksl  & ~(masked | maskmix);
            masksh  = masksh  & ~(masked | maskmix);

            % now census in each region
            % first deep water region
            runs.water.off.deep(tt:tend) = full(sum( ...
                bsxfun(@times, maskoff, regdp .* dV),1));
            runs.water.sl.deep(tt:tend) = full(sum( ...
                bsxfun(@times, masksl, regdp .* dV),1));
            runs.water.sh.deep(tt:tend) = full(sum( ...
                bsxfun(@times, masksh, regdp .* dV),1));
            runs.water.edd.deep(tt:tend) = full(sum( ...
                bsxfun(@times, masked, regdp .* dV),1));
            runs.water.mix.deep(tt:tend) = full(sum( ...
                bsxfun(@times, maskmix, regdp .* dV),1));

            % now slope region
            runs.water.off.slope(tt:tend) = full(sum( ...
                bsxfun(@times, maskoff, regsl .* dV),1));
            runs.water.sl.slope(tt:tend) = full(sum( ...
                bsxfun(@times, masksl, regsl .* dV),1));
            runs.water.sh.slope(tt:tend) = full(sum( ...
                bsxfun(@times, masksh, regsl .* dV),1));
            runs.water.edd.slope(tt:tend) = full(sum( ...
                bsxfun(@times, masked, regsl .* dV),1));
            runs.water.mix.slope(tt:tend) = full(sum( ...
                bsxfun(@times, maskmix, regsl .* dV),1));

            % now shelf region
            runs.water.off.shelf(tt:tend) = full(sum( ...
                bsxfun(@times, maskoff, regsh .* dV),1));
            runs.water.sl.shelf(tt:tend) = full(sum( ...
                bsxfun(@times, masksl, regsh .* dV),1));
            runs.water.sh.shelf(tt:tend) = full(sum( ...
                bsxfun(@times, masksh, regsh .* dV),1));
            runs.water.edd.shelf(tt:tend) = full(sum( ...
                bsxfun(@times, masked, regsh .* dV),1));
            runs.water.mix.shelf(tt:tend) = full(sum( ...
                bsxfun(@times, maskmix, regsh .* dV),1));

            % statistics of shelf water on the slope
            runs.water.sh.xslope(tt:tend) = full(sum( ...
                bsxfun(@times, masksh, regsl .* xsp .*dV),1)) ./ ...
                (runs.water.sh.slope(tt:tend));
            runs.water.sh.yslope(tt:tend) = full(sum( ...
                bsxfun(@times, masksh, regsl .* ysp .*dV),1)) ./ ...
                (runs.water.sh.slope(tt:tend));
            runs.water.sh.zslope(tt:tend) = full(sum( ...
                bsxfun(@times, masksh, regsl .* zsp .*dV),1)) ./ ...
                (runs.water.sh.slope(tt:tend));

            % statistics of eddy/mix water on shelf.
            % - x location, y-location, z-location
            runs.water.eddmix.xshelf(tt:tend) = full(sum( ...
                bsxfun(@times, (maskmix | masked), regsh .* xsp .*dV),1)) ./ ...
                (runs.water.edd.shelf(tt:tend) + runs.water.mix.shelf(tt:tend));

            runs.water.eddmix.yshelf(tt:tend) = full(sum( ...
                bsxfun(@times, (maskmix | masked), regsh .* ysp .*dV),1)) ./ ...
                (runs.water.edd.shelf(tt:tend) + runs.water.mix.shelf(tt:tend));

            runs.water.eddmix.zshelf(tt:tend) = full(sum( ...
                bsxfun(@times, (maskmix | masked), regsh .* zsp .*dV),1)) ./ ...
                (runs.water.edd.shelf(tt:tend) + runs.water.mix.shelf(tt:tend));

            % how uniform is the "plume" in the vertical - i.e.,
            % baroclinicity - RMS (tracer) / RMS (depth avg tracer)



        end
        toc(ticstart);

        time = runs.time/86400;

        runs.water.comment = [''];

        water = runs.water;


        %%
        % calculate total classified volume
        if debug
            water.totvol = sum(dV(:));
            masses = fieldnames(runs.water);
            classvol = zeros(size(runs.water.mix.deep));
            for ii=1:length(masses)
                if strcmpi(masses{ii},'eddmix'), continue; end
                try
                    regions = fieldnames(runs.water.(masses{ii}));
                    for jj=1:length(regions)
                        if regexp(regions{jj},'^[xyz]'), continue; end
                        classvol = classvol + runs.water.(masses{ii}).(regions{jj});
                    end
                catch ME
                end
            end
            water.classvol = classvol;
            figure;
            plot(water.totvol - water.classvol);
        end
        save([runs.dir '/watermass.mat'], 'water');
    end

    % eddy bulk properties - integrated PV, RV, volume, energy
    function [] = eddy_bulkproperties(runs)
        %%
        slab = 10; % read 10 at a time
        nt = size(runs.zeta, 3);

        sz4dfull = [fliplr(size(runs.rgrid.z_r))-[2 2 0] slab];
        sz4dsp = [prod(sz4dfull(1:3)) slab];
        sz4dspend = [sz4dsp(1) mod(nt,slab)];
        sz3dsp = [sz4dsp(1) 1];

        szpvfull = [fliplr(size(runs.rgrid.z_r))-[2 2 1] slab];
        szpvsp = [prod(szpvfull(1:3)) slab];
        szpvspend = [szpvsp(1) mod(nt,slab)];
        szpv3dsp = [szpvsp(1) 1];

        dVsp = reshape(runs.rgrid.dV(2:end-1,2:end-1,:), sz3dsp);
        dVpvsp = reshape( avg1(runs.rgrid.dV(2:end-1, 2:end-1,:),3), szpv3dsp);

        % store variables to optimize parfor loop
        rgrid = runs.rgrid;
        vormask = runs.eddy.vormask;
        eddname = runs.eddname;
        dirname = runs.dir;
        thresh = runs.eddy_thresh;
        N = runs.rgrid.N;
        zr = permute(runs.rgrid.z_r(:, 2:end-1, 2:end-1), [3 2 1]);

        pvname = [runs.dir '/ocean_vor.nc'];
        if exist(pvname,'file')
            dopv = 1;
        else
            dopv = 0;
        end

        % background density field
        if runs.bathy.axis == 'y'
            tback = permute( dc_roms_read_data(dirname, 'temp', [1 1], ...
                    {'x' 1 1; 'y' 2 sz4dfull(2)+1}, ...
                    [], rgrid), [3 1 2]);
        else
            error('Not implemented for N-S isobaths');
        end

        ticstart = tic;
        for mm=1:ceil(size(runs.zeta,3)/slab)
            tt = (mm-1)*slab + 1;
            disp([' mm= ' num2str(mm) '/' num2str( ...
                    ceil(size(runs.zeta,3)/slab))]);
            tend = tt + slab - 1;
            if tend > nt
                tend = nt;
                sz = sz4dspend;
                szpv = szpvspend;
            else
                sz = sz4dsp;
                szpv = szpvsp;
            end
            eddye = dc_roms_read_data(dirname, eddname, ...
                    [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
                    [],rgrid); %#ok<*PROP>

            masked  = sparse(reshape(eddye > thresh, sz));
            maskvor = sparse(reshape( repmat( ...
                    permute(vormask(:,:,tt:tend), [1 2 4 3]), ...
                    [1 1 N 1]), sz));

            %vol{tt} = runs.domain_integratesp(masked.*maskvor, dVsp);
            % calculate total volume
            volcell{mm} = full(nansum( bsxfun(@times, masked.*maskvor, dVsp)));

            % integrated energies
            % yes, using sz4dfull(1)+1 IS correct. sz4dfull has interior
            % RHO point counts
            u = avg1(dc_roms_read_data(dirname, 'u', ...
                    [tt tend],{'y' 2 sz4dfull(2)+1}, ...
                    [],rgrid),1) - runs.params.bg.ubt; %#ok<*PROP>
            v = avg1(dc_roms_read_data(dirname, 'v', ...
                    [tt tend],{'x' 2 sz4dfull(1)+1}, ...
                    [],rgrid),2) - runs.params.bg.vbt; %#ok<*PROP>

            temp = dc_roms_read_data(dirname, 'temp', ...
                    [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
                    [], rgrid);

            pe = - runs.params.phys.TCOEF* bsxfun(@times, ...
                        bsxfun(@minus, temp, tback), zr)  ...
                    .* runs.params.phys.g;

            intpe{mm} = full(nansum( bsxfun(@times, ...
                        masked.*maskvor.*reshape(pe, sz), dVsp)));

            intke{mm} = full(nansum( bsxfun(@times, ...
                    masked.*maskvor.*reshape(0.5 * (u.^2 + v.^2), sz), dVsp)));

            % integrated PV, RV
            if dopv
                disp('Reading pv, rv');
                pv = ncread(pvname, 'pv',[1 1 1 tt],[Inf Inf Inf tend-tt+1]); %#ok<*PROP>
                rv = avg1(avg1(ncread(pvname, 'rv',[1 1 1 tt], ...
                            [Inf Inf Inf tend-tt+1]),1),2);
                % pv,rv are at N-1 levels in vertical, so we need
                % to calculate masks again
                masked  = sparse(reshape(avg1(eddye,3) > thresh, szpv));
                maskvor = sparse(reshape( repmat( ...
                            permute(vormask(:,:,tt:tend), [1 2 4 3]), ...
                            [1 1 N-1 1]), szpv));

                intpv{mm} = full(nansum( bsxfun(@times, ....
                        masked.*maskvor.*reshape(pv, szpv), dVpvsp))) ...
                        ./ volcell{mm};
                intrv{mm} = full(nansum( bsxfun(@times, ....
                    masked.*maskvor.*reshape(rv, szpv), dVpvsp))) ...
                    ./ volcell{mm};
            end
        end
        toc(ticstart);

        % save data to structure
        eddy = runs.eddy;
        eddy.vol = cell2mat(volcell);
        eddy.PV = cell2mat(intpv);
        eddy.RV = cell2mat(intrv);
        eddy.KE = cell2mat(intke);
        eddy.PE = cell2mat(intpe);

        save([runs.dir '/eddytrack.mat'],'eddy');
    end

    % domain integration for sparse matrix input
    function [out] = domain_integratesp(runs,in, dV)

        if ~exist('dV','var') % not good idea
            dV = reshape(runs.rgrid.dV, runs.streamer.sz3dsp);
        end

        out = full(nansum( bsxfun(@times, in, dV)));
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

    function [] = distrib_csz_old(runs)
        % lets subtract out mean at each z-level to account for near
        % surface and near bottom upwelling.
        % This has to be done after interpolating to constant z-level
        % because you can't take a constant z-level mean otherwise
        yend = find_approx(runs.rgrid.y_rho(:,1),100*1000);
        t0 = runs.eddy.trevind;
        read_start = [1 1 1 t0];
        read_count = [Inf yend Inf 35];

        zdye = ncread(runs.out_file,runs.zdname, ...
                        read_start,read_count);
        csdye = ncread(runs.out_file,runs.csdname, ...
                        read_start, read_count)/1000;
        w = ncread(runs.out_file,'w', ...
                        read_start, read_count);
        %asdye = double(ncread(runs.out_file,runs.asdname, ...
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
            tic;
            runs.update_zeta(hz,ii);
            runs.update_eddy_contour(he,ii);
            runs.update_title(ht,titlestr,ii);
            runs.video_update();
            toc;
            pause(0.03);
        end
        runs.video_write();
    end

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

    function [] = animate_3d(runs)
        stride = [1 1 1 1];

        xrmat = repmat(runs.rgrid.xr(1:stride(1):end,1:stride(2):end)', ...
                        [1 1 runs.rgrid.N]);
        yrmat = repmat(runs.rgrid.yr(1:stride(1):end,1:stride(2):end)', ...
                        [1 1 runs.rgrid.N]);
        zrmat = permute(runs.rgrid.zr(1:stride(1):end,1:stride(2):end,:),[2 1 3]);

        %eddye = roms_read_data(runs.dir,runs.eddname,[1 1 1 1],[Inf Inf Inf Inf], ...
        %            stride);

        tic;csdye = ncread(runs.out_file,runs.csdname);toc;
        tic;eddye = ncread(runs.out_file,runs.eddname);toc;
        csdye = permute(csdye,[2 1 3 4]);
        eddye = permute(eddye,[2 1 3 4]);
        mask = zeros(size(eddye,2),size(eddye,1),size(eddye,4));
        mask(2:end-1,2:end-1,:)=runs.eddy.vormask;
        mask = 1 + zeros(size(mask));

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
        % WRONGGGGGGGGGGGGGGGGG!
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
        if ~exist('tind','var')
            tind = 1;
        end
%             if ~exist([runs.dir '/ocean_vor.nc'],'file')
%                 dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
%             end

        runs.video_init('vor');

        debug = 0;

        %%
        znew = linspace(min(runs.rgrid.z_r(:)), max(runs.rgrid.z_r(:)), 50)';
        % prepare grids for differentiation
        xvor = avg1(avg1(runs.rgrid.xr,1),2);
        yvor = avg1(avg1(runs.rgrid.yr,1),2);
        N = runs.rgrid.N;
        Nnew = length(znew);

        % setup grid
        [sx sy] = size(runs.rgrid.x_rho');
        gridu.xmat = repmat(runs.rgrid.x_u',[1 1 Nnew]);
        gridu.ymat = repmat(runs.rgrid.y_u',[1 1 Nnew]);
        gridu.zmat = permute(runs.rgrid.z_u,[3 2 1]);
        gridu.znew = repmat(permute(znew,[ 3 2 1]),[sx-1 sy 1]);
        %gridu.s = runs.rgrid.s_rho;
        %gridu.zw = runs.rgrid.z_w;
        %gridu.s_w = runs.rgrid.s_w;

        gridv.xmat = repmat(runs.rgrid.x_v',[1 1 Nnew]);
        gridv.ymat = repmat(runs.rgrid.y_v',[1 1 Nnew]);
        gridv.zmat = permute(runs.rgrid.z_v,[3 2 1]);
        gridv.znew = repmat(permute(znew,[ 3 2 1]),[sx sy-1 1]);
        %gridv.s = runs.rgrid.s_rho;
        %gridv.zw = runs.rgrid.z_w;
        %gridv.s_w = runs.rgrid.s_w;

        gridr.xmat = repmat(runs.rgrid.x_rho',[1 1 Nnew]);
        gridr.ymat = repmat(runs.rgrid.y_rho',[1 1 Nnew]);
        gridr.zmat = permute(runs.rgrid.z_r,[3 2 1]);
        gridr.znew = repmat(permute(znew,[3 2 1]),[sx sy 1]);
        %gridr.s = runs.rgrid.s_rho;
        %gridr.zw = runs.rgrid.z_r;
        %gridr.s_w = runs.rgrid.s_w;

        gridw.xmat = repmat(runs.rgrid.x_rho',[1 1 Nnew]);
        gridw.ymat = repmat(runs.rgrid.y_rho',[1 1 Nnew]);
        gridw.zmat = permute(runs.rgrid.z_w,[3 2 1]);
        gridw.znew = repmat(permute(znew,[3 2 1]),[sx sy 1]);
        %gridw.s = runs.rgrid.s_w;
        %gridw.zw = runs.rgrid.z_w;
        %gridw.s_w = runs.rgrid.s_w;

        gridrv.xmat = repmat(xvor,[1 1 Nnew]);
        gridrv.ymat = repmat(yvor,[1 1 Nnew]);
        gridrv.zmat = avg1(avg1(avg1(permute(runs.rgrid.z_r,[3 2 1]),1),2),3);
        gridrv.znew = repmat(permute(znew,[3 2 1]),[sx-1 sy-1 1]);
        gridrv.s = avg1(runs.rgrid.s_rho);
        gridrv.zw = avg1(avg1(avg1(permute(runs.rgrid.z_w,[3 2 1]),1),2),3);
        gridrv.s_w = avg1(runs.rgrid.s_w);

        % for depth integration
        h = runs.bathy.h(2:end-1,2:end-1);
        csr = runs.rgrid.Cs_r(2:end-1);
        csw = runs.rgrid.Cs_w(2:end-1);

        xavg = avg1(avg1(xvor,1),2)/1000; yavg = avg1(avg1(yvor,1),2)/1000;

        depthRange = [100 -max(runs.bathy.h(:))];
        trange = tind:2:size(runs.zeta,3);

        disp(['starting from t instant = ' num2str(trange(1))]);
        runs.vorbudget.hadv = nan(length(trange)-1);
        runs.vorbudget.vadv = runs.vorbudget.hadv;
        runs.vorbudget.tilt = runs.vorbudget.hadv;
        runs.vorbudget.str  = runs.vorbudget.hadv;
        runs.vorbudget.beta = runs.vorbudget.hadv;
        %runs.vorbudget.sol = runs.vorbudget.hadv;
        runs.vorbudget.budget = runs.vorbudget.hadv;
        %runs.vorbudget.conthis = runs.vorbudget.hadv;
        %%
        for kk=1:length(trange)-1
            tt = trange(kk);
            zeta = runs.zeta(2:end-1,2:end-1,tt);

            % read data
            %fname = [runs.dir '/ocean_his.nc.new2'];
            %fname = runs.out_file;
            %u1 = double(ncread(fname,'u',[1 1 1 tt],[Inf Inf Inf 2]));
            %v1 = double(ncread(fname,'v',[1 1 1 tt],[Inf Inf Inf 2]));
            %w  = double(ncread(fname,'w',[1 1 1 tt],[Inf Inf Inf 1]));
            %zeta = double(ncread(fname,'zeta',[1 1 tt],[Inf Inf 1]));

            u = dc_roms_read_data(runs.dir,'u',tt,{},[],runs.rgrid);
            v = dc_roms_read_data(runs.dir,'v',tt,{},[],runs.rgrid);
            w = dc_roms_read_data(runs.dir,'w',tt,{},[],runs.rgrid);
            rho = dc_roms_read_data(runs.dir,'rho',tt,{},[],runs.rgrid);

            % interpolate to znew depths
            disp('interpolating variables');
            u = interpolate(u, gridu.zmat, znew);
            v = interpolate(v, gridv.zmat, znew);
            w = interpolate(w, gridw.zmat, znew);
            rho = interpolate(rho, gridr.zmat, znew);

            ux = diff(u,1,1)./diff(gridu.xmat,1,1);
            uy = diff(u,1,2)./diff(gridu.ymat,1,2);
            uz = diff(u,1,3)./diff(gridu.znew,1,3);

            vx = diff(v,1,1)./diff(gridv.xmat,1,1);
            vy = diff(v,1,2)./diff(gridv.ymat,1,2);
            vz = diff(v,1,3)./diff(gridv.znew,1,3);

            wx = diff(w,1,1)./diff(gridw.xmat,1,1);
            wy = diff(w,1,2)./diff(gridw.ymat,1,2);
            wz = diff(w,1,3)./diff(gridw.znew,1,3);

            rx = diff(rho,1,1)./diff(gridr.xmat,1,1);
            ry = diff(rho,1,2)./diff(gridr.ymat,1,2);
            rz = diff(rho,1,3)./diff(gridr.znew,1,3);

            % tendency term code - not really needed since it is probably a
            % bad estimate when using daily snapshots .
%             if debug
%                 u1 = interpolate(u1, gridu.zmat, znew);
%                 v1 = interpolate(v1, gridv.zmat, znew);
%                 v1x = diff(v1,1,1)./diff(gridv.xmat,1,1);
%                 u1y = diff(u1,1,2)./diff(gridu.ymat,1,2);
%                 rv1 = v1x-u1y;
%             end
            rv = vx-uy;
            rvx = diff(rv,1,1)./diff(gridrv.xmat,1,1);
            rvy = diff(rv,1,2)./diff(gridrv.ymat,1,2);
            rvz = diff(rv,1,3)./diff(gridrv.znew,1,3);

            str = avg1(-1 * avg1(avg1(bsxfun(@plus,rv,avg1(avg1(runs.rgrid.f',1),2)),1),2) ...
                            .* (ux(:,2:end-1,:,:) + vy(2:end-1,:,:)),3);
            tilt = -1 * avg1(avg1( avg1(avg1(wx,2),3) .* avg1(vz,1) + ...
                    avg1(avg1(wy,1),3) .* avg1(uz,2) ,1),2);
            beta = avg1(avg1(runs.params.phys.beta * v(2:end-1,:,:),2),3);
            hadv = avg1( avg1(u(:,2:end-1,:),1) .* avg1(rvx,2) + ...
                    avg1(v(2:end-1,:,:),2) .* avg1(rvy,1),3);
            vadv = avg1(avg1( avg1(avg1(avg1(w,1),2),3) .* rvz ,1),2);

          %  sol = -runs.params.phys.g/runs.params.phys.rho0 .* ...
          %          ( avg1(rx,2) .* avg1(zy,1) - avg1(ry,1) .* avg1(zx,2));

            zint = avg1(znew);
            RV   = trapz(zint, repnan(rv,0), 3);
            STR  = trapz(zint, repnan(str,0), 3);
            TILT = trapz(zint, repnan(tilt,0), 3);
            BETA = trapz(zint, repnan(beta,0), 3);
            HADV = trapz(zint, repnan(hadv,0), 3);
            VADV = trapz(zint, repnan(vadv,0), 3);

            BUD = trapz(zint, repnan( str+tilt - beta - hadv -vadv,0), 3);

            limy = [0 150];
            titlestr = 'Depth avg rvor';
            % plot
            if kk == 1
                figure; maximize();
                ax(1) = subplot(2,4,[1:2]);
                hvor = pcolor(xavg,yavg,RV); hold on; shading flat;
                axis image;
                ht = runs.set_title('Depth avg rvor',tt);
                he(1) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                shading flat
                caxis([-1 1] * max(abs(RV(:)))); colorbar;
                ylim(limy);

                ax(2) = subplot(2,4,3);
                hbet = pcolor(xavg,yavg,-BETA); colorbar; shading flat;
                he(2) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');
                caxis([-1 1] * max(abs(BETA(:))));
                title('- \beta V');

                ax(3) = subplot(2,4,4); cla
                xran = 1:6:size(xavg,1); yran = 1:4:size(yavg,2);
                hquiv = quiver(xavg(xran,yran),yavg(xran,yran), ...
                        ubar(xran,yran), vbar(xran,yran),1.5);
                title('(ubar,vbar)');
                he(3) = runs.plot_eddy_contour('contour',tt);
                hbathy = runs.plot_bathy('contour','k');

%                 ax(4) = subplot(2,4,5);
%                 htend = pcolor(xavg,yavg,TEND); colorbar; shading flat;
%                 he(4) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 caxis([-1 1] * max(abs(TEND(:))));
%                 title('d\xi/dt');

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
                spaceplots(0.06*ones([1 4]),0.05*ones([1 2]))
                linkaxes(ax,'xy');
                runs.video_update();
                pause();
            else
                set(hvor ,'cdata',RV);
                set(hadv ,'cdata',-ADV);
                set(hbet ,'cdata',-BETA);
                set(hstr ,'cdata',STR);
                set(htilt,'cdata',TILT);
                %set(htend,'cdata',TEND);
                try
                    set(hquiv,'udata',ubar(xran,yran),'vdata',vbar(xran,yran));
                catch ME
                end

                runs.update_eddy_contour(he,tt);
                runs.update_title(ht,titlestr,tt);
                runs.video_update();
                pause(0.01);
            end
        end
        runs.video_write();

        runs.vorbudget.time = runs.time(trange(1:end-1));
        plot(runs.vorbudget.time,-runs.vorbudget.hadv,'r'); hold on
        plot(runs.vorbudget.time,-runs.vorbudget.vadv,'g');
        plot(runs.vorbudget.time,runs.vorbudget.tilt,'b');
        plot(runs.vorbudget.time,runs.vorbudget.str,'c');
        plot(runs.vorbudget.time,runs.vorbudget.sol,'m');
        plot(runs.vorbudget.time,-runs.vorbudget.beta,'y');
        plot(runs.vorbudget.time,runs.vorbudget.budget,'k');
        title('signs so that all terms are on RHS and tendency is LHS');
        legend('hadv','vadv','tilt','str','sol','beta','budget');


        vorbudget = runs.vorbudget;
        save([runs.dir '/vorbudget.mat'],'vorbudget');

    end

    function [] = animate_vorbudget_deprecated(runs,tind)
            if ~exist('tind','var')
                tind = 1;
            end
%             if ~exist([runs.dir '/ocean_vor.nc'],'file')
%                 dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
%             end

        runs.video_init('vor');
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


        gridr.xmat = repmat(runs.rgrid.x_rho',[1 1 N]);
        gridr.ymat = repmat(runs.rgrid.y_rho',[1 1 N]);
        gridr.zmat = permute(runs.rgrid.z_r,[3 2 1]);
        gridr.s = runs.rgrid.s_rho;
        gridr.zw = runs.rgrid.z_r;
        gridr.s_w = runs.rgrid.s_w;

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
        trange = tind:2:size(runs.zeta,3);

        disp(['starting from t instant = ' num2str(trange(1))]);
        runs.vorbudget.hadvtot = nan(length(trange)-1);
        runs.vorbudget.vadvtot = runs.vorbudget.hadvtot;
        runs.vorbudget.tilttot = runs.vorbudget.hadvtot;
        runs.vorbudget.strtot  = runs.vorbudget.hadvtot;
        runs.vorbudget.betatot = runs.vorbudget.hadvtot;
        runs.vorbudget.soltot = runs.vorbudget.hadvtot;
        runs.vorbudget.budgettot = runs.vorbudget.hadvtot;
        runs.vorbudget.conthistot = runs.vorbudget.hadvtot;
        for kk=1:length(trange)-1
            tt = trange(kk);
            zeta = runs.zeta(2:end-1,2:end-1,tt);

            % read data
            %fname = [runs.dir '/ocean_his.nc.new2'];
            %fname = runs.out_file;
            %u1 = double(ncread(fname,'u',[1 1 1 tt],[Inf Inf Inf 2]));
            %v1 = double(ncread(fname,'v',[1 1 1 tt],[Inf Inf Inf 2]));
            %w  = double(ncread(fname,'w',[1 1 1 tt],[Inf Inf Inf 1]));
            %zeta = double(ncread(fname,'zeta',[1 1 tt],[Inf Inf 1]));

            u1 = dc_roms_read_data(runs.dir,'u',[tt tt+1],{},[],runs.rgrid);
            v1 = dc_roms_read_data(runs.dir,'v',[tt tt+1],{},[],runs.rgrid);
            w =  dc_roms_read_data(runs.dir,'w',tt,{},[],runs.rgrid);
            rho = dc_roms_read_data(runs.dir,'rho',tt,{},[],runs.rgrid);

            u = u1(:,:,:,1); v = v1(:,:,:,1);
            u1(:,:,:,1) = []; v1(:,:,:,1) = [];

            % get Hz
            Hz  = diff(set_depth(runs.rgrid.Vtransform, runs.rgrid.Vstretching, ...
                    runs.rgrid.theta_s, runs.rgrid.theta_b, runs.rgrid.hc, ...
                    runs.rgrid.N, 5, runs.rgrid.h', runs.zeta(:,:,tt),0),1,3);

            try
                % ROMS outputs omega as Hz.*Ds/Dt. rescale to get Ds/Dt
                omega = avg1(dc_roms_read_data(runs.dir,'omega',tt,{},[],runs.rgrid),3) ...
                            ./ Hz;
                %omega = double(ncread(fname,'omega',[1 1 1 tt],[Inf Inf Inf 1]));
            catch ME
                udzdx = avg1(u,1) .* diff(gridu.zmat,1,1)./diff(gridu.xmat,1,1);
                vdzdy = avg1(v,2) .* diff(gridv.zmat,1,2)./diff(gridv.ymat,1,2);
                % this is a good estimate - problem areas are in the sponge
                % FACTOR OF HZ?
                omega = avg1(w(2:end-1,2:end-1,:),3)  ...
                        - udzdx(:,2:end-1,:) - vdzdy(2:end-1,:,:);
                % this is a good estimate - problem areas are in the sponge
                w2 = udzdx(:,2:end-1,:) + vdzdy(2:end-1,:,:) + ...
                     omega;
            end

            %rvor1 = double(ncread(runs.out_file,'rvorticity',[1 1 1 tt], ...
            %    [Inf Inf Inf 2]));
            %rvor = rvor1(:,:,:,1);

            %% z co-ordinate
            % differentiate


        gridrv.xmat = repmat(xvor,[1 1 N-1]);
        gridrv.ymat = repmat(yvor,[1 1 N-1]);
        gridrv.zmat = avg1(avg1(avg1(permute(runs.rgrid.z_r,[3 2 1]),1),2),3);
        gridrv.s = avg1(runs.rgrid.s_rho);
        gridrv.zw = avg1(avg1(avg1(permute(runs.rgrid.z_w,[3 2 1]),1),2),3);
        gridrv.s_w = avg1(runs.rgrid.s_w);

%
%             ux = diff_cgrid(gridu,u,1); uy = diff_cgrid(gridu,u,2);
%                 uz = diff_cgrid(gridu,u,3);
%             vx = diff_cgrid(gridv,v,1); vy = diff_cgrid(gridv,v,2);
%                 vz = diff_cgrid(gridv,v,3);
%             wx = avg1(diff_cgrid(gridw,w,1),3); wy = avg1(diff_cgrid(gridw,w,2),3);
%                 wz = avg1(diff_cgrid(gridw,w,3),3);
%
%             % calculate relative vorticity
%             v1x = diff_cgrid(gridv,v1,1); u1y = diff_cgrid(gridu,u1,2);
%             rvor = vx-uy; rv1 = v1x-u1y;
%             rvx = diff_cgrid(gridrv,rvor,1); rvy = diff_cgrid(gridrv,rvor,2);
%                 rvz = diff_cgrid(gridrv,rvor,3);
%
%             % check continuity
%             % THIS DOESN't WORK with HISTORY FILES but does average files
%             % better than CONTHIS
%             cont = ux(:,2:end-1,:) + vy(2:end-1,:,:) + wz(2:end-1,2:end-1,:);
%             CONT = sum(cont .* avg1(runs.rgrid.dV(2:end-1,2:end-1,:),3),3) ./ ...
%                      sum(avg1(runs.rgrid.dV(2:end-1,2:end-1,:),3),3);
%
%             % form terms - avg to interior RHO points
%             % in z co-ordinates
%             adv = avg1( avg1(avg1(u(:,:,2:end-1),1),2) .* rvx,2) + ...
%                     avg1( avg1(avg1(v(:,:,2:end-1),1),2) .* rvy,1) + ...
%                         avg1(avg1( avg1(avg1(avg1(w(:,:,2:end-1),1),2),3) ...
%                                .* rvz    ,1),2);
%             str = avg1(avg1( ...
%                     avg1(   bsxfun(@plus,rvor,avg1(avg1(runs.rgrid.f',1),2)) ...
%                             .* avg1(avg1(wz,1),2)   ,1) ...
%                                 ,2),3);
%
%             bet = avg1(beta * v(2:end-1,:,2:end-1),2);
%
%             tilt = avg1( avg1(avg1(avg1(wy,1).*avg1(uz,2) - ...
%                         avg1(wx,2).*avg1(vz,1),1),2),3);
%
%             tend = avg1(avg1( ...
%                     avg1(rv1-rvor,3)./diff(runs.time(1:2)) ,1),2);
%
%
%             % depth integrate
%             [ubar,vbar] = uv_barotropic(u,v,Hz);
%             [rint,ravg] = roms_depthIntegrate(avg1(avg1(rvor,1),2), ...
%                             csr,csw,h,zeta,depthRange);
%             [~,ADV] = roms_depthIntegrate(adv ,csr,csw,h,zeta,depthRange);
%             [~,STR] = roms_depthIntegrate(str ,csr,csw,h,zeta,depthRange);
%             [~,BET] = roms_depthIntegrate(bet ,csr,csw,h,zeta,depthRange);
%             [~,TILT] = roms_depthIntegrate(tilt,csr,csw,h,zeta,depthRange);
%             [~,TEND] = roms_depthIntegrate(tend,csr,csw,h,zeta,depthRange);
%             rplot = ravg;
%
%             budget = tend+adv+bet-tilt-str;
%             BUD = TEND+ADV+BET-TILT-STR;

           %% in s  co-ordinates
            % this works on history file OUTSIDE THE SPONGE
            %- not so well when I estimate omega from w
            duHzdx = diff(u .* avg1(Hz,1),1,1)./diff(gridu.xmat,1,1);
            dvHzdy = diff(v .* avg1(Hz,2),1,2)./diff(gridv.ymat,1,2);
            doHzds = diff(omega .* Hz,1,3);
            try
                conthis = avg1(duHzdx(:,2:end-1,:)  ...
                        + dvHzdy(2:end-1,:,:),3) + doHzds(2:end-1,2:end-1,:);
                CONTHIS =  sum((conthis .* avg1(runs.rgrid.dV(2:end-1,2:end-1,:),3)),3) ./ ...
                        sum(runs.rgrid.dV(2:end-1,2:end-1,:),3);
            catch ME
                conthis = avg1(duHzdx(:,2:end-1,:)  ...
                        + dvHzdy(2:end-1,:,:),3) + doHzds;
                CONTHIS =  sum((conthis .* avg1(runs.rgrid.dV(2:end-1,2:end-1,:),3)),3) ./ ...
                    avg1(sum(runs.rgrid.dV(2:end-1,2:end-1,:),3),3);
            end

            ux = diff(u,1,1)./diff(gridu.xmat,1,1);
            uy = diff(u,1,2)./diff(gridu.ymat,1,2);
            us = bsxfun(@rdivide,diff(u,1,3), ...
                    diff(permute(gridu.s',[3 2 1]),1,3));

            u1y = diff(u1,1,2)./diff(gridu.ymat,1,2);
            v1x = diff(v1,1,1)./diff(gridv.xmat,1,1);

            vx = diff(v,1,1)./diff(gridv.xmat,1,1);
            vy = diff(v,1,2)./diff(gridv.ymat,1,2);
            vs = bsxfun(@rdivide,diff(v,1,3), ...
                diff(permute(gridv.s',[3 2 1]),1,3));

            ox = diff(omega,1,1)./diff(gridr.xmat,1,1);
            oy = diff(omega,1,2)./diff(gridr.ymat,1,2);
           % os = bsxfun(@rdivide,diff(omega,1,3), ...
           %         diff(permute(gridr.s',[3 2 1]),1,3));

            rhox = diff(rho,1,1)./diff(gridr.xmat,1,1);
            rhoy = diff(rho,1,2)./diff(gridr.ymat,1,2);

            zx = diff(gridr.zmat,1,1)./diff(gridr.xmat,1,1);
            zy = diff(gridr.zmat,1,2)./diff(gridr.ymat,1,2);

            rv = vx-uy; rv1 = v1x-u1y;

            %if kk == 1
                gridrv.xmat(:,:,end+1) = gridrv.xmat(:,:,1);
                gridrv.ymat(:,:,end+1) = gridrv.ymat(:,:,1);
                gridrv.s = runs.rgrid.s_rho;
            %end
            rvx = diff(rv,1,1)./diff(gridrv.xmat,1,1);
            rvy = diff(rv,1,2)./diff(gridrv.ymat,1,2);
            rvs = bsxfun(@rdivide,diff(rv,1,3), ...
                    diff(permute(gridrv.s',[3 2 1]),1,3));

            tend = (rv1-rv)./diff(runs.time(1:2));

            hadv = avg1(avg1(u(:,2:end-1,:),1) .* avg1(rvx,2) + ...
                    avg1(v(2:end-1,:,:),2) .* avg1(rvy,1),3);

            vadv = avg1(avg1(avg1(omega,1),2),3) .* rvs;

            adv = hadv+avg1(avg1(vadv,1),2);

            str = -1 .* avg1(avg1(bsxfun(@plus,rv,avg1(avg1(runs.rgrid.f',1),2)),1),2) ...
                        .* (ux(:,2:end-1,:) + vy(2:end-1,:,:));

            tilt = avg1(avg1(oy,1),3) .* avg1(us,2) -  ...
                    avg1(avg1(ox,2),3) .* avg1(vs,1);

            bet = beta * avg1(v(2:end-1,:,:),2);

            sol = -runs.params.phys.g/runs.params.phys.rho0 .* ...
                    ( avg1(rhox,2) .* avg1(zy,1) - avg1(rhoy,1) .* avg1(zx,2));

            budget = avg1(avg1(avg1(tend - sol,1),2)+bet-str,3) + ...
                        hadv + avg1(avg1(vadv - tilt,1),2);


            runs.vorbudget.hadvtot(:,kk) = sum(hadv(:));
            runs.vorbudget.vadvtot(:,kk) = sum(vadv(:));
            runs.vorbudget.betatot(:,kk) = sum(bet(:));
            runs.vorbudget.soltot(:,kk) = sum(sol(:));
            runs.vorbudget.strtot(:,kk) = sum(str(:));
            runs.vorbudget.tilttot(:,kk) = sum(tilt(:));
            runs.vorbudget.budgettot(:,kk) = sum(budget(:));
            runs.vorbudget.conthistot(:,kk) = sum(conthis(:));

%            ubar = avg1(ubar(:,2:end-1),1);
%            vbar = avg1(vbar(2:end-1,:),2);

            limy = [0 150];
            titlestr = 'Depth avg rvor';
            % plot
%             if kk == 1
%                 figure; maximize();
%                 ax(1) = subplot(2,4,[1:2]);
%                 hvor = pcolor(xavg,yavg,rplot); hold on; shading flat;
%                 axis image;
%                 ht = runs.set_title('Depth avg rvor',tt);
%                 he(1) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 shading flat
%                 caxis([-1 1] * max(abs(rplot(:)))); colorbar;
%                 ylim(limy);
%
%                 ax(2) = subplot(2,4,3);
%                 hbet = pcolor(xavg,yavg,-BET); colorbar; shading flat;
%                 he(2) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 caxis([-1 1] * max(abs(BET(:))));
%                 title('- \beta V');
%
%                 ax(3) = subplot(2,4,4); cla
%                 xran = 1:6:size(xavg,1); yran = 1:4:size(yavg,2);
%                 hquiv = quiver(xavg(xran,yran),yavg(xran,yran), ...
%                         ubar(xran,yran), vbar(xran,yran),1.5);
%                 title('(ubar,vbar)');
%                 he(3) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%
%                 ax(4) = subplot(2,4,5);
%                 htend = pcolor(xavg,yavg,TEND); colorbar; shading flat;
%                 he(4) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 caxis([-1 1] * max(abs(TEND(:))));
%                 title('d\xi/dt');
%
%                 ax(5) = subplot(2,4,6);
%                 hadv = pcolor(xavg,yavg,-ADV); colorbar; shading flat;
%                 he(5) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 caxis([-1 1] * max(abs(ADV(:))));
%                 title('-Advection');
%
%                 ax(6) = subplot(2,4,7);
%                 htilt = pcolor(xavg,yavg,TILT); colorbar; shading flat;
%                 he(6) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 caxis([-1 1] * max(abs(TILT(:))));
%                 title('Tilting');
%
%                 ax(7) = subplot(2,4,8);
%                 hstr = pcolor(xavg,yavg,STR); colorbar; hold on; shading flat;
%                 %hquiv = quiverclr(xavg(xran,yran),yavg(xran,yran), ...
%                 %    ubar(xran,yran),vbar(xran,yran),0.3,STR(xran,yran), ...
%                 %    [-1 1]*1e-11);
%                 %set(gca,'color',[0 0 0]);
%                 he(7) = runs.plot_eddy_contour('contour',tt);
%                 hbathy = runs.plot_bathy('contour','k');
%                 caxis([-1 1] * max(abs(STR(:))));
%                 title('Stretching = (f+\xi)w_z')
%                 spaceplots(0.06*ones([1 4]),0.05*ones([1 2]))
%                 linkaxes(ax,'xy');
%                 runs.video_update();
%                 pause();
%             else
%                 set(hvor ,'cdata',rplot);
%                 set(hadv ,'cdata',-ADV);
%                 set(hbet ,'cdata',-BET);
%                 set(hstr ,'cdata',STR);
%                 set(htilt,'cdata',TILT);
%                 set(htend,'cdata',TEND);
%                 try
%                     set(hquiv,'udata',ubar(xran,yran),'vdata',vbar(xran,yran));
%                 catch ME
%                 end
%
%                 runs.update_eddy_contour(he,tt);
%                 runs.update_title(ht,titlestr,tt);
%                 runs.video_update();
%                 pause(0.01);
%             end
        end
%        runs.video_write();

        runs.vorbudget.time = runs.time(trange(1:end-1));
        plot(runs.vorbudget.time,-runs.vorbudget.hadvtot,'r'); hold on
        plot(runs.vorbudget.time,-runs.vorbudget.vadvtot,'g');
        plot(runs.vorbudget.time,runs.vorbudget.tilttot,'b');
        plot(runs.vorbudget.time,runs.vorbudget.strtot,'c');
        plot(runs.vorbudget.time,runs.vorbudget.soltot,'m');
        plot(runs.vorbudget.time,-runs.vorbudget.betatot,'y');
        plot(runs.vorbudget.time,runs.vorbudget.budgettot,'k');
        title('signs so that all terms are on RHS and tendency is LHS');
        legend('hadv','vadv','tilt','str','sol','beta','budget');


        vorbudget = runs.vorbudget;
        save([runs.dir '/vorbudget.mat'],'vorbudget');
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
            if isempty(runs.eddye)
                runs.eddye = dc_roms_read_data(runs.dir,runs.eddname, ...
                    [],{'z' runs.rgrid.N runs.rgrid.N},[],runs.rgrid);
            end
            if isempty(runs.csdye)
                runs.csdye = dc_roms_read_data(runs.dir,runs.csdname, ...
                    [],{'z' runs.rgrid.N runs.rgrid.N},[],runs.rgrid);
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
                dc_roms_read_data(runs.dir,runs.eddname,i,{},[],runs.rgrid), ...
                depth,grdr);
            csdye = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,runs.csdname,i,{},[],runs.rgrid), ...
                depth,grdr);
            u = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,'u',i,{},[],runs.rgrid), depth,grdu);
            v = dc_roms_zslice_var( ...
                dc_roms_read_data(runs.dir,'v',i,{},[],runs.rgrid), depth,grdv);
            % get on interior RHO points
            u = avg1(u(:,2:end-1),1);
            v = avg1(v(2:end-1,:),2);
            % decimate for quiver
            u = u(1:dxi:end,1:dyi:end);
            v = v(1:dxi:end,1:dyi:end);
        end

        % get scale for u,v
        if ~isempty(runs.usurf)
            uref = max(max(abs(runs.usurf(:,:,1))));
        else
            uref = ncread(runs.out_file,'u',[1 1 40 1],[Inf Inf 1 1]);
            uref = max(abs(uref(:)));
        end
        if ~isempty(runs.vsurf)
            vref = max(max(abs(runs.vsurf(:,:,1))));
        else
            vref = ncread(runs.out_file,'v',[1 1 40 1],[Inf Inf 1 1]);
            vref = max(abs(vref(:)));
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
        set(he,'LineColor',[1 0 0],'LineWidth',1)

        hcsdye = pcolor(runs.rgrid.xr/1000,runs.rgrid.yr/1000, ...
                    fillnan((csdye/1000 < clim_csd(2)) .* csdye/1000,0));
        shading flat;
        caxis(clim_csd);colormap(cmcsd); freezeColors;
        hcb2 = colorbar; cbunits(hcb2,'Cross-shore dye');cbfreeze(hcb2);

        hq = quiver(runs.eddy.xr(1:dxi:end,1:dyi:end)/1000,runs.eddy.yr(1:dxi:end,1:dyi:end)/1000, ...
                    u./uref, v./vref);
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
                    dc_roms_read_data(runs.dir,runs.eddname,i), depth,grdr);
                csdye = dc_roms_zslice_var( ...
                    dc_roms_read_data(runs.dir,runs.csdname,i), depth,grdr);
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
            %handle.CData = runs.zeta(:,:,tt);
            set(handle,'CData',runs.zeta(:,:,tt));
        catch ME
            %handle.ZData = runs.zeta(:,:,tt);
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