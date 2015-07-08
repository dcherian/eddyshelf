classdef runs < handle
properties
    % dir & file names
    name; dir; out_file; ltrans_file; flt_file; givenFile; tracpy_file; ...
        fpos_file;
    % data
    zeta; temp; usurf; vsurf; vorsurf; csdsurf; ubot; vbot; eddsurf; ...
        rhosurf; edcsdyesurf;
    rbacksurf; % background density at surface
    % dimensional and non-dimensional time
    time; ndtime; tscale; tscaleind;
    % barotropic vel (geostrophic)
    ubar; vbar;
    % dyes
    csdye; asdye; zdye; eddye; % cross-shore, along-shore, z dyes, eddy dye
    % dye names
    csdname; asdname; zdname; eddname;
    % velocity names
    csvelname; asvelname;
    % vorticity budget
    vorbudget;
    % grid & bathymetry
    rgrid; bathy
    % float data
    roms; ltrans; tracpy;
    % sponge details
    sponge; % actual sponge mask
    spng; % contains sponge indices
    % eddy track data
    eddy; noeddy;
    % bottom torque calculations
    bottom; angmom; pbot;
    % tanh trajectory fit
    traj; res;
    % wnoise metric
    wmetric;
    % along-shore jet properties
    jet;
    % threshold values
    eddy_thresh = 0.7;
    % initial params
    params
    % fluxes - cross-shore & along-shore; - energy fluxes
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
methods(Static)
    [y0,X,y1,x0,yc] = tanh_fit(x, y, plot_flag, test);
end
methods
    % constructor
    function [runs] = runs(dir, reset,  do_all)

        read_zeta = 0;
        if ~exist('reset','var')
            reset = 0;
        end

        if ~exist('do_all', 'var')
            do_all = 0;
        end

        if isdir(dir)
            runs.dir = dir;
            files = roms_find_file(dir,'his');
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
        runs.flt_file = [runs.dir '/ocean_flt.nc'];
        runs.ltrans_file = [runs.dir '/ltrans.nc'];
        runs.tracpy_file = [runs.dir '/tracks/tracpy.nc'];
        runs.fpos_file = [runs.dir '/' roms_find_file(dir, 'fpos')];

        % get grid
        zeta0 = double(ncread(runs.out_file,'zeta',[1 1 1],[Inf Inf 1]));
        runs.rgrid = roms_get_grid(runs.out_file,runs.out_file, ...
                        zeta0',1);
        runs.rgrid.xr = runs.rgrid.x_rho';
        runs.rgrid.yr = runs.rgrid.y_rho';
        runs.rgrid.z_r = single(runs.rgrid.z_r);
        runs.rgrid.z_u = single(runs.rgrid.z_u);
        runs.rgrid.z_v = single(runs.rgrid.z_v);
        runs.rgrid.z_w = single(runs.rgrid.z_w);
        %runs.rgrid.zr = permute(runs.rgrid.z_r,[3 2 1]);
        runs.rgrid.z_uw = [];
        runs.rgrid.z_vw = [];
        runs.rgrid.zeta = [];
        runs.rgrid.dx = mean(1./runs.rgrid.pm(:));
        runs.rgrid.dy = mean(1./runs.rgrid.pn(:));

        % remove needless rgrid matrices
        runs.rgrid.angle = [];
        runs.rgrid.x_psi = [];
        runs.rgrid.y_psi = [];
        if runs.rgrid.nomask
            runs.rgrid.mask_rho = [];
                runs.rgrid.mask_rho_nan = [];
            runs.rgrid.mask_psi = [];
            runs.rgrid.mask_u = [];
            runs.rgrid.mask_v = [];
        end
        if runs.rgrid.nolatlon
            runs.rgrid.lon_rho = [];
            runs.rgrid.lon_psi = [];
            runs.rgrid.lon_u = [];
            runs.rgrid.lon_v = [];

            runs.rgrid.lat_rho = [];
            runs.rgrid.lat_psi = [];
            runs.rgrid.lat_u = [];
            runs.rgrid.lat_v = [];
        end

        runs.rbacksurf = ncread(runs.out_file, 'rho', [1 1 runs.rgrid.N ...
                            1], [1 1 1 1]);

        % read zeta
        if ~runs.givenFile
            if read_zeta
                runs.zeta = dc_roms_read_data(dir,'zeta',[],{},[],runs.rgrid, ...
                                          'his', 'single');
            end
            runs.time = dc_roms_read_data(dir,'ocean_time',[],{}, ...
                                          [],runs.rgrid, 'his', 'single');
            %try
            %    runs.csdye  = roms_read_data(dir,runs.csdname, ...
            %        [1 1 runs.rgrid.N 1],[Inf Inf 1 Inf]);
            %catch ME
            %end
        else
            if read_zeta
                runs.zeta = (ncread(runs.out_file,'zeta'));
            end
            runs.time = (ncread(runs.out_file,'ocean_time'));
        end

        runs.rgrid.ocean_time = runs.time;

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
        runs.params.misc = roms_load_misc(runs.out_file);

        if isnan(runs.params.bg.ubt)
            runs.params.bg.ubt = 0;
        end
        if isnan(runs.params.bg.vbt)
            runs.params.bg.vbt = 0;
        end

        % sometimes I forget to change T0 in *.in file
        % but they aren't stored in the history file output, only
        % R0 is
        % runs.params.phys.T0 = ncread(runs.out_file, 'T0');
        % runs.params.phys.R0 = ncread(runs.out_file, 'R0');
        % runs.params.phys.S0 = ncread(runs.out_file, 'S0');

        % fill bathy
        [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = ...
                        find_shelfbreak(runs.out_file);
        [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = ...
                        find_shelfbreak(runs.out_file,'slope');
        runs.bathy.h = runs.rgrid.h';

        % remove background zeta
        if read_zeta
            if runs.bathy.axis == 'x'
                runs.zeta = bsxfun(@minus, runs.zeta, runs.zeta(:,1, ...
                                                                1));
            else
                runs.zeta = bsxfun(@minus, runs.zeta, runs.zeta(1,:, ...
                                                                1));
            end
        end

        % read in sponge
        runs.sponge = ncread(runs.out_file, 'visc2_r') > 0;

        % find sponge indices
        sz = ceil(size(runs.sponge)/2); % mid-points
        runs.spng.sx1 = find(runs.sponge(1:sz(1),sz(2)) == 0, 1, 'first');
        runs.spng.sx2 = sz(1) + find(runs.sponge(sz(1):end,sz(2)) == 1, 1, ...
                             'first') - 2;
        runs.spng.sy1 = find(runs.sponge(sz(1),1:sz(2)) == 0, 1, 'first');
        runs.spng.sy2 = sz(2) + find(runs.sponge(sz(1), sz(2):end) == 1, 1, ...
                                     'first') - 2;

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

        if runs.bathy.axis == 'y'
            runs.asvelname = 'u';
            runs.csvelname = 'v';
        else
            runs.asvelname = 'v';
            runs.csvelname = 'u';
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

        % extra processing of eddy track
        if ~runs.noeddy
            % rerun track_eddy if not new enough
            if ~isfield(runs.eddy,'vor')
                runs.eddy = track_eddy(dir);
            end

            % if gaussian profile then track_eddy fits Lz2. copy to
            % Lgauss for backwards compatibility
            if runs.params.flags.vprof_gaussian && ~runs.noeddy
                if ~isnan(runs.eddy.Lz2)
                    runs.eddy.Lgauss = runs.eddy.Lz2;
                    runs.eddy.Lz2 = nan(size(runs.eddy.Lz2));
                end
            end

            % cvx, cvy for vorticity contour
            runs.eddy.cvx = [0 diff(runs.eddy.vor.cx)./ ...
                diff(runs.eddy.t*86400)];
            runs.eddy.cvy = [0 diff(runs.eddy.vor.cy)./ ...
                diff(runs.eddy.t*86400)];

            runs.params.nondim.eddy.Bu = (runs.params.phys.f0 * ...
                                          runs.params.eddy.dia/2 / ...
                                          runs.params.eddy.depth).^2 / ...
                runs.params.phys.N2;

            % scale time by eddy translation
            if runs.bathy.axis == 'y'
                runs.eddy.tscaleind = find_approx(runs.eddy.my, ...
                                                  runs.bathy.xsl, 1);
            else
                runs.eddy.tscaleind = find_approx(runs.eddy.mx, ...
                                                  runs.bathy.xsl, 1);
            end
            runs.eddy.tscale = runs.eddy.t(runs.eddy.tscaleind) .* ...
                86400;

            % find time when southern edge crosses slopebreak
            runs.eddy.edgtscaleind = find_approx(runs.eddy.vor.se, ...
                                                 runs.bathy.xsl, 1);
            runs.eddy.edgtscale = runs.eddy.t(runs.eddy.edgtscaleind) ...
                * 86400;

            if ~isfield(runs.eddy, 'tend')
                runs.eddy.tend = length(runs.eddy.dia);
            else
                % truncate time vector just in case
                if runs.eddy.tend == length(runs.eddy.dia)
                    runs.eddy.t = runs.eddy.t(1:runs.eddy.tend);
                end
            end

            runs.eddy.Bu = runs.params.phys.N2 .* runs.eddy.Lgauss.^2 ...
                ./ runs.params.phys.f0^2 ./ (runs.eddy.vor.dia/2).^2;

            % eddy turnover time scale
            runs.eddy.turnover = (runs.eddy.vor.dia(1)/2)./runs.eddy.V(1);

            % non-dimensionalized time
            runs.ndtime = runs.time ./ runs.eddy.turnover;

            % save memory by converting masks to logical
            runs.eddy.mask = logical(repnan(runs.eddy.mask, 0));
            runs.eddy.vormask = logical(repnan(runs.eddy.vormask,0));

            % drhothresh based on ssh mask if it doesn't exist
            if ~isfield(runs.eddy, 'drhothreshssh')
                rs = ncread(runs.out_file, 'rho', [1 1 runs.rgrid.N ...
                                    1], [Inf Inf 1 1]);
                rs = rs(2:end-1,2:end-1) - rs(1,1);
                runs.eddy.drhothreshssh = squeeze(nanmax(nanmax(rs .* ...
                                         fillnan(runs.eddy.mask(:,:,1),0), [], 1), [], 2));
            end
            % drhothresh based on ssh mask if it doesn't exist
            if ~isfield(runs.eddy, 'drhothresh')
                rs = ncread(runs.out_file, 'rho', [1 1 runs.rgrid.N ...
                                    1], [Inf Inf 1 1]);
                rs = rs(2:end-1,2:end-1) - rs(1,1);
                runs.eddy.drhothresh = squeeze(nanmax(nanmax(rs .* ...
                                         fillnan(runs.eddy.vormask(:,:,1),0), [], 1), [], 2));
            end
            if isfield(runs.eddy,'cvx')
                if runs.eddy.cvx(1) == 0 || runs.eddy.cvy(1) == 0
                    runs.eddy.cvx(1) = NaN;
                    runs.eddy.cvy(1) = NaN;
                end
            end

            % proximity to shelfbreak
            if runs.bathy.axis == 'y'
                edge = runs.eddy.vor.se;
            else
                edge = runs.eddy.vor.we;
            end
            runs.eddy.prox = (edge-runs.bathy.xsb);

            % time of reversal
            try
                runs.eddy.trevind = find(runs.eddy.cvx < 0,1,'first');
                runs.eddy.trev = runs.time(runs.eddy.trevind);
            catch ME
                disp('Eddy did not reverse direction');
                runs.eddy.trev = nan;
            end
            if isempty(runs.eddy.trev), runs.eddy.trev = NaN; end

            % Early et al (2011) estimates for zonal, meridional
            % velocities
            A = runs.eddy.amp(1);
            Vr = runs.params.phys.beta * runs.rrdeep^2;
            if ~isfield(runs.params.eddy, 'Ldef')
                runs.params.eddy.Ldef = runs.rrdeep;
            end
            Nqg = runs.params.phys.f0 * runs.params.eddy.Ldef / ...
                  runs.params.phys.g * Vr;
            runs.eddy.Vest_zonal = Vr * (Nqg/A - 1);
            runs.eddy.Vest_mer = Vr * Nqg/A;

            % estimate southward vel.
            % (beta * Lr^2)^2 *1/2 * 1/amp * NH/g
            %runs.eddy.Vy = -(runs.params.phys.beta*(runs.params.eddy.dia(1)/2)^2)^2/2 ...
            %                        /runs.eddy.amp(1) * ...
            %                sqrt(runs.params.phys.N2)/runs.params.phys.g*max(runs.bathy.h(:));

%                            -(run3.params.phys.beta*(run3.params.eddy.dia(1)/2)^2)^2 ...
%                                        *1/2 * 1/run3.eddy.amp(1) * ...
%                                sqrt(run3.params.phys.N2)/run3.params.phys.g*max(run3.bathy.h(:))

            % water depth at eddy center
            h = runs.bathy.h(2:end-1,2:end-1);
            ix = vecfind(runs.eddy.xr(:,1), runs.eddy.mx);
            iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.my);
            runs.eddy.hcen = h(sub2ind(size(runs.eddy.xr),ix,iy));

            % f at eddy center
            f = runs.rgrid.f(2:end-1,2:end-1)';
            runs.eddy.fcen = f(sub2ind(size(runs.eddy.xr),ix,iy))';

            if runs.bathy.axis == 'y'
                iy = vecfind(runs.eddy.yr(1,:)',runs.eddy.vor.se);
                runs.eddy.hedge = h(1,iy);
            end

            % remove needless h-matrix
            runs.eddy.h = [];

            if isfield(runs.eddy.vor, 'Ro')
                runs.eddy.Ro = runs.eddy.vor.Ro;
            end

            % calculate βL/f
            dA = 1./runs.rgrid.pm' .* 1./runs.rgrid.pn';
            betal = fillnan(runs.eddy.vormask, 0) .* ...
                    abs(bsxfun(@minus, runs.rgrid.f(2:end-1,2:end-1)', ...
                       permute(runs.eddy.fcen, [3 2 1])))/2;
            betaldA = bsxfun(@times, betal, dA(2:end-1,2:end-1));

            runs.eddy.betahat = squeeze(nansum(nansum(betaldA, 1), 2))' ./ ...
                runs.eddy.vor.area ./ runs.eddy.fcen';
            runs.eddy.Rh = runs.eddy.Ro ./ runs.eddy.betahat;

            if isfield(runs.eddy, 'KE')
                if size(runs.eddy.KE,1) == 1
                    runs.eddy.KE = runs.eddy.KE';
                    runs.eddy.PE = runs.eddy.PE';
                    runs.eddy.vol = runs.eddy.vol';
                end
            end
        end

        if do_all == 1
            runs.asfluxes;
            runs.csfluxes;
            %runs.water_census;
            %runs.jetdetect;
            %runs.eddy_bulkproperties;
        end

        % Process floats
        try
            runs.roms = floats('roms', runs.dir, runs.rgrid, ...
                               runs.bathy.xsb, runs.fpos_file);
            runs.filter_cross_sb('roms');
        catch ME
            disp('Reading ROMS floats failed.');
            disp(ME);
        end
        try
            runs.tracpy = floats('tracpy',runs.tracpy_file,runs.rgrid, ...
                                 runs.bathy.xsb);
            runs.filter_cross_sb('tracpy');
        catch ME
            disp('Reading tracpy floats failed.');
            disp(ME);
        end
        if exist(runs.ltrans_file,'file')
            try
                runs.ltrans = floats('ltrans',runs.ltrans_file, ...
                                     runs.rgrid, runs.bathy.xsb);
            catch ME
                warning('LTRANS data not read in');
            end
            runs.filter_cross_sb('ltrans');
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

        % load vorticity budget data
          % load water mass data
        if exist([dir '/vorbudget.mat'],'file') && reset ~= 1
            disp('Loading vorticity budget');
            vorbudget = load([dir '/vorbudget.mat'], 'vorbudget');
            runs.vorbudget = vorbudget.vorbudget;
            clear vorbudget
        end

        % load fluxes if the file exists
        if exist([dir '/fluxes.mat'],'file') && reset ~= 1
            disp('Loading fluxes');
            data = load([dir '/fluxes.mat']);
            if isfield(data, 'csflux'), runs.csflux = data.csflux; end
            if isfield(data, 'asflux'), runs.asflux = data.asflux; end
            clear data

            % find time when flux is increasing up
            try
                trans  = runs.csflux.west.itrans.shelf(:,1);
                mtrans = max(abs(trans));
                runs.csflux.tscaleind = find_approx(trans, 0.1*mtrans, 1);
                runs.csflux.tscale = runs.csflux.time(runs.csflux.tscaleind);
            catch ME
                warning(['Couldn''t calculate flux based ' ...
                         'timescale']);
            end

            % if asflux.iy hasn't been saved default to between
            % shelfbreak and edge of northern sponge.
            try
                runs.asflux.iy(1);
            catch ME
                runs.asflux.iy = [runs.bathy.isb runs.spng.sy2];
            end

            % create hmat for plotting depth averaged fluxes
            runs.asflux.hmat = repmat(...
                runs.bathy.h(1,runs.asflux.iy(1):runs.asflux.iy(2))', ...
                [1 length(runs.eddy.t)]);
        end

        % load jet diagnostics if the file exists
        if exist([dir '/jet.mat'],'file') && reset ~= 1
            disp('Loading jet diagnostics');
            data = load([dir '/jet.mat']);
            runs.jet = data.jet;
            clear data
        end

        % load bottom torque diagnostics if the file exists
        if exist([dir '/bottom.mat'],'file') && reset ~= 1
            disp('Loading bottom torque diagnostics');
            data = load([dir '/bottom.mat']);
            runs.bottom = data.bottom;
            clear data
        end

        % load angular momentum diagnostics if the file exists
        if exist([dir '/angmom.mat'],'file') && reset ~= 1
            disp('Loading angular momentum diagnostics');
            data = load([dir '/angmom.mat']);
            runs.angmom = data.angmom;
            clear data
        end

        % load trajectory fit
        if exist([dir '/traj.mat'],'file') && reset ~= 1
            disp('Loading trajectory');
            data = load([dir '/traj.mat']);
            runs.traj = data.traj;
            clear data
        end

        % load eddy energy diagnostics if the file exists
        if exist([dir '/energy.mat'],'file') && reset ~= 1
            disp('Loading energy diagnostics');
            data = load([dir '/energy.mat']);
            runs.eddy.energy = data.energy;
            clear data
        end

        % set time-scale for normalization
        if isfield(runs.csflux, 'tscale')
            runs.tscale = runs.csflux.tscale;
            runs.tscaleind = runs.csflux.tscaleind;
        else
            warning(['Using eddy center based time-scale instead ' ...
                     'of flux in ' runs.name]);
            runs.tscale = runs.eddy.tscale;
            runs.tscaleind = runs.eddy.tscaleind;
        end
    end

    function [] = info(runs)
        roms_info(runs.dir);
    end

    function [] = plot_test1(runs)

        if isempty(runs.zeta), runs.read_zeta; end

        imx = vecfind(runs.rgrid.x_rho(1,:), runs.eddy.mx);
        for tt=1:size(runs.zeta,3)
            zetamat(:,tt) = runs.zeta(imx(tt),:,tt);
        end
        zy = (diff(zetamat,1,1));
        animate(zetamat'); center_colorbar;
        linex(runs.traj.tind);

        figure;
        plot(max(zy,[],1)); hold all;
        plot(-1*min(zy,[],1));
        linex(runs.traj.tind);
    end

    function [] = plot_btrq(runs)

        t0 = 5;

        figure;
        hold all;
        plot(runs.ndtime(t0:end), runs.bottom.byu(t0:end));
        plot(runs.ndtime(t0:end), runs.bottom.btrq(t0:end));
        title(runs.name);
        legend('\beta \psi', 'dH/dy * p_{bot}');
        liney(0); linex(runs.ndtime(runs.traj.tind-t0));
    end

    function [] = plot_velprofiles(runs)

        [~,~,tind] = runs.locate_resistance();

        it = [1 tind];

        ix = vecfind(runs.rgrid.x_rho(1,:), runs.eddy.mx(it));
        iy = vecfind(runs.rgrid.y_rho(:,1), ...
                     runs.eddy.my(it) - runs.eddy.vor.lmin(it)/3);

        corder_backup = get(0, 'DefaultAxesColorOrder');
        set(0, 'DefaultAxesLineStyleorder','-');
        set(0, 'DefaultAxesColorOrder', brighten(cbrewer('seq','Reds',length(it)), ...
                                                 -0.5));
        hf = figure; hold all
        t0 = 1;
        for tt=1:length(it)
            tt
            u = dc_roms_read_data(runs.dir, 'u', it(tt), {'x' ix(tt) ...
                                ix(tt); 'y' iy(tt) iy(tt)}, [], ...
                                  runs.rgrid, 'his');
            figure(hf);
            plot(abs(u./u(end)), ...
                 runs.rgrid.z_u(:,iy(tt),ix(tt)) ./ runs.eddy.Lgauss(1));
        end

        z = [-3:0.01:0];
        prof = (1-erf(abs(z)));
        hplt = plot(prof./prof(end), z, 'b-');
        linex(0); liney(runs.eddy.hcen(tind)./runs.eddy.Lgauss(t0)*-1);
        legend(hplt, '1 - erf(|z|/Lz)', 'Location', 'SouthEast');
        ylabel('z/L_z'); xlabel('U(z)/U(0)'); title(runs.name);

        runs.animate_field('u', [], it(end), 1);
        plot(runs.rgrid.x_rho(1,ix(end))/1000, ...
             runs.rgrid.y_rho(iy(end),1)/1000, 'kx');

        set(0, 'DefaultAxesColorOrder', corder_backup);
        set(0,'DefaultAxesLineStyleOrder',{'-','--','-.'});
    end

    function [] = plot_dEdt(runs)

        annostr = 'plot_dEdt';

        beta = runs.params.phys.beta;
        dt = diff(runs.eddy.t*86400);
        TE = (runs.eddy.PE + runs.eddy.KE)./abs(1);
        dEdt = bsxfun(@rdivide, diff(TE), dt');
        dhdt = bsxfun(@rdivide, diff(runs.eddy.Lgauss), dt');

        % velocity
        if runs.bathy.axis == 'y'
            v = avg1(runs.eddy.mvy);
        else
            v = avg1(runs.eddy.mvx);
        end

        tvec = avg1(runs.eddy.t*86400)./runs.eddy.turnover;

        nanmask = ~(dt == 0) &  ~isnan(v);

        figure;
        plot(tvec, avg1(TE));

        % mask out
        dEdt = dEdt(nanmask,:);
        dhdt = dhdt(nanmask);
        tvec = tvec(nanmask);
        v = v(nanmask);

        nsmooth = 20*dt(1)/runs.eddy.turnover;
        tind = 1:ceil(50*runs.eddy.turnover/86400);

        for ii=1:size(dEdt,2)
            dEdtsmth(:,ii) = smooth(dEdt(:,ii), nsmooth);
        end

        figure; insertAnnotation(annostr);
        plot(tvec, -bsxfun(@rdivide, dEdtsmth(:,1), min(dEdtsmth(:,1))));
        hold all
        %plot(tvec, -v./min(v));
        %plot(tvec, smooth(dhdt./dhdt(1), 5));
        %legend('dE/dt','cvy','dh/dt');
        ylabel('dE/dt');
        title(runs.name);
        hl =liney(0); uistack(hl,'bottom');
        xlabel('Time / turnover time');
        maximize(); pause(1);
        beautify;

        %stop

        param = (runs.eddy.Lgauss)./abs(runs.eddy.hcen - runs.eddy.Lgauss) ...
                .* runs.params.phys.beta .* runs.eddy.V;

        figure; insertAnnotation(annostr);
        hold all;
        scatter(smooth(beta * v(tind),nsmooth), smooth(dEdt(tind),nsmooth), ...
                24, tvec(tind), 'filled')
        colormap(brighten(cbrewer('seq','Reds',24), -0.75));
        title([runs.name [' | Scatter plot color coded by time | ' ...
                          'nsmooth = '] num2str(nsmooth)]);
        xlabel('v'); ylabel('dE/dt');
        beautify;

        %[c,lags] = xcorr(v(tind),dEdt(tind), 'coeff');
        %figure;
        %plot(lags,c);

    end

    function [] = plot_nof1998(runs)
        figure;
        ind = runs.csflux.tscaleind;
        H = runs.eddy.Lgauss(ind:end);
        beta = runs.params.phys.beta;
        Ldef = runs.params.eddy.Ldef;
        t = runs.eddy.t(ind:end) * 86400;
        t = t - t(1);

        Hnof = H(1)./(1 + sqrt(2) * beta .* Ldef .*t/18).^2;

        plot(t./runs.eddy.turnover, H);
        hold all
        plot(t./runs.eddy.turnover, Hnof);

    end

    function [] = plot_flierl1984(runs)
        beta = runs.params.phys.beta;
        Ldef = runs.params.eddy.Ldef;
        f0 = runs.params.phys.f0;

        cg = - beta .* Ldef^2;

        fl = f0 ./ beta ./ nanmean(runs.eddy.mvy./86.4);

        TE = runs.eddy.KE + runs.eddy.PE;
        time = runs.eddy.t * 86400;

        figure;
        plot(avg1(TE) ./ (diff(TE)./diff(time)));
        liney(fl)
    end

    function [] = plot_vvelnoise(runs)

        figure;
        insertAnnotation([runs.name '.plot_vvelnoise()']);
        % find sponge edges
        sz = size(runs.sponge);
        sx1 = find(runs.sponge(1:sz(1)/2,sz(2)/2) == 0, 1, 'first');
        sx2 = sz(1)/2 + find(runs.sponge(sz(1)/2:end,sz(2)/2) == 1, 1, ...
                             'first') - 2;

        sy2 = find(runs.sponge(sz(1)/2, :) == 1, 1, 'first') - 1;

        v1 = dc_roms_read_data(runs.dir, 'v', [], {'x' sx1 sx1; 'y' ...
                            sy2 sy2; 'z' 10 10}, [], runs.rgrid);
        v2 = dc_roms_read_data(runs.dir, 'v', [], {'x' sx2 sx2; 'y' ...
                            sy2 sy2; 'z' 10 10}, [], runs.rgrid);
        t = runs.time/86400;

        r1 = dc_roms_read_data(runs.dir, 'rho', [], {'x' sx1 sx1; 'y' ...
                            sy2 sy2; 'z' 10 10}, [], runs.rgrid);
        r2 = dc_roms_read_data(runs.dir, 'rho', [], {'x' sx2 sx2; 'y' ...
                            sy2 sy2; 'z' 10 10}, [], runs.rgrid);

        subplot(211)
        plot(t,v1,t,v2);
        linex([runs.eddy.edgtscale runs.eddy.tscale]./ ...
              runs.eddy.turnover);
        ylabel('v (m/s)');
        title(runs.name);
        legend(['(' num2str(sx1) ',' num2str(sy2) ')'], ...
               ['(' num2str(sx2) ',' num2str(sy2) ')']);
        beautify;

        subplot(212)
        plot(t,r1,t,r2);
        ylabel('\rho - 1000 (kg/m^3)');
        xlabel('Time (days)');
        beautify;
    end

    function [] = plot_cg(runs)

        vec = runs.csflux.ikefluxxt(:,300,end)';
        plot(vec);

        %        gamma = 3; beta = 2;
        %fs = morsespace(gamma, beta, length(vec));
        %w = wavetrans(vec', {gamma, beta, fs});

        m2 = 0;1./runs.rrdeep^2
        k2 = (2*pi/210/1000)^2
        l2 = (2*pi/150/1000)^2
        beta = runs.params.phys.beta;

        cgx = beta .* (k2 - (l2 + m2)) ./ (k2+l2+m2)^2
        cpx = -1 * beta ./ (k2+l2+m2)

        figure; maximize(); pause(0.1);

        pcolorcen(runs.rgrid.x_rho(1,2:end-1)'/1000, runs.eddy.t, ...
                  runs.csflux.ikefluxxt(:,:,end)');
        center_colorbar;
        hold all

        % eddy center
        plot(runs.eddy.cx/1000, runs.eddy.t);

        % long rossby wave
        c = - runs.params.phys.beta * runs.rrdeep^2
        x = xlim;
        x = [x(1):1:x(2)];
        plot(x, 150 + (x-510)./(c*86.4));
        plot(x, 200 + (x-510)./(cpx*86.4));

        legend('energy flux', 'eddy center', 'long wave speed', ...
               'estimated speed', 'Location', 'SouthWest');
    end

    function [] = predict_waterdepth(runs)
        Lz = runs.eddy.Lgauss(1);
        beta = runs.params.phys.beta;
        beta_t = runs.bathy.sl_slope * runs.params.phys.f0 ./ ...
                 Lz;

        % 1 - erf(H/Lz) = m (β/β_t) + c
        m = 1.38;
        c = 0.046;

        H = erfinv( 1 - (m*beta/beta_t + c)) * Lz
    end

    function [] = plot_sponge_enflux(runs)

        figure; maximize(); pause(0.1);
        insertAnnotation([runs.name '.plot_sponge_enflux()']);

        pcolorcen(runs.csflux.ikefluxxt(:,:,end)');
        xlabel('X (index)');
        ylabel('Time (index)');
        center_colorbar;
        beautify;

        % y (ax2) limits - between shelfbreak & edge of sponge
        sz = size(runs.sponge);
        sy1 = runs.bathy.isb;
        sy2 = find(runs.sponge(sz(1)/2, :) == 1, 1, 'first') - 1;

        ymat = repmat(runs.rgrid.y_rho(sy1:sy2,1), [1 length(runs.eddy.t)])/1000;
        tmat = repmat(runs.eddy.t,[length(runs.asflux.yvec) 1]);

        hmat = repmat(runs.rgrid.h(1,runs.asflux.iy(1):runs.asflux.iy(2))', ...
                      [1 length(runs.eddy.t)]);

        figure; maximize(); pause(0.1);
        insertAnnotation([runs.name '.plot_sponge_enflux()']);

        subplot(131);
        pcolorcen(tmat, ymat, runs.asflux.ikefluxyt(:,:,2)./hmat);
        center_colorbar; clim = caxis;
        title('Edge of sponge');
        ylabel('Y (km)');
        xlabel('Time (days)');
        beautify;
        subplot(132);
        pcolorcen(tmat, ymat, runs.asflux.ikefluxyt(:,:,5)./hmat);
        center_colorbar; caxis(clim);
        title('Middle of sponge');
        xlabel('Time (days)');
        beautify;
        subplot(133);
        pcolorcen(tmat, ymat, runs.asflux.ikefluxyt(:,:,4)./hmat);
        center_colorbar; caxis(clim);
        title('Open boundary');
        xlabel('Time (days)');
        beautify;
        spaceplots(0.05*ones([1 4]),0.05*ones([1 2]));
    end

    function [] = plot_asflux(runs)

        asindex = [1 2];

        n = length(asindex);

        % y (ax2) limits
        sy1 = runs.asflux.iy(1);
        sy2 = runs.asflux.iy(2);

        ymat = repmat(runs.rgrid.y_rho(sy1:sy2,1), [1 length(runs.eddy.t)])/1000;
        tmat = repmat(runs.eddy.t,[length(runs.asflux.yvec) 1]);
        hmat = runs.asflux.hmat;

        figure; maximize(); pause(0.1);
        insertAnnotation([runs.name '.plot_asflux()']);

        for ii=1:length(asindex)
            subplot(1,n,ii);
            pcolorcen(tmat, ymat, runs.asflux.ipefluxyt(:,:,asindex(ii))./hmat);
            center_colorbar; clim = caxis;
            title(['x = ' num2str(runs.asflux.x(asindex(ii))/1000) ...
                   ' km']);
            liney(ymat(runs.bathy.isl-sy1 + 1,1), 'slopebreak');
            liney(ymat(runs.bathy.isb-sy1 + 1,1), 'shelfbreak');
            ylabel('Y (km)');
            xlabel('Time (days)');
            beautify;
        end


        figure; maximize(); pause(0.1);
        insertAnnotation([runs.name '.plot_asflux()']);

        for ii=1:length(asindex)
            subplot(1,n,ii);
            pcolorcen(tmat, ymat, runs.asflux.eddy.ipefluxyt(:,:,asindex(ii))./hmat);
            center_colorbar; clim = caxis;
            title(['x = ' num2str(runs.asflux.x(asindex(ii))/1000) ...
                   ' km']);
            liney(ymat(runs.bathy.isl-sy1 + 1,1), 'slopebreak');
            liney(ymat(runs.bathy.isb-sy1 + 1,1), 'shelfbreak');
            ylabel('Y (km)');
            xlabel('Time (days)');
            beautify;
        end

        % spaceplots(0.05*ones([1 4]),0.05*ones([1 2]));
    end

    function [] = plot_asflux_budget(runs)
        beta = runs.params.phys.beta;
        cp = abs(min(runs.eddy.mvy) * 1000/86400);
        L =  runs.eddy.Ls(1);
        f0 = runs.params.phys.f0;

        disp(['Timescale = ' num2str(f0./beta./cp/86400/360)]);
        tind = [runs.eddy.tscaleind runs.csflux.tscaleind];

        E = runs.eddy.energy.intTE;3

        stop
        % energy budget within sponge layers
        figure;
        hold all

        plot(runs.csflux.ikeflux(:,1) + runs.csflux.ipeflux(:,1));
        plot(- runs.csflux.ikeflux(:,2) - runs.csflux.ipeflux(:,2));

        plot(runs.asflux.ikeflux(:,2) + runs.asflux.ipeflux(:,2));
        plot(- runs.asflux.ikeflux(:,3) - runs.asflux.ipeflux(:,3));
        linex(tind);
        legend('Shelfbreak', 'north', 'west', 'east');

        asindex = 3;
        figure;
        ax(1) = subplot(3,1,[1 2]);
        pcolorcen(runs.asflux.ikefluxyt(:,:,asindex));
        linex(tind);
        center_colorbar;
        colorbar('off');
        ax(2) = subplot(3,1,3);
        plot(runs.asflux.ikeflux(:,asindex));
        liney(0);
        linex(tind);
        linkaxes(ax, 'x');
    end

    function [] = read_pbot(runs)
        fname = [runs.dir '/mombudget.mat'];


        imnx = runs.spng.sx1+2; imxx = runs.spng.sx2-2;
        imny = 2; imxy = runs.spng.sy2-2;
        slbot = diff(runs.rgrid.h',1,2)./diff(runs.rgrid.y_rho',1, ...
                                              2);
        slbot = avg1(slbot(imnx:imxx, imny-1:imxy),2);

        if exist(fname, 'file')
            temp = load(fname, 'pbot');
            runs.pbot = bsxfun(@times, temp.pbot, slbot > 0.65*max(slbot(:)));
        else
            error(['No mombudget.mat file. Run bottom_torque ' ...
                   'first']);
        end
    end

    function [] = read_zeta(runs, t0, ntimes)

        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [t0 t0+1];
        else
            tind = [1 length(runs.time)];
        end

        if isempty(runs.zeta) | (max(tind) > size(runs.zeta,3)) | ...
                any(isempty(runs.zeta(:,:,tind(1):tind(2))))
            if ~runs.givenFile
                runs.zeta(:,:,tind(1):tind(2)) = ...
                    dc_roms_read_data(runs.dir,'zeta',tind,{},[], ...
                                      runs.rgrid, 'his', ...
                                      'single');
            else
                runs.zeta = (ncread(runs.out_file,'zeta'));
            end
        end
    end

    % read eddy-dye at surface and save it
    function [] = read_eddsurf(runs, t0, ntimes)
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [t0 t0+1];
        else
            tind = [1 length(runs.time)];
        end

        sz = [size(runs.rgrid.x_rho') length(runs.time)];
        if isempty(runs.eddsurf) | (t0 > size(runs.eddsurf,3))
            %| any(isnan(fillnan(runs.eddsurf(:,:,tind(1):tind(2)),0)))
            if ~runs.givenFile
                runs.eddsurf(:,:,tind(1):tind(2)) = ...
                    dc_roms_read_data(runs.dir, runs.eddname, [tind], ...
                                      {'z' runs.rgrid.N runs.rgrid.N}, ...
                                      [], runs.rgrid, 'his', 'single');
            else
                runs.eddsurf(:,:,tind) = single(squeeze(ncread(runs.out_file,runs.eddname, ...
                                                     [1 1 runs.rgrid.N ...
                                    1], [Inf Inf 1 Inf])));
            end
        end
    end

    function [] = read_csdsurf(runs, t0, ntimes)
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [t0 t0+1];
        else
            tind = [1 length(runs.time)];
        end

        if isempty(runs.csdsurf) | (t0 > size(runs.csdsurf,3)) | ...
                any(isnan(fillnan(runs.csdsurf(:,:,tind(1):tind(2)),0)))
            if ~runs.givenFile
                runs.csdsurf(:,:,tind(1):tind(2)) = ...
                    dc_roms_read_data(runs.dir, runs.csdname, [tind], {'z' ...
                                    runs.rgrid.N runs.rgrid.N}, [], runs.rgrid, ...
                                             'his');
            else
                runs.csdsurf = ...
                    single(squeeze(ncread(runs.out_file,runs.csdname, ...
                                          [1 1 runs.rgrid.N 1], ...
                                          [Inf Inf 1 Inf])));
            end
        end
    end

    function [] = read_rhosurf(runs, t0, ntimes)
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [t0 t0+1];
        else
            tind = [1 length(runs.time)];
        end

        if isempty(runs.rhosurf) | ...
                any(isempty(runs.rhosurf(:,:,tind(1):tind(2))))
            if ~runs.givenFile
                runs.rhosurf(:,:,tind(1):tind(2)) = ...
                    dc_roms_read_data(runs.dir, 'rho', [tind], ...
                                      {'z' runs.rgrid.N runs.rgrid.N}, ...
                                      [], runs.rgrid, 'his', 'single');
            else
                runs.rhosurf = single(squeeze(ncread(runs.out_file,'rho',[1 1 runs.rgrid.N ...
                                    1], [Inf Inf 1 Inf])));
            end
        end
    end

    % plot velocity / ζ sections through the eddy center
    function [] = plot_eddsec(runs, times)

        if ~exist('times', 'var') || isempty(times)
            times = linspace(runs.tscale/86400, runs.time(end)/86400, ...
                             6);
        end

        opt = 'zy';
        velname = opt(1);
        axname = opt(2); % name of axis for plot;

        figure;
        subplot(2,1,1);
        cmap = brighten(cbrewer('seq','Greys',length(times)+3),0);
        cmap = cmap(3:end,:,:); % chuck out lightest colors
        hold all;

        for ii=1:length(times)
            tind = find_approx(runs.time/86400, times(ii), 1);

            if axname == 'y'
                ref = runs.eddy.my(tind); % 0 for x-axis
                loc = runs.eddy.mx(tind);
                locax = 'x';
            else
                ref = runs.eddy.cx(tind);
                loc = runs.bathy.xsb; %runs.eddy.cy(tind);
                locax = 'y';
            end

            % free-surface ζ
            if velname == 'z'
                units = 'm';
                if axname == 'y'
                    eval(['axvec = runs.rgrid.' axname '_rho(:,1);']);
                else
                    eval(['axvec = runs.rgrid.' axname '_rho(1,:);']);
                end
                runs.read_zeta;
                vel = (dc_roms_read_data(runs, 'zeta', tind, {locax ...
                                    num2str(loc) num2str(loc); ...
                                    'z' 1 1}));
            else  % velocities
                units = 'm/s';
                if axname == 'y'
                    eval(['axvec = runs.rgrid.' axname '_' velname '(:,1);']);
                else
                    eval(['axvec = runs.rgrid.' axname '_' velname '(1,:)'';']);
                end
                vel = (dc_roms_read_data(runs, velname, tind, {locax ...
                                    num2str(loc) num2str(loc); ...
                                    'z' runs.rgrid.N ...
                                    runs.rgrid.N}));


            end

            if velname == 'u'
                [vmax, indmax] = min(vel(:));
            else
                [vmax, indmax] = max(vel(:));
            end

            if velname == 'z'
                % for ζ, the max. is the center so, this line in
                % the else condition doesn't work
                axscale = runs.eddy.vor.dia(tind)/2;
            else
                axscale = abs(axvec(indmax) - ref);
            end

            % normalized vectors to plot
            xvec = (axvec - ref)/axscale;
            vvec = vel./vmax;
            %vvec = vvec - vvec(1);
            % plot
            if velname == 'z'
                hplt = plot(xvec, vvec, '-', 'Color', cmap(ii,:));
            else
                hplt = plot(xvec, abs(vvec), '-', 'Color', cmap(ii,:));
            end

            addlegend(hplt, [num2str(times(ii)) ' | ' ...
                             num2str(vvec(runs.bathy.isb)) ' ' units]);
            if axname == 'y'
                plot(xvec(runs.bathy.isb), vvec(runs.bathy.isb), ...
                     'b*');
            end
        end

        linex([-1 0 1]); liney([0], [], 'k');
        xlim([-1 1]*6);

        limx = xlim;
        xvec = linspace(limx(1), limx(2), 60);

% $$$         a = 6;
% $$$         vel = -1*diff(exp(-abs(xvec).^(a))*(a-1)/a)./diff(xvec);
        vel = avg1(exp(-(xvec-1).^2/2));
        [vmax,indmax] = max(abs(vel));
        %plot(avg1(xvec)./xvec(indmax), vel./vmax, 'r');

        ylabel([velname ' ./ max(' velname ')'])
        xlabel([upper(axname) ' Distance from center / (radius)']);
        title([velname ' | ' runs.name]);
        beautify;

        subplot(2,1,2)
        plot(runs.csflux.time/86400, runs.csflux.west.shelf(:,1));
        limy = ylim; ylim([0 limy(2)]); linex(times);
        ylabel('Flux'); xlabel('Time (days)');
        beautify;
    end

    % read surface velocities for animate_pt & surf vorticity plot
    function [] = read_velsurf(runs, t0, ntimes)
        disp('Reading surface velocity fields...');
        start = [1 1 runs.rgrid.N 1];
        count = [Inf Inf 1 Inf];
        stride = [1 1 1 1];

        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        if ntimes == 1
            tind = [t0 t0+1];
        else
            tind = [1 length(runs.time)];
        end

        if runs.givenFile
            runs.usurf = double(squeeze(ncread(runs.out_file, ....
                'u',start,count,stride)));
        else
            runs.usurf(:,:,tind(1):tind(2)) = ...
                       dc_roms_read_data(runs.dir,'u', ...
                                         [tind],{'z' runs.rgrid.N runs.rgrid.N}, ...
                                         [],runs.rgrid, 'his', 'single');
        end
        if runs.givenFile
            runs.vsurf = double(squeeze(ncread(runs.out_file, ....
                'v',start,count,stride)));
        else
            runs.vsurf(:,:,tind(1):tind(2)) = ...
                dc_roms_read_data(runs.dir,'v', ...
                                  [tind],{'z' runs.rgrid.N runs.rgrid.N}, ...
                                  [],runs.rgrid, 'his', 'single');
        end
    end

    function [] = read_velbar(runs)
        disp('Reading depth-averaged velocity fields')

        if runs.givenFile
            runs.ubar = double(squeeze(ncread(runs.out_file, ....
                'ubar')));
        else
            runs.ubar = dc_roms_read_data(runs.dir, 'ubar', [], {}, ...
                                          [], runs.rgrid, 'his', 'single');
        end
        if runs.givenFile
            runs.vbar = double(squeeze(ncread(runs.out_file, 'vbar')));
        else
            runs.vbar = dc_roms_read_data(runs.dir, 'vbar', [], {}, ...
                                          [], runs.rgrid, 'his', 'single');
        end
    end

    % read surface velocities for animate_pt & surf vorticity plot
    function [] = read_velbot(runs)
        disp('Reading bottom velocity fields...');
        start = [1 1 1 1];
        count = [Inf Inf 1 Inf];
        stride = [1 1 1 1];

        if runs.givenFile
            runs.ubot = double(squeeze(ncread(runs.out_file, ....
                'u',start,count,stride)));
        else
            runs.ubot = dc_roms_read_data(runs.dir,'u', ...
                [],{'z' 1 1},[],runs.rgrid, ...
                                           'his', 'single');
        end
        if runs.givenFile
            runs.vbot = double(squeeze(ncread(runs.out_file, ....
                'v',start,count,stride)));
        else
            runs.vbot = dc_roms_read_data(runs.dir,'v', ...
                [],{'z' 1 1},[],runs.rgrid, ...
                                           'his', 'single');
        end
    end

   %% floats
    function [] = compare_floats(runs)
        ltransc = floats('ltrans',[runs.dir '/ltrans-compare.nc'],runs.rgrid);
        runs.roms.plot_stats;
        ltransc.plot_stats;
    end

    % create initial seed file for ltrans
    function [] = create_ltrans(runs)
        if isempty(runs.zeta)
            runs.read_zeta;
        end

        % call tools function
        ltrans_create(runs.rgrid,runs.zeta,runs.eddy, [runs.dir ...
                            '/ltrans_init.txt']);
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

   %% analysis

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

        tind = runs.eddy.tscaleind; find(runs.time == runs.eddy.trev);

        % choose plots
        trackflag = 0
        watermassflag = 0
        plumeflag = 0

        if trackflag
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
        end

        %% water mass plots
        time = runs.time/86400;
        if watermassflag
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
        end

        %% study plumes
        % offshore plume on shelf &
        % shelf plume on slope
        if plumeflag
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
        end

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
        %if isfield(runs.eddy.vor,'vol')
            figure;
            subplot(211)
            %plot(runs.time/86400,runs.eddy.vor.vol);
            %ylabel('Volume (m^3)');
            plot(runs.csflux.time/86400, ...
                 runs.csflux.west.shelf(:,1)/1e6);
            hold all
            plot(runs.csflux.time/86400, ...
                 runs.csflux.east.slope(:,1)/1e6);
            liney(0);
            ylabel('Transport (Sv)');
            legend('West - shelf water', 'East - slope water');

            subplot(212)
            hold all
            try
                plot(runs.time/86400,abs(runs.eddy.vor.zdcen));
                plot(runs.time/86400,abs(runs.eddy.vor.zcen));
            catch ME
            end
            plot(eddy.t,eddy.hcen/2);
            plot(eddy.t, ...
                runs.params.phys.f0 / sqrt(runs.params.phys.N2) * runs.eddy.Ls*2);
            plot(eddy.t,eddy.Lgauss);
            xlabel('Time (days)');
            ylabel('Z-scale (m)');
            linex(tind);
            suplabel(runs.dir,'t');
            packrows(2,1);
            legend('z-centroid','zdye-centroid', ...
                'H_{center}/2','f*dia/N','vertical (Gaussian) scale');
            %end
        end

    % make plots like dewar & hogg - looking for hydraulic jump
    function [] = tempvelsec(runs)

        t0 = 50;
        [temp, xt, yt, ~] = dc_roms_read_data(runs.dir, 'temp', [t0 Inf], ...
                                 {runs.bathy.axis runs.bathy.isb runs.bathy.isl; ...
                                  'z' 1 1}, [], runs.rgrid, 'avg', 'single');
        [asvel, xas, yas, ~] = dc_roms_read_data(runs.dir, runs.asvelname, [t0 Inf], ...
                                  {runs.bathy.axis runs.bathy.isb runs.bathy.isl; ...
                                  'z' 1 1}, [], runs.rgrid, 'avg', 'single');

        xas = xas(:,:,1); yas = yas(:,:,1);
        xt = xt(:,:,1); yt = yt(:,:,1);

        vscale = max(abs(asvel(:)));

        figure;
        tt = 1;
        tind = t0 + tt - 1;
        hgvel = pcolorcen(xas/1000, yas/1000, asvel(:,:,tt));
        caxis([-1 1]*vscale); cbfreeze; freezeColors;
        hold on;
        [~,hgtemp] = contour(xt/1000, yt/1000, temp(:,:,tt), 40, ...
                             'k');
        hbathy = runs.plot_bathy('contour', 'b');
        hgt = title(['tt = '  num2str(runs.time(t0 + tt - 1)/86400) ...
                     ' days']);
        set(gca, 'ydir', 'reverse');
        % mark eddy extents
        hglinex = linex([runs.eddy.vor.we(tind) runs.eddy.vor.ee(tind)]/1000);
        hgliney = liney([runs.eddy.vor.se(tind) runs.eddy.vor.ne(tind)]/1000);
        pause(0.01);
        for tt = 2:size(temp,3)
            tind = t0 + tt - 1;
            set(hgvel, 'CData', asvel(:,:,tt));
            set(hgtemp, 'ZData', temp(:,:,tt));
            set(hgt, 'String', ['tt = '  num2str(runs.time(t0 + tt ...
                                                           - 1)/86400) ...
                                ' days']);
            set(hglinex(1), 'XData', [1 1]*runs.eddy.vor.we(tind)/ ...
                            1000);
            set(hglinex(2), 'XData', [1 1]*runs.eddy.vor.ee(tind)/ ...
                            1000);
            set(hgliney(1), 'YData', [1 1]*runs.eddy.vor.ne(tind)/ ...
                            1000);
            set(hgliney(2), 'YData', [1 1]*runs.eddy.vor.se(tind)/ ...
                            1000);
            pause(0.01);
        end
    end

    % comparison plots
    function [] = compare_plot(runs,num)
        eddy = runs.eddy;
        % 86400 since eddy.t is already in days
        eddy.t = eddy.t./ (eddy.tscale/86400);
        ii = num;

        % line styles, markers & colors
        colors = cbrewer('qual', 'Dark2', 8);
        linestyle = {'-','--','-.','-'};
        markers = {'none','none','none','.'};

        aa = 6; bb = aa*2;
        tloc = [1:0.5:floor(max(eddy.t))];
        tind = vecfind(eddy.t, tloc);

        % plot eddy tracks
        % background velocity displacement
        if ~isfield(runs.eddy, 'bgvel')
            runs.eddy_bgflow();
        end
        displace = cumtrapz(runs.time, runs.eddy.bgvel);
        plotx = (eddy.mx - displace - eddy.mx(1))/1000;
        ploty = (eddy.my - eddy.my(1))/1000;
        figure(1)
        hold on
        %subplot(aa,2,bb);
        %subplot(aa,2,1); hold all
        %pcolorcen(runs.rgrid.xr/1000,runs.rgrid.yr/1000,runs.bathy.h);            colorbar
        xlabel('X (km)'); ylabel('Y (km)');
        he = plot(plotx, ploty, 'Color',colors(ii,:),'LineWidth',2);
        hold all;
        addlegend(he, runs.name, 'NorthEast');
        try
            plot(plotx(tind), ploty(tind),'*', ...
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
        %{
        if isfield(runs.csflux,'west')
            ftime = runs.csflux.time/eddy.tscale;
            figure(4);
            subplot(4,1,1);
            hold on;
            plot(ftime, runs.csflux.west.shelf(:,1), 'Color', colors(ii,:));
            ylabel('Shelf water flux - sb');
            title('West');
            xlim(limx);

            subplot(4,1,2);
            hold on;
            plot(ftime, runs.csflux.west.slope(:,1), 'Color', colors(ii,:));
            ylabel('Slope water flux - sb');
            xlim(limx);

            try
                subplot(4,1,3);
                hold on;
                plot(ftime, runs.csflux.west.pv(:,1), 'Color', colors(ii,:));
                ylabel('PV flux');
                xlim(limx);
            catch ME
            end

            try
                subplot(4,1,4);
                hold on;
                plot(ftime, runs.csflux.west.rv(:,1), 'Color', colors(ii,:));
                ylabel('RV flux');
                xlim(limx);
            catch ME
            end

            figure(5);
            subplot(4,1,1);
            hold on;
            plot(ftime, runs.csflux.east.shelf(:,1), 'Color', colors(ii,:));
            ylabel('Shelf water flux - sb');
            title('East');
            xlim(limx);

            subplot(4,1,2);
            hold on;
            plot(ftime, runs.csflux.east.slope(:,1), 'Color', colors(ii,:));
            ylabel('Slope water flux - sb');
            xlim(limx);

            try
                subplot(4,1,3);
                hold on;
                plot(ftime, runs.csflux.east.pv(:,1), 'Color', colors(ii,:));
                ylabel('PV flux');
                xlim(limx);

                subplot(4,1,4);
                hold on;
                plot(ftime, runs.csflux.east.rv(:,1), 'Color', colors(ii,:));
                ylabel('RV flux');
                xlim(limx);
            catch ME
            end
        end
        %}
        time = eddy.t;

        %% plot water masses
        if isfield(runs.water, 'off')
            % normalize volumes by initial eddy volume
            if isfield(runs.eddy, 'vol')
                evol0 = 1;runs.eddy.vol(runs.eddy.tscaleind);
            else
                evol0 = 1;
            end
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

        % background flow velocity estimates
        %{
        figure(6)
        subplot(211); hold on;
        hbg = plot(time, runs.eddy.bgvel, 'Color', colors(ii,:));
        ylabel('mean(vel. at x=eddy center)');
        addlegend(hbg, runs.name, 'NorthWest');
        subplot(212); hold on;
        if runs.bathy.axis == 'y'
            plot(time, squeeze(mean(runs.ubar(3,:,:),2)), 'Color', ...
                 colors(ii,:));
        else
            plot(time, squeeze(mean(runs.vbar(:,3,:),1)), 'Color', ...
                 colors(ii,:));
        end
        xlabel('Time');
        ylabel('mean(inflow 2d vel)');
        %}

        % jet diagnostics
        if isfield('jet', runs)
            jtime = runs.time/86400;
            figure(7);
            subplot(3,2,1)
            plot(jtime, runs.jet.vscale, 'Color', colors(ii,:));
            ylabel('Max. velocity');
            subplot(3,2,2)
            plot(jtime, runs.jet.zscale, 'Color', colors(ii,:));
            hold on
            plot(jtime, -1*runs.jet.h, 'Color', colors(ii,:), 'LineStyle', ...
                 '-.');
            ylabel('z-loc of max vel');
            subplot(323)
            plot(jtime, runs.jet.bc, 'Color', colors(ii,:));
            liney(0);
            ylabel('Baroclinicity measure');
            subplot(324)
            plot(time, runs.jet.yscale, 'Color', colors(ii,:));
            ylabel('y-loc of max. vel');
        end

        figure(8);
        set(gcf, 'Renderer', 'painters');
        subplot(2,1,1)
        hold on
        plot(runs.csflux.time/eddy.tscale, ...
             runs.csflux.west.shelf(:,1)/1e6, 'Color', colors(ii,:));
        plot(runs.csflux.time/eddy.tscale, ...
             runs.csflux.east.slope(:,1)/1e6, 'Color', colors(ii,:), ...
             'LineStyle', linestyle{2});
        liney(0);
        ylabel('Transport (Sv)');
        legend('West - shelf water', 'East - slope water');

        subplot(2,1,2)
        hold on
        hline = plot(eddy.t,eddy.hcen, 'Color', colors(ii,:));
        xlabel('Time (days)');
        ylabel('center-isobath');
        %linex(tind);
        suplabel(runs.dir,'t');
        %        packrows(2,1);
        addlegend(hline, runs.name);

        % shelf water envelope
        if isfield(runs.csflux.west.shelfwater, 'envelope')
            normtrans = sum(runs.csflux.west.shelfwater.itrans);

            figure(9)
            subplot(2,1,1)
            hold on
            hline = plot(runs.csflux.time/86400, ...
                         runs.csflux.west.shelfwater.envelope/ runs.rrshelf, ...
                         'Color', colors(ii,:));
            addlegend(hline, runs.name);
            xlabel('Time');
            ylabel({'Location of water parcel farthest from shelfbreak' ...
                    'in terms of shelfbreak rossby radius'})

            subplot(2,1,2)
            hold on;
            plot(runs.csflux.west.shelfwater.bins/runs.rrshelf, ...
                 runs.csflux.west.shelfwater.itrans./normtrans, 'color', colors(ii,:));
            ylabel('Total volume transported');
            xlabel('Bin = location / RR_{shelf} ');
        end

        % shelf water vorticity budget
        if isfield(runs, 'vorbudget')
            figure(10)
            subplot(2,1,1)
            hold all
            hline = plot(runs.csflux.time/86400, runs.csflux.west.shelf, ...
                         'Color', colors(ii,:));
            addlegend(hline, runs.name);

            subplot(2,1,2)
            hold all
            hline = plot(runs.vorbudget.time/86400, ...
                         runs.vorbudget.shelf.str, 'Color', ...
                         colors(ii,:));
        end
    end

    % calculate surface vorticity field
    function [] = calc_vorsurf(runs)
        if isempty(runs.usurf) || isempty(runs.vsurf)
            runs.read_velsurf;
        end

        if isempty(runs.vorsurf)
            vx =  bsxfun(@rdivide,diff(runs.vsurf,1,1), ...
                         diff(runs.rgrid.x_v',1,1));

            uy = bsxfun(@rdivide,diff(runs.usurf,1,2), ...
                        diff(runs.rgrid.y_u',1,2));

            runs.vorsurf = vx - uy;

            runs.rgrid.xvor = avg1(avg1(runs.rgrid.xr,1),2);
            runs.rgrid.yvor = avg1(avg1(runs.rgrid.yr,1),2);
        end
    end

    % check time vectors
    function [] = check_time(runs)
        figure; hold all;
        try
            plot(runs.time/86400);
            plot(runs.eddy.t);
            plot(runs.csflux.time/86400);
        catch ME
        end
        legend('time', 'eddy.time', 'csflux.time');
    end

    % calculate geostrophically balanced barotropic velocities
    function [] = calc_ubarg(runs)
        runs.ubarg = -1 * 9.81 .* bsxfun(@rdivide,diff(runs.zeta,1,2), ...
                            avg1(runs.rgrid.f',2).*diff(runs.rgrid.yr,1,2));

        runs.vbarg =      9.81 .* bsxfun(@rdivide,diff(runs.zeta,1,1), ...
                            avg1(runs.rgrid.f',1).*diff(runs.rgrid.xr,1,1));
    end

    % let's try to estimate background flow acting on eddy
    function [] = eddy_bgflow(runs)
        if runs.bathy.axis == 'y'
            if isempty(runs.ubar)
                runs.ubar = dc_roms_read_data(runs.dir, 'ubar', [], ...
                                              {},  [], runs.rgrid);
            end
            bg = runs.ubar;
            cind = vecfind(runs.rgrid.x_rho(1,:), runs.eddy.cx);
            edgeind = vecfind(runs.rgrid.y_rho(:,1), runs.eddy.ne);
            for ii = 1:size(bg,3)
                runs.eddy.bgvel(ii) = mean(runs.ubar(cind(ii), ...
                                                     edgeind(ii):end, ii));
            end
        else
            if isempty(runs.vbar)
                runs.vbar = dc_roms_read_data(runs.dir, 'vbar', [], ...
                                              {}, [], runs.rgrid);
            end
            bg = runs.vbar;
            cind = vecfind(runs.rgrid.y_rho(:,1), runs.eddy.cy);
            edgeind = vecfind(runs.rgrid.y_rho(:,1), runs.eddy.ee);
            for ii = 1:size(bg,3)
                runs.eddy.bgvel(ii) = mean(runs.vbar(edgeind(ii):end, ...
                                                     cind(ii), ii));
            end
        end
    end

    % check eddy vertical scale estimations
    function [] = plot_eddy_vscale(runs)

        c = hypot(runs.eddy.cvx, runs.eddy.cvy) / 86.4;
        c = nanmean(c(1:50));

        % read in initial velocity field
        u0 = ncread(runs.out_file, 'u', [1 1 1 1], [Inf Inf Inf ...
                            1]);
        v0 = ncread(runs.out_file, 'v', [1 1 1 1], [Inf Inf Inf ...
                            1]);
        U = hypot(avg1(u0(:,2:end-1,:), 1), avg1(v0(2:end-1,:,:), ...
                                                 2));

        % volume of eddy that satisfies U/c criterion
        dV = bsxfun(@times, runs.rgrid.dV(2:end-1, 2:end-1,:) ...
                    .* (U > c), runs.eddy.vormask(:,:,1));
        runs.eddy.Ucvol = nansum(dV(:));

        if ~isfield(runs.eddy, 'zT') || isempty(runs.eddy.zT)
            for ii=1:size(runs.eddy.T, 1)
                ix = vecfind(runs.rgrid.x_rho(1,:), ...
                             runs.eddy.vor.cx(ii));
                iy = vecfind(runs.rgrid.y_rho(:,1), ...
                             runs.eddy.vor.cy(ii));
                runs.eddy.zT(ii,:) = squeeze(runs.rgrid.z_r(:, iy, ix))';
            end
            runs.eddy.tmat = repmat(runs.time', [1 size(runs.eddy.T, ...
                                                     2)]);
            eddy = runs.eddy;
            save([runs.dir '/eddytrack.mat'], 'eddy');
        end

        figure;
        if isfield(runs.eddy, 'dyecen')
            subplot(211)
            pcolorcen(runs.eddy.tmat, runs.eddy.zT, ...
                     runs.eddy.dyecen);
            %   colormap(flipud(colormap('bone')));
            caxis([0 1]);
            colorbar;
            hold all
            plot(runs.time/86400, -1*runs.eddy.Lz2, 'c');
            plot(runs.time/86400, -1*runs.eddy.Lgauss, 'm');

            tcen = find_approx(runs.eddy.my, runs.bathy.xsl, ...
                               1);
            tse = find_approx(runs.eddy.se, runs.bathy.xsl, 1);
            linex([tse tcen], [], 'r');
            title(['Eddy dye profiles | ' runs.name]);
            subplot(212)
        end
        pcolorcen(runs.eddy.tmat, runs.eddy.zT, ...
                 runs.eddy.T./max(runs.eddy.T(1,:)));
        %        colormap(flipud(colormap('bone')));
        caxis([0 1]);
        hold all
        plot(runs.time/86400, -1*runs.eddy.Lz2, 'c');
        plot(runs.time/86400, -1*runs.eddy.Lgauss, 'm');
        ylabel('Z (m)'); xlabel('Time (days)');
        title('Scaled temp anomaly at eddy center');
        colorbar;
        legend('Scaled temp anomaly', 'sine fit', 'Gaussian fit', 'Location', ...
               'SouthEast');
        contour(runs.eddy.tmat, runs.eddy.zT, runs.eddy.T, [0], ...
                'LineWidth', 2,'Color', 'k');
    end

    % filter floats that cross shelfbreak
    function [] = filter_cross_sb(runs, name)

        tic;
        if ~exist('name', 'var') || isempty(name), name = 'tracpy'; end

        winds = []; einds = []; reenter_inds = [];

        % assign float structure
        eval(['flt = runs.' name ';']);
        if isempty(flt); return; end

        xsb = runs.bathy.xsb;

        % floats that start shoreward of shelfbreak and end
        % offshore of shelfbreak
        %ind = find(flt.y(end,:) > xsb & flt.init(:,2)' < xsb);
        ind = find(flt.init(:,2)' < xsb);
        x1 = flt.x(:,ind); y1 = flt.y(:,ind);

        % how many floats re-enter shelf
        reenter = 0;

        % loop through filtered floats
        for ff=1:size(y1,2)
            % now to check when they actually cross the shelfbreak
            tind1 = find(y1(:,ff) < xsb, 1, 'last');
            tind2 = find(y1(:,ff) > xsb, 1, 'first');

            % oops, this one didn't cross
            if isempty(tind2), continue; end

            % find appropriate time index in eddy time series
            tfloat = flt.time(tind2);
            teddy = find_approx(runs.eddy.t * 86400, tfloat, 1);

            % make sure they end off the shelf
            if y1(end, ff) > xsb
                % these crossed to the west off the eddy
                if x1(tind2, ff) < runs.eddy.vor.ee(teddy)
                    winds = [winds; ind(ff)];
                else
                    einds = [einds; ind(ff)];
                end
            else
                % these didn't end off the shelf

                % These have crossed the shelfbreak to the WEST of
                % the eddy and then crossed back in
                if x1(tind2, ff) < runs.eddy.vor.ee(teddy)
                    reenter = reenter + 1;
                    reenter_inds = [reenter_inds; ff];
                end
            end
        end

        eval(['runs.', name, '.winds = winds;']);
        eval(['runs.', name, '.einds = einds;']);
        eval(['runs.', name, '.reenter = reenter;']);
        eval(['runs.', name, '.reenter_inds = reenter_inds;']);
        toc;

    end

    % filter floats that start same time, place as ROMS deployment
    function [] = filter_floats_start_roms(runs, type)

        type = 'tracpy';

        flt = runs.(type);
        roms = runs.roms;

        roms_inds = [];
        for ii=1:size(roms.init,1)
            ff = find(abs(flt.init(:,1) - roms.init(ii,1)) < 1e1 & ...
                      abs(flt.init(:,2) - roms.init(ii,2)) < 1e1 & ...
                      abs(flt.init(:,3) - roms.init(ii,3)) < 1e-3 & ...
                      abs(flt.init(:,4) - roms.init(ii,4)) < 0.5*86400);
            roms_inds = [roms_inds; ff];
        end

    end

    % check w-noise level
    function [] = wnoise(runs)
        w = dc_roms_read_data(runs.dir, 'w', [], {'y' runs.bathy.isl-5 ...
                            runs.bathy.isl+5}, [], runs.rgrid, ...
                              'avg');

        % reduce
        mask = ~permute(repnan(runs.eddy.mask(:,runs.bathy.isl-5: ...
                                      runs.bathy.isl+5,:), 0), [1 2 4 3]);
        sz = size(w);
        runs.wmetric = max(abs(reshape(bsxfun(@times, w(2:end-1,:,:,:), mask), ...
                             [(sz(1)-2)*sz(2)*sz(3) sz(4)])), [], 1);
        plot(runs.wmetric);

    end

    % calculate ctw dispersion relation
    function [] = ctwdisprel(runs)

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

    %% deprecated functions

    % study vorticity export onto shelf
    % now in fluxes
    function [] = deprecated_vorexport(runs)
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

    function [] = deprecated_transport(runs)
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

    function [] = plot_shelfvelprofiles(runs, tt)
    % read cross-shore velocities on shelf
        cind = find_approx(runs.rgrid.x_rho(1,:), runs.eddy.vor.cx(tt) ...
                           , 1);
        eind = find_approx(runs.rgrid.x_rho(1,:), runs.eddy.vor.ee(tt) ...
                           , 1);
        wind = find_approx(runs.rgrid.x_rho(1,:), runs.eddy.vor.we(tt) ...
                           , 1);

        v = dc_roms_read_data(runs.dir, 'v', tt, {'y' 1 runs.bathy.isb+3}, ...
                              [], runs.rgrid, 'his', 'single');
        %%
        ind = [5:10:runs.bathy.isb-1 runs.bathy.isb];

        scale = 5e3;

        figure;
        subplot(1,3,1)
        hold on
        for ii=1:length(ind)
            plot(runs.rgrid.y_rho(ind(ii),1)/1000  + squeeze(v(wind, ...
                                                              ind(ii),:))*scale, ...
                 runs.rgrid.z_r(:,ind(ii),1), 'k');
        end
        linex(runs.bathy.xsb/1000);
        title(['West | day ' num2str(runs.time(tt)/86400)]);

        subplot(1,3,2);
        hold on
        for ii=1:length(ind)
            plot(runs.rgrid.y_rho(ind(ii),1)/1000  + squeeze(v(cind, ...
                                                              ind(ii),:))*scale, ...
                 runs.rgrid.z_r(:,ind(ii),1), 'k');
        end
        linex(runs.bathy.xsb/1000);
        title(['Center | day ' num2str(runs.time(tt)/86400)]);

        subplot(1,3,3);
        hold on
        for ii=1:length(ind)
            plot(runs.rgrid.y_rho(ind(ii),1)/1000  + squeeze(v(eind, ...
                                                              ind(ii),:))*scale, ...
                 runs.rgrid.z_r(:,ind(ii),1), 'k');
        end
        linex(runs.bathy.xsb/1000);
        title(['East | day ' num2str(runs.time(tt)/86400)]);

    end

    function [] = slope_parameter(runs)
        u = dc_roms_read_data(runs.dir, 'u', [], {'z' 1 2}, [], ...
                              runs.rgrid, 'his', 'single');
        v = dc_roms_read_data(runs.dir, 'v', [], {'z' 1 2}, [], ...
                              runs.rgrid, 'his', 'single');

        zu = permute(runs.rgrid.z_u, [3 2 1]);
        zv = permute(runs.rgrid.z_v, [3 2 1]);

        uz = bsxfun(@times, diff(u,1,3), diff(zu,1,3));
        vz = bsxtun(@times, diff(v,1,3), diff(zv,1,3));

        alpha = runs.params.bathy.sl_slope;

        stat u
        stat uz
        stat v
        stat vz
    end

    function [] = deprecated_distrib_csz(runs)
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

    function [] = disprel(runs)

        beta = runs.params.phys.beta;
        Ldef = sqrt(runs.params.phys.N2) * runs.bathy.hsb / ...
               runs.params.phys.f0;

        k = 2*pi./[0:0.05:50]/1000;
        figure;
        hold all
        for n = 0:10
            c = - beta ./ (k.^2 + (n*pi/Ldef)^2);
            hgplt = plot(k, c);
            addlegend(hgplt, num2str(n));
        end
        tind = find_approx(runs.eddy.t*86400 / runs.eddy.tscale, ...
                           1.5, 1);
        liney(mean(runs.eddy.mvx(tind:end)) * 1000/86400);
    end

    function [] = csvel_hov(runs, iz, loc)
        if ~exist('loc', 'var') || isempty(loc)
            loc = [-30 -20 -10 0] * 1000 + runs.bathy.xsb;
            locu = [200 250 300 350]*1000;
        end

        flag_v = 1; % csvel
        flag_z = 0; % zeta
        flag_u = 0; % csvel y-z

        if ~exist('iz', 'var') || isempty(iz)
            iz = runs.rgrid.N;
        end

        if iz == runs.rgrid.N
            if isempty(runs.vsurf)
                runs.read_velsurf;
            end
            u = runs.usurf;
            v = runs.vsurf;
            titlestr = 'surface';
        else
            if isempty(runs.vbot)
                runs.read_velbot;
            end
            u = runs.ubot;
            v = runs.vbot;
            titlestr = 'bottom';
        end

        % step topography phase speed estimate
        cp = runs.params.phys.f0 * runs.bathy.xsb * ...
             (runs.bathy.hsl - runs.bathy.hsb)./runs.bathy.hsl;

        if flag_v
            figure;
            for ii=1:length(loc)
                ax(ii) = subplot(2,2,ii);

                ind = find_approx(runs.rgrid.y_v(:,1), loc(ii), 1);

                xmat = repmat(runs.rgrid.x_v(1,:)', [1 size(runs.vsurf, 3)]);
                tmat = repmat(runs.time, [size(runs.vsurf, 1) 1]);
                pcolorcen(xmat/1000, tmat./runs.eddy.tscale, squeeze(runs.vsurf(:,ind,:)));
                colorbar; center_colorbar;
                xlabel('X (km)');
                ylabel('Time (non-dimensional)');
                title([titlestr ' v (m/s) at y = ' num2str(loc(ii)) ' km']);

                hold on
                plot(runs.eddy.vor.cx/1000, runs.eddy.t*86400 / runs.eddy.tscale);
                plot(runs.eddy.vor.ee/1000, runs.eddy.t*86400 / runs.eddy.tscale);
                plot(runs.eddy.vor.we/1000, runs.eddy.t*86400 / ...
                     runs.eddy.tscale);

                limx = xlim;
                xvec = limx(1):10:limx(2);
                tvec = 0.5 + 1./0.02 .* (xvec * 1000)./runs.eddy.tscale;
                plot(xvec, tvec);

                beautify;
            end

            suplabel(runs.name, 't');
            linkaxes(ax, 'xy');
        end

        %%%%%%%%%%%%%% ζ
        % along-shelfbreak pressure gradient at shelfbreak
        if flag_z
            %zx = bsxfun(@rdivide, squeeze(diff(runs.zeta(:, runs.bathy.isb, :), ...
            %                          1, 1)), diff(runs.rgrid.x_rho(1,:)', ...
            %1, 1));
            loc2 = [-20 0 10 20] + runs.bathy.isb;
            runs.read_zeta;
            figure;
            for ii=1:length(loc2)
                ax(ii) = subplot(2,2,ii);
                pcolorcen(xmat/1000, tmat/runs.eddy.tscale, ...
                          squeeze(runs.zeta(:, loc2(ii), :)));
                colorbar; center_colorbar;
                hold on
                plot(runs.eddy.vor.cx/1000, runs.eddy.t*86400 / runs.eddy.tscale);
                plot(runs.eddy.vor.ee/1000, runs.eddy.t*86400 / runs.eddy.tscale);
                plot(runs.eddy.vor.we/1000, runs.eddy.t*86400 / ...
                     runs.eddy.tscale);
                xlabel('X (km)');
                ylabel('Time (non-dimensional)');
                title(['zeta (m) at y = ' ...
                       num2str(runs.rgrid.y_rho(loc2(ii),1))]);
                beautify;
            end

            %suplabel(runs.name, 't');
            linkaxes(ax, 'xy');
        end

        %%%%%%%%%%%%%% u
        if flag_u
            figure;
            for ii=1:length(locu)
                ax(ii) = subplot(2,2,ii);

                ind = find_approx(runs.rgrid.x_v(1,:), locu(ii), 1);

                ymat = repmat(runs.rgrid.y_u(:,1), [1 size(runs.usurf, 3)]);
                tmat = repmat(runs.time, [size(runs.usurf, 2) 1]);
                contourf(tmat./86400, ymat./1000, squeeze(runs.usurf(ind,:,:)));
                colorbar; center_colorbar;
                ylabel('Y (km)');
                xlabel('Time (non-dimensional)');
                ylim([0 runs.bathy.xsb/1000]);
                title([titlestr ' u (m/s) at x = ' num2str(locu(ii)/1000) ' km']);

                beautify;
            end
        end
    end

    %% animation functions

    function [] = animate_sbvel(runs, t0)

        if ~exist('t0', 'var'), t0 = 1; end

        ftype = 'his';

        usb = dc_roms_read_data(runs.dir, 'u', ...
                                [t0 Inf], {runs.bathy.axis runs.bathy.isb runs.bathy.isb; ...
                            'z' runs.rgrid.N runs.rgrid.N}, ...
                                [], runs.rgrid, ftype, 'single');

        vsb = squeeze(avg1(dc_roms_read_data(runs.dir, 'v', ...
                                             [t0 Inf], {runs.bathy.axis runs.bathy.isb-1 runs.bathy.isb; ...
                            'z' runs.rgrid.N runs.rgrid.N}, ...
                                             [], runs.rgrid, ftype, ...
                                             'single'), 2));

        csd = dc_roms_read_data(runs.dir, runs.csdname, ...
                                [t0 Inf], {runs.bathy.axis runs.bathy.isb runs.bathy.isb; ...
                            'z' runs.rgrid.N runs.rgrid.N}, ...
                                [], runs.rgrid, ftype, 'single');


        usb = avg1(usb,1);
        vsb = vsb(2:end-1, :);
        shelfmask = csd(2:end-1, :) < runs.bathy.xsb;
        ii=1;

        figure;
        %quiver(runs.rgrid.x_r(1,2:end-1),
    end

    function [] = animate_zeta(runs, hax, t0, ntimes)
        if ~exist('hax', 'var'), hax = []; end
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = length(runs.time); end

        runs.animate_field('zeta', hax, t0, ntimes);
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

    function [] = animate_vorsurf(runs, hax, t0, ntimes)
        if isempty(runs.vorsurf)
            runs.calc_vorsurf();
        end

        if ~exist('hax', 'var'), hax = []; end
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = length(runs.time); end

        runs.video_init('surfvorcsd');

        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        % read in dye values
        if isempty(runs.csdsurf)
            runs.csdsurf = dc_roms_read_data(runs.dir, runs.csdname, [], {'z' ...
                                runs.rgrid.N runs.rgrid.N}, [], runs.rgrid, ...
                                             'his');
        end

        csdlevels = runs.bathy.xsb + [0 -20] * 1000;

        tt = t0;
        if isempty(hax), figure(); else axis(hax); end
        maximize();
        vorsurf = bsxfun(@rdivide, runs.vorsurf, ...
                         avg1(avg1(runs.rgrid.f', 1), 2));
        if ntimes ~= 1
            vormax = max(abs(vorsurf(:)));
        else
            vormax = max(max(abs(vorsurf(ix,iy,t0))));
        end

        levels = linspace(-vormax,vormax,20);
        [hh] = pcolor(runs.rgrid.xvor(ix,iy)/1000, ...
                      runs.rgrid.yvor(ix,iy)/1000, ...
                      vorsurf(ix,iy,tt));
        caxis([-1 1] * vormax); hcbar = colorbar; shading interp;
        hcbar.Label.String = '\zeta/f';
        hold on
        % bathy contour
        runs.plot_bathy('contour', [1 1 1]*0.7);

        % eddy contours
        hedd = runs.plot_eddy_contour('contour', tt);
        % hssh = runs.plot_eddy_sshcontour('contour', tt);

        % dye-contours (current)
        %[~,hcsd] = contour(runs.rgrid.x_rho'/1000, runs.rgrid.y_rho'/1000, ...
        %                   runs.csdsurf(:,:,tt), csdlevels, 'Color', [1 1 1]*0.5, ...
        %                   'LineWidth', 2);
        % dye contours (initial)
        % contour(runs.rgrid.x_rho'/1000, runs.rgrid.y_rho'/1000, ...
        %        runs.csdsurf(:,:,tt), csdlevels(1), 'Color', [1 1 1]*0.9, ...
        %        'LineWidth', 2)
        xlabel('X (km)'); ylabel('Y (km)');
        axis image;
        ht = title(['Surface vorticity @ t = ' ...
                    num2str(runs.time(tt)/86400) ' days']);
        ax = gca;
        htext = text(0.05,0.9, ...
                     ['t = ' num2str(runs.time(tt)/86400) ' days'], ...
                     'Units', 'normalized');
        colormap(flipud(cbrewer('div','RdBu',20)));
        beautify([18 18 20]);

        for tt = 2:3:ntimes
            set(hh,'ZData', double(vorsurf(ix,iy,tt)));
            %set(hcsd, 'ZData', double(runs.csdsurf(ix,iy,tt)));
            %set(h0vor, 'ZData', double(runs.vorsurf(:,:,tt)));
            shading flat;
            runs.update_eddy_contour(hedd, tt);
            % runs.update_eddy_sshcontour(hssh, tt);
            set(ht,'String',['Surface vorticity | t = ' ...
                             num2str(runs.time(tt)/86400) ' ' ...
                             'days']);
            set(htext,'String', ['t = ' num2str(runs.time(tt)/86400) ...
                                ' days']);
            runs.video_update();
            pause(0.05);
        end

        runs.video_write();
    end

    function [] = create_floats_init_file(runs, type)
        if ~strcmpi(type, 'roms')
            eval(['init = runs.' type '.init;']);
        else
            init = runs.roms.actual_init;
        end

        fname = [runs.dir '/' type '_init.mat'];
        disp(['Writing to ' , fname]);
        save(fname, 'init');
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

    function [] = plot_streamerwidth(runs)
        cxi = runs.eddy.vor.ee;
        westmask = bsxfun(@times, ...
                          bsxfun(@lt, runs.eddy.xr(:,1), cxi), ...
                          ~runs.sponge(2:end-1,1));
        eastmask = bsxfun(@times, 1 - westmask, ...
                          ~runs.sponge(2:end-1,1));

        shelfxt = runs.csflux.shelfxt;

        for tt=runs.eddy.tscaleind:size(shelfxt,2)
            clf;
            plot((runs.eddy.xr(:,1) - runs.eddy.vor.ee(tt))/1000, shelfxt(: ...
                                                              ,tt));
            ylim([min(shelfxt(:)) max(shelfxt(:))]);
            xlim([-1 1]*200);
            lmin = runs.eddy.vor.lmin(tt)/1000;
            lmaj = runs.eddy.vor.lmaj(tt)/1000;
            linex([0 -lmin -lmaj], {'east edge'; 'minor axis'; ...
                                'major axis'}, []);
            title([runs.name ' | ' num2str(tt)]);
            pause(0.05);
        end
    end

    function [] = plot_shelfvorbudget(runs)

        if ~isfield(runs.vorbudget, 'time')
            error('Vorbudget terms haven''t been calculated');
        end

        time = runs.vorbudget.time/86400;

        figure;
        subplot(2,1,1)
        [ax,h1,h2] = plotyy(time, runs.vorbudget.shelf.rv, runs.csflux.time/86400, ...
                            runs.csflux.west.shelf/1e6);
        set(ax(2), 'XTick', []);
        ylabel(ax(2), 'Cross-shelfbreak Transport (Sv)')

        beautify;
        limx = xlim;

        xlabel('Time (days)');
        ylabel('Volume averaged relative vorticity (shelf water)');
        liney(0, [], [1 1 1]*0.5);

        subplot(2,1,2)
        hold all
        plot(time, runs.vorbudget.shelf.str, 'Color', [0.68 0.85 ...
                            0.90]);
        plot(time, runs.vorbudget.shelf.hadv + runs.vorbudget.shelf.vadv, ...
             time, runs.vorbudget.shelf.tilt, ...
             time, runs.vorbudget.shelf.bfric, ...
             time, runs.vorbudget.shelf.beta);
        plot(time, runs.vorbudget.shelf.hadv, 'Color',[1 1 1]*0.6);
        plot(time, runs.vorbudget.shelf.vadv, 'Color',[1 1 1]*0.6);
        plot(time, smooth(runs.vorbudget.shelf.str, 3), 'Color', [0 0 1]);
        xlim(limx);
        liney(0, [], [1 1 1]*0.5);
        xlabel('Time (days)');
        ylabel('sec^{-2}');
        legend('str', 'adv', 'tilt', 'bfric', 'beta', 'hadv', ...
               'vadv', 'Location', 'NorthWest');
        beautify
    end

    function [] = animate_vorbudget(runs,tind, plotflag)

        vorbudgetstart = tic;

        runs.vorbudget = [];

        if ~exist('tind','var')
            tind = 1;
        end
%             if ~exist([runs.dir '/ocean_vor.nc'],'file')
%                 dc_roms_vorticity(runs.dir,tind,'ocean_vor.nc');
%             end

        runs.video_init('vor');

        if ~exist('plotflag', 'var')
            plotflag = 1;
        end

        debug = 0;

        %%
        zmin = -1 * runs.bathy.hsb;min(runs.rgrid.z_r(:));
        zmax = max(runs.rgrid.z_r(:,end,end));
        zwnew = unique([linspace(zmin, -1*runs.bathy.hsb, 70) ...
                        linspace(-1*runs.bathy.hsb, zmax-0.01, 36)]');
        zrnew = avg1(zwnew);
        % for integrating quantities later
        zint = avg1(zrnew);

        % prepare grids for differentiation
        xvor = avg1(avg1(runs.rgrid.xr,1),2);
        yvor = avg1(avg1(runs.rgrid.yr,1),2);
        N = runs.rgrid.N;
        Nnew = length(zrnew);

        % setup grid
        [sx sy] = size(runs.rgrid.x_rho');
        gridu.xmat = repmat(runs.rgrid.x_u',[1 1 Nnew]);
        gridu.ymat = repmat(runs.rgrid.y_u',[1 1 Nnew]);
        gridu.zmat = permute(runs.rgrid.z_u,[3 2 1]);
        gridu.znew = repmat(permute(zrnew,[3 2 1]),[sx-1 sy 1]);
        %gridu.s = runs.rgrid.s_rho;
        %gridu.zw = runs.rgrid.z_w;
        %gridu.s_w = runs.rgrid.s_w;

        gridv.xmat = repmat(runs.rgrid.x_v',[1 1 Nnew]);
        gridv.ymat = repmat(runs.rgrid.y_v',[1 1 Nnew]);
        gridv.zmat = permute(runs.rgrid.z_v,[3 2 1]);
        gridv.znew = repmat(permute(zrnew,[3 2 1]),[sx sy-1 1]);
        %gridv.s = runs.rgrid.s_rho;
        %gridv.zw = runs.rgrid.z_w;
        %gridv.s_w = runs.rgrid.s_w;

        gridr.xmat = repmat(runs.rgrid.x_rho',[1 1 Nnew]);
        gridr.ymat = repmat(runs.rgrid.y_rho',[1 1 Nnew]);
        gridr.zmat = permute(runs.rgrid.z_r,[3 2 1]);
        gridr.znew = repmat(permute(zrnew,[3 2 1]),[sx sy 1]);
        %gridr.s = runs.rgrid.s_rho;
        %gridr.zw = runs.rgrid.z_r;
        %gridr.s_w = runs.rgrid.s_w;

        gridw.xmat = repmat(runs.rgrid.x_rho',[1 1 Nnew+1]);
        gridw.ymat = repmat(runs.rgrid.y_rho',[1 1 Nnew+1]);
        gridw.zmat = permute(runs.rgrid.z_w,[3 2 1]);
        gridw.znew = repmat(permute(zwnew,[3 2 1]),[sx sy 1]);
        %gridw.s = runs.rgrid.s_w;
        %gridw.zw = runs.rgrid.z_w;
        %gridw.s_w = runs.rgrid.s_w;

        gridrv.xmat = repmat(xvor,[1 1 Nnew]);
        gridrv.ymat = repmat(yvor,[1 1 Nnew]);
        gridrv.zmat = avg1(avg1(avg1(permute(runs.rgrid.z_r,[3 2 1]),1),2),3);
        gridrv.znew = repmat(permute(zrnew,[3 2 1]),[sx-1 sy-1 1]);
        gridrv.s = avg1(runs.rgrid.s_rho);
        gridrv.zw = avg1(avg1(avg1(permute(runs.rgrid.z_w,[3 2 1]),1),2),3);
        gridrv.s_w = avg1(runs.rgrid.s_w);

        % for depth integration - banas code
        h = runs.bathy.h(2:end-1,2:end-1);
        %        csr = runs.rgrid.Cs_r(2:end-1);
        %csw = runs.rgrid.Cs_w(2:end-1);

        % for depth-averaging
        hmax = max(abs(zint)); % max. depth of integration
        hmat = h .* (h <= hmax) + hmax .* (h > hmax);
        % add in sponge mask
        hmat = hmat .* fillnan(~runs.sponge(2:end-1, 2:end-1), 0);

        % for bottom friction I need to mask out the area that
        % doesn't touch the bottom
        hbfric = fillnan(hmat .* (hmat == runs.bathy.h(2:end-1,2:end-1)), ...
                         0);

        xavg = avg1(avg1(xvor,1),2)/1000; yavg = avg1(avg1(yvor,1),2)/1000;

        % AREA AVERAGING - for bottom friction terms
        dA = 1./runs.rgrid.pm(2:end-1,2:end-1)' .* 1./runs.rgrid.pn(2:end-1, ...
                                                          2:end-1)';
        dA = dA .* ~runs.sponge(2:end-1, 2:end-1);
        area = sum(dA(:));

        % VOLUME AVERAGING
        % 2D array - water column volume for each (x,y) - masked
        %dVxy = 1./runs.rgrid.pm(2:end-1,2:end-1)' .* 1./runs.rgrid.pn(2:end-1, 2:end-1)' ...
        %.* hmat;

        % 3D array - cell volumes for each (x,y,z)
        % nansum(dV(:))  ~= nansum(dVxy(:)) since,
        % i'm integrating to a level just
        % above the bottom.
        zmat = repmat(permute(zrnew, [3 2 1]), [size(hmat,1) ...
                            size(hmat,2)]);
        zmat(bsxfun(@lt, zmat, -1 * hmat)) = NaN;
        dV = bsxfun(@times, ...
                    bsxfun(@times, dA, diff(zmat, 1, 3)), ...
                    ~isnan(hmat));
        vol = nansum(dV(:));
        %disp(['error in volumes = ' num2str((vol - nansum(dVxy(:)))./vol ...
        %                                    * 100) ' percent']);

        % time range and file reading parameters
        slab = 12;
        stride = 1;

        timehis = dc_roms_read_data(runs.dir, 'ocean_time', [], {}, ...
                                    [], runs.rgrid, 'his');
        trange = tind:stride:length(timehis);
        disp(['starting from t instant = ' num2str(trange(1))]);

        % save vorticity budget for whole domain
        runs.vorbudget.hadv = nan([length(trange) 1]);
        runs.vorbudget.vadv = runs.vorbudget.hadv;
        runs.vorbudget.tilt = runs.vorbudget.hadv;
        runs.vorbudget.str  = runs.vorbudget.hadv;
        runs.vorbudget.beta = runs.vorbudget.hadv;
        runs.vorbudget.bfric = runs.vorbudget.hadv;
        %runs.vorbudget.sol = runs.vorbudget.hadv;
        %runs.vorbudget.budget = runs.vorbudget.hadv;

        % vorticity budget for shelf water
        runs.vorbudget.shelf.hadv = runs.vorbudget.hadv;
        runs.vorbudget.shelf.vadv = runs.vorbudget.hadv;
        runs.vorbudget.shelf.str = runs.vorbudget.hadv;
        runs.vorbudget.shelf.tilt = runs.vorbudget.hadv;
        runs.vorbudget.shelf.beta = runs.vorbudget.hadv;
        runs.vorbudget.shelf.bfric = runs.vorbudget.hadv;

        runs.vorbudget.comment = ['hadv + vadv + beta = str + tilt  ' ...
                            '+ bfric'];
        %runs.vorbudget.conthis = runs.vorbudget.hadv;
        %%

        for kk=1:slab:length(trange)
            tt = trange(kk);
            disp(['kk = ' num2str(kk/slab) '/' num2str(length(trange)/slab) ...
                  ' | tt = ' num2str(tt/2) ' days | plotflag = ' ...
                  num2str(plotflag) ' | run = ' runs.name]);
            %zeta = runs.zeta(2:end-1,2:end-1,tt);
            % read data
            %fname = [runs.dir '/ocean_his.nc.new2'];
            %fname = runs.out_file;
                        %w  = double(ncread(fname,'w',[1 1 1 tt],[Inf Inf Inf 1]));
            %zeta = double(ncread(fname,'zeta',[1 1 tt],[Inf Inf 1]));

            if stride ~= 1
                error('stride does not work');
            end

            % read in history file data
            tindices = [tt tt+stride*slab-1]
            if tt+stride*slab-1 > trange(end)
                tindices(end) = trange(end)
            end
            uh = dc_roms_read_data(runs.dir,'u',tindices,{},[],runs.rgrid, ...
                                  'his', 'single');
            vh = dc_roms_read_data(runs.dir,'v',tindices,{},[],runs.rgrid, ...
                                  'his', 'single');
            wh = dc_roms_read_data(runs.dir,'w',tindices,{},[],runs.rgrid, ...
                                   'his', 'single');
            csdye = dc_roms_read_data(runs.dir, runs.csdname, tindices, ...
                                      {}, [], runs.rgrid, 'his', 'single');
            %zeta = dc_roms_read_data(runs.dir, 'zeta', tt, {}, [], ...
            %                         runs.rgrid, 'his', 'single');

            %rhoh = dc_roms_read_data(runs.dir,'rho',tt,{},[],runs.rgrid, ...
            %                         'his');

            %ubar = dc_roms_read_data(runs.dir, 'ubar', tt, {}, [], ...
            %                         runs.rgrid, 'his');
            %vbar = dc_roms_read_data(runs.dir, 'vbar', tt, {}, [], ...
            %                         runs.rgrid, 'his');

            % interpolate to znew depths
            disp('interpolating variables');
            u = single(interpolate(uh, gridu.zmat, zrnew));
            v = single(interpolate(vh, gridv.zmat, zrnew));
            w = single(interpolate(wh, gridw.zmat, zwnew));
            csd = single(interpolate(csdye, gridr.zmat, zrnew));
            % rho = interpolate(rhoh, gridr.zmat, zrnew);

            ux = bsxfun(@rdivide, diff(u,1,1), diff(gridu.xmat,1,1));
            uy = bsxfun(@rdivide, diff(u,1,2), diff(gridu.ymat,1,2));
            uz = bsxfun(@rdivide, diff(u,1,3), diff(gridu.znew,1,3));

            vx = bsxfun(@rdivide, diff(v,1,1), diff(gridv.xmat,1,1));
            vy = bsxfun(@rdivide, diff(v,1,2), diff(gridv.ymat,1,2));
            vz = bsxfun(@rdivide, diff(v,1,3), diff(gridv.znew,1,3));

            wx = bsxfun(@rdivide, diff(w,1,1), diff(gridw.xmat,1,1));
            wy = bsxfun(@rdivide, diff(w,1,2), diff(gridw.ymat,1,2));
            wz = bsxfun(@rdivide, diff(w,1,3), diff(gridw.znew,1,3));

            %rx = diff(rho,1,1)./diff(gridr.xmat,1,1);
            %ry = diff(rho,1,2)./diff(gridr.ymat,1,2);
            %rz = diff(rho,1,3)./diff(gridr.znew,1,3);

            %cont = ux(:, 2:end-1, :) + vy(2:end-1, :, :) + ...
            %       wz(2:end-1, 2:end-1, :);

            % check cont
            %ix = 150; iy = 164;
            %ix = 240; iy = 164
            %figure; hold all;
            %hold all;
            %plot(squeeze(ux(ix, iy+1,:)) + squeeze(vy(ix+1, iy,:)), zrnew)
            %plot(-1*squeeze(wz(ix+1, iy+1,:)), zrnew)
            %plot(squeeze(cont(ix, iy, :)), zrnew);
            %legend('ux + vy', 'wz', 'ux + vy + wz');

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
            rvx = bsxfun(@rdivide, diff(rv,1,1), diff(gridrv.xmat,1,1));
            rvy = bsxfun(@rdivide, diff(rv,1,2), diff(gridrv.ymat,1,2));
            rvz = bsxfun(@rdivide, diff(rv,1,3), diff(gridrv.znew,1,3));

            rvavg = avg1(avg1(avg1(rv, 1), 2), 3);

            if debug
                u1h = double(ncread(runs.dir,'u',[1 1 1 tt+1],[Inf Inf Inf 1]));
                v1h = double(ncread(runs.dir,'v',[1 1 1 tt+1],[Inf Inf ...
                                    Inf 1]));
                t1 = double(ncread(runs.dir, 'ocean_time'));

                u1 = interpolate(u1h, gridu.zmat, zrnew);
                v1 = interpolate(v1h, gridv.zmat, zrnew);

                u1y = diff(u1,1,2)./diff(gridu.ymat,1,2);
                v1x = diff(v1,1,1)./diff(gridv.xmat,1,1);

                rv1 = v1x-u1y;

                % calculate term and average to agree with 'budget' size
                drvdt = avg1(avg1(avg1( ...
                    (rv1-rv)./(t1(tt+1)-t1(tt)), 1), 2), 3);

                DRVDT = trapz(zint, repnan(drvdt, 0), 3)./hmat;
            end

            str = avg1(avg1(avg1(bsxfun(@plus, rv, ...
                                             avg1(avg1(runs.rgrid.f',1),2)),1) ...
                                 ,2) .* -1 .* (ux(:,2:end-1,:,:) + vy(2:end-1,:,:,:)),3);
            %wz(2:end-1,2:end-1,:), 3);%

            tilt = -1 * avg1(avg1( avg1(wx(:,:,2:end-1,:),2) .* avg1(vz,1) + ...
                    avg1(wy(:,:,2:end-1,:),1) .* avg1(uz,2) ,1),2);
            beta = avg1(avg1(runs.params.phys.beta * v(2:end-1,:,:,:),2),3);
            hadv = avg1( avg1(u(:,2:end-1,:,:),1) .* avg1(rvx,2) + ...
                    avg1(v(2:end-1,:,:,:),2) .* avg1(rvy,1),3);
            vadv = avg1(avg1( avg1(avg1(w(:,:,2:end-1,:),1),2) .* rvz ...
                              ,1),2);

            budget = str + tilt - hadv - vadv - beta;

            % shelf water budget
            % shelf water mask defined with csdye + I remove sponge
            % region based on filtering already done in hmat
            shelfmask = bsxfun(@times, (avg1(csd(2:end-1, 2:end-1, :, :),3) < ...
                         runs.bathy.xsb), ~isnan(hmat));
            %shelfmaskrv = bsxfun(@times, avg1(avg1(csd,1),2) < ...
            %                            runs.bathy.xsb, ~isnan(hmat));
            %  sol = -runs.params.phys.g/runs.params.phys.rho0 .* ...
            %          ( avg1(rx,2) .* avg1(zy,1) - avg1(ry,1) .* avg1(zx,2));

            % need bottom vorticity for bfric calculation
            rvbot = nan(size(squeeze(rvavg(:,:,1,:))));
            shelfmaskbot = rvbot;

            if runs.params.misc.rdrg ~= 0
                tic;
                disp('calculating bottom vorticity');
                if kk == 1
                    % valid cells never change, so save mask
                    % (botmask) that when multiplied with field
                    % gives me the bottom values.
                    botmask  = nan(size(rvavg(:,:,:,1)));
                    for kkk = 1:size(rvbot, 3)
                        for iii = 1:size(rvbot,1)
                            for jjj = 1:size(rvbot,2)
                                % locate first 0  since z=1 is bottom
                                zind = find(isnan(squeeze(rvavg(iii,jjj,:,kkk))) ...
                                            == 0, 1, 'first');
                                if ~isempty(zind)
                                    botmask(iii,jjj,zind) = 1;
                                end
                            end
                        end
                    end
                end
                rvbot = squeeze(nansum(bsxfun(@times, rvavg, botmask),3));
                shelfmaskbot = squeeze(nansum(bsxfun(@times, shelfmask, ...
                                                     botmask),3));
                ubot = squeeze(nansum(bsxfun(@times, avg1(avg1(u(:, ...
                                                                 2: ...
                                                                 end-1,:,:),3),1), botmask), 3)) .* shelfmaskbot;
                ;
                vbot = squeeze(nansum(bsxfun(@times, avg1(avg1(v(2: ...
                                                                 end-1,:,:,:),3),2), botmask), 3)) .* shelfmaskbot;;

                toc;
            end

            % depth INTEGRATED QUANTITIES
            RV   = avg1(avg1(trapz(zrnew, repnan(rv,0), 3),1), 2);
            %RVSHELF = trapz(zrnew, repnan(rvavg,0) .* ...
            %                shelfmask, 3);

            % depth - AVERAGED quantities for plotting
            STR  = squeeze(bsxfun(@rdivide, trapz(zint, repnan(str,0),  3), hmat));
            TILT = squeeze(bsxfun(@rdivide, trapz(zint, repnan(tilt,0), 3), hmat));
            BETA = squeeze(bsxfun(@rdivide, trapz(zint, repnan(beta,0), 3), hmat));
            HADV = squeeze(bsxfun(@rdivide, trapz(zint, repnan(hadv,0), 3), hmat));
            VADV = squeeze(bsxfun(@rdivide, trapz(zint, repnan(vadv,0), 3), hmat));
            ADV = HADV + VADV;

            % FRICTION only when integrating to bottom surface
            BFRIC = bsxfun(@times, bsxfun(@rdivide, -runs.params.misc.rdrg .* rvbot, ...
            hmat), hmat == h);
            BFRICSHELF = BFRIC .* shelfmaskbot;
            bfric = bsxfun(@times, -runs.params.misc.rdrg .* rvbot, ...
                           hmat == h);
            bfricshelf = bfric .* shelfmaskbot;

            % BUDGET = TEND = d(RV)/dt
            %BUD = STR + BFRIC + TILT - BETA - ADV;

            % ubar, vbar calculated for depth averaged interval
            % only
            % ubar = bsxfun(@rdivide, trapz(zrnew, repnan(avg1(u(:,2:end-1,:,:),1),0), 3) ...
            %            ,hmat);
            % vbar = bsxfun(@rdivide, trapz(zrnew, repnan(avg1(v(2:end-1,:,:,:),2),0), ...
            %              3), hmat);
            if debug
                BUD = BUD - DRVDT;
                imagesc(BUD');
            end
            %BUD = trapz(zint, repnan( str+tilt - beta - hadv
            %-vadv,0), 3);

            % volume of shelfwater - shelfmask has sponge taken out
            shelfvol = bsxfun(@times, shelfmask, dV);
            shelfvol = squeeze(nansum(nansum(nansum(shelfvol, 1), ...
                                             2), 3));

            % area of shelfwater in contact with bottom
            shelfarea = bsxfun(@times, shelfmaskbot, dA);
            shelfarea = squeeze(nansum(nansum(shelfarea, 1), 2));

            % reshape for volume averaging
            sz4d = size(rvavg);
            if length(sz4d) == 3
                sz4d(4) = 1;
            end
            sz2d = [sz4d(1)*sz4d(2)*sz4d(3) sz4d(4)];

            % calculate vorticity eqn terms - with shelfmask -
            % volume averaged
            indices = [tindices(1):tindices(end)] - trange(1) + 1;
            runs.vorbudget.shelf.vol(indices) = shelfvol;
            %            runs.vorbudget.shelf.area(indices) = shelfarea;
            runs.vorbudget.shelf.rv(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, rvavg .* shelfmask, dV),1), 2), 3)) ...
                ./ shelfvol;
            runs.vorbudget.shelf.str(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, str .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
            runs.vorbudget.shelf.tilt(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, tilt .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
            runs.vorbudget.shelf.hadv(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, hadv .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
            runs.vorbudget.shelf.vadv(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, vadv .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
            runs.vorbudget.shelf.beta(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, beta .* shelfmask, dV),1), 2), 3)) ./ shelfvol;
            runs.vorbudget.shelf.bfric(indices) = squeeze(nansum(nansum( ...
                bsxfun(@times, bfricshelf, dA), 1), 2)) ./ ...
                shelfvol;

            bfricold = squeeze(nansum(nansum( ...
                bsxfun(@times, BFRICSHELF, dA), 1), 2)) ./ shelfarea;


            % save volume averaged quantities for whole domain
            runs.vorbudget.rv(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, rvavg, dV),1), 2), 3)) ...
                ./ vol;
            runs.vorbudget.str(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, str, dV),1), 2), 3)) ./ vol;
            runs.vorbudget.tilt(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, tilt, dV),1), 2), 3)) ./ vol;
            runs.vorbudget.hadv(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, hadv, dV),1), 2), 3)) ./ vol;
            runs.vorbudget.vadv(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, vadv, dV),1), 2), 3)) ./ vol;
            runs.vorbudget.beta(indices) = squeeze(nansum(nansum(nansum( ...
                bsxfun(@times, beta, dV),1), 2), 3)) ./ vol;
            runs.vorbudget.bfric(indices) = squeeze(nansum(nansum( ...
                bsxfun(@times, bfric, dA), 1), 2)) ./ vol;

            if plotflag
                limc = [-1 1] * nanmax(abs(ADV(:)));
                limy = [0 150];
                limx = [xvor(find(~runs.sponge(:,1) == 1, 1, 'first'),1) ...
                        xvor(find(~runs.sponge(:,1) == 1, 1, 'last'),1)]/1000;
                titlestr = 'Depth integrated rvor';
                % plot
                if kk == 1
                    figure; maximize();
                    ax(1) = subplot(2,4,[1:2]);
                    hvor = pcolor(xavg, yavg, RV); hold on; shading flat;
                    axis image;
                    ht = runs.set_title('Depth int rvor', ceil(tt/2));
                    he(1) = runs.plot_eddy_contour('contour',ceil(tt/2));
                    hbathy = runs.plot_bathy('contour','k');
                    shading flat
                    caxis([-1 1] * nanmax(abs(RV(:))));
                    colorbar;
                    ylim(limy); xlim(limx);

                    ax(2) = subplot(2,4,3);

                    if runs.params.misc.rdrg == 0
                        hbet = pcolor(xavg,yavg,-BETA);
                        title('- \beta V');
                    else
                        hbet = pcolor(xavg, yavg, BFRIC);
                        title('Bottom Friction');
                    end
                    colorbar; shading flat;
                    he(2) = runs.plot_eddy_contour('contour', ceil(tt/2));
                    hbathy = runs.plot_bathy('contour','k');
                    caxis(limc); %caxis([-1 1] * nanmax(abs(BETA(:))));

                    ylim(limy); xlim(limx);

                    ax(3) = subplot(2,4,4); cla
                    xran = 1:6:size(xavg,1); yran = 1:4:size(yavg,2);
                    hquiv = quiver(xavg(xran,yran),yavg(xran,yran), ...
                                   ubar(xran,yran), vbar(xran,yran),1.5);
                    title('(ubar,vbar)');
                    he(3) = runs.plot_eddy_contour('contour', ceil(tt/2));
                    hbathy = runs.plot_bathy('contour','k');
                    ylim(limy); xlim(limx);

                    %                 ax(4) = subplot(2,4,5);
                    %                 htend = pcolor(xavg,yavg,TEND); colorbar; shading flat;
                    %                 he(4) = runs.plot_eddy_contour('contour', ceil(tt/2));
                    %                 hbathy = runs.plot_bathy('contour','k');
                    %                 caxis([-1 1] * max(abs(TEND(:))));
                    %                 title('d\xi/dt');

                    ax(5) = subplot(2,4,7);
                    hgadv = pcolor(xavg,yavg,-ADV); colorbar; shading flat;
                    he(5) = runs.plot_eddy_contour('contour', ceil(tt/2));
                    hbathy = runs.plot_bathy('contour','k');
                    caxis(limc); %caxis([-1 1] * max(abs(ADV(:))));
                    title('-Advection');

                    ax(6) = subplot(2,4,8);
                    htilt = pcolor(xavg,yavg,TILT); colorbar; shading flat;
                    he(6) = runs.plot_eddy_contour('contour', ceil(tt/2));
                    hbathy = runs.plot_bathy('contour','k');
                    caxis(limc/10); %caxis([-1 1] * max(abs(TILT(:))));
                    title('Tilting');

                    ax(7) = subplot(2,4,[5 6]);
                    hstr = pcolor(xavg,yavg,STR); colorbar; hold on; shading flat;
                    %hquiv = quiverclr(xavg(xran,yran),yavg(xran,yran), ...
                    %    ubar(xran,yran),vbar(xran,yran),0.3,STR(xran,yran), ...
                    %    [-1 1]*1e-11);
                    %set(gca,'color',[0 0 0]);
                    he(7) = runs.plot_eddy_contour('contour',ceil(tt/2));
                    hbathy = runs.plot_bathy('contour','k');
                    caxis(limc); %caxis([-1 1] * max(abs(STR(:))));
                    title('Stretching = (f+\xi)w_z')
                    spaceplots(0.06*ones([1 4]),0.05*ones([1 2]))
                    linkaxes(ax,'xy');
                    runs.video_update();
                    pause();
                else
                    set(hvor ,'cdata',RV);
                    set(hgadv ,'cdata',-ADV);
                    if runs.params.misc.rdrg == 0
                        set(hbet ,'cdata',-BETA);
                    else
                        set(hbet, 'cdata', BFRIC);
                    end
                    set(hstr ,'cdata',STR);
                    set(htilt,'cdata',TILT);
                    %set(htend,'cdata',TEND);
                    try
                        set(hquiv,'udata',ubar(xran,yran),'vdata',vbar(xran,yran));
                    catch ME
                    end

                    runs.update_eddy_contour(he,ceil(tt/2));
                    runs.update_title(ht,titlestr,ceil(tt/2));
                    runs.video_update();
                    pause(0.01);
                end
            end
        end

        if plotflag
            runs.video_write();
        end
        runs.vorbudget.time = timehis(trange);

        figure;
        plot(runs.vorbudget.time,-runs.vorbudget.hadv,'r'); hold on
        plot(runs.vorbudget.time,-runs.vorbudget.vadv,'g');
        plot(runs.vorbudget.time,runs.vorbudget.tilt,'b');
        plot(runs.vorbudget.time,runs.vorbudget.str,'c');
        %plot(runs.vorbudget.time,runs.vorbudget.sol,'m');
        plot(runs.vorbudget.time,-runs.vorbudget.beta,'y');
        %       plot(runs.vorbudget.time,runs.vorbudget.budget,'k');
        title('signs so that all terms are on RHS and tendency is LHS');
        legend('hadv','vadv','tilt','str','beta');

        vorbudget = runs.vorbudget;
        vorbudget.hash = githash;
        save([runs.dir '/vorbudget.mat'],'vorbudget');

        toc(vorbudgetstart);
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
        vorbudget.hash = githash;
        save([runs.dir '/vorbudget.mat'],'vorbudget');
    end

    function [] = animate_center(runs)
        runs.video_init('center');
        eddy = runs.eddy;
        xvec = runs.rgrid.xr(:,1);
        yvec = runs.rgrid.yr(1,:)';

        varname = 'rho';

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
            temper = dc_roms_read_data(runs.dir,varname,[t0 st Inf], ...
                            {'y' iymin iymax},stride,runs.rgrid, ...
                                       'his', 'single');
            strat = dc_roms_read_data(runs.dir,varname,[1 1], ...
                            {'y' Inf Inf},stride, runs.rgrid, 'his', ...
                                      'single');

            temper = bsxfun(@minus,temper,permute(strat,[1 3 2]));
            %temper = roms_read_data(runs.dir,varname,[1 iymin 1 t0], ...
            %                  ceil([Inf iymax-iymin+1 Inf Inf]./stride), stride);
            %              toc;
            %strat  = roms_read_data(runs.dir,varname,[Inf 1 1 1], ...
            %                  ceil([1 1 Inf 1]./stride),stride);
            %              toc

        else
            stride = [1 sxy sz st];
            temper = dc_roms_read_data(runs.dir,varname,[t0 st Inf], ...
                            {'x' ixmin ixmax},stride, runs.rgrid, ...
                                       'his', 'single');
            strat = dc_roms_read_data(runs.dir,varname,[1 1], ...
                            {'y' Inf Inf},stride, runs.rgrid, 'his', ...
                                      'single');
            temper = bsxfun(@minus,temper,permute(strat,[3 1 2]));
                        %temper = roms_read_data(runs.dir,varname,[ixmin 1  1 t0], ...
            %                ceil([ixmax-ixmin+1 Inf Inf Inf]./stride),stride);
            %            toc;
            %strat  = roms_read_data(runs.dir,varname,[1 1 1 1], ...
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
        runs.read_zeta;
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

        range = ['runs.spng.sx1:runs.spng.sx2,' ...
                 'runs.spng.sy1:runs.spng.sy2'];


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
        zend = runs.zeta(:,:,end);
        if ~isnan(zend)
            caxis([min(zend(:)) max(zend(:))]);
        end
        center_colorbar;
    end
    function update_zeta(runs,handle,tt)
        try
            %handle.CData = runs.zeta(:,:,tt);
            set(handle,'CData',double(runs.zeta(:,:,tt)));
        catch ME
            %handle.ZData = runs.zeta(:,:,tt);
            set(handle,'ZData',double(runs.zeta(:,:,tt)));
        end
    end

    function [hplot] = plot_surf(runs,varname,plottype,tt)
        if ~exist('tt','var'), tt = 1; end

        if ~strcmpi(varname, 'pbot')
            range = ['runs.spng.sx1:runs.spng.sx2,' ...
                     'runs.spng.sy1:runs.spng.sy2'];
            axrange  = range;
            eval(['zend = runs.' varname '(' range ',end);']);
        else
            range = ':,:';
            axrange = ['runs.spng.sx1+2:runs.spng.sx2-2,' ...
                       'runs.spng.sy1+1:runs.spng.sy2-2'];
            zend = runs.pbot(:,:,end);
        end

        if ~isnan(zend)
            crange = double([min(zend(:)) max(zend(:))]);
            caxis(crange);
        end

        if strcmpi(plottype,'pcolor')
            eval(['hplot = pcolor(runs.rgrid.xr(' axrange ')/1000,' ...
                  'runs.rgrid.yr(' axrange ')/1000,' ...
                  'double(runs.' varname '(' range ',tt)));']);
            if runs.makeVideo
                shading interp;
            else
                shading flat
            end
        else
            if strcmpi(plottype,'contourf') || strcmpi(plottype,'contour')
                eval(['[cc,hplot] = ' plottype ...
                      '(runs.rgrid.xr(' range ')/1000,' ...
                      'runs.rgrid.yr(' range ')/1000,'...
                      'double(runs.' varname '(' range ',tt)));']);
                hplot.LevelList = linspace(crange(1),crange(2), 10);
                shading flat
            end
        end

        if strcmpi(varname, 'eddye')
            center_colorbar;
        end
    end
    function update_surf(runs,varname,handle,tt)

        if ~strcmpi(varname, 'pbot')
            range = ['runs.spng.sx1:runs.spng.sx2,' ...
                     'runs.spng.sy1:runs.spng.sy2'];
            axrange  = range;
            eval(['zend = runs.' varname '(' range ',end);']);
        else
            range = ':,:';
            axrange = ['runs.spng.sx1+2:runs.spng.sx2-2,' ...
                       'runs.spng.sy1+1:runs.spng.sy2-2'];
            zend = runs.pbot(:,:,end);
        end

        try
            eval(['set(handle,''CData'',double(runs.' varname '(' ...
                  range ',tt)))']);
        catch ME
            levels = handle.LevelList;
            eval(['set(handle,''ZData'',double(runs.' varname '(' ...
                  range ',tt)))']);
            handle.LevelList = levels;
        end
    end


    function [hplot] = plot_eddy_contour(runs,plottype,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2],1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        hold on;
        [~,hplot] = contour(runs.eddy.xr(ix,iy)/1000,runs.eddy.yr(ix,iy)/1000, ...
                    runs.eddy.vormask(ix,iy,tt),'Color','k','LineWidth',1);
    end

    function update_eddy_contour(runs,handle,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2],1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        for ii=1:length(handle)
            try
                set(handle(ii),'ZData', ...
                               double(runs.eddy.vormask(ix,iy,tt)));
            catch ME
                set(handle(ii),'CData', ...
                               double(runs.eddy.vormask(ix,iy,tt)));
            end
        end
    end


    function [hplot] = plot_rho_contour(runs,plottype,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2],1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        mask = ((runs.rhosurf(2:end-1,2:end-1,tt) - runs.rbacksurf) < ...
               runs.eddy.drhothresh(1)); % .* runs.eddy.vormask(:,:,tt);
        hold on;

        [~,hplot] = contour(runs.eddy.xr(ix,iy)/1000, ...
                            runs.eddy.yr(ix,iy)/1000, ...
                            mask(ix,iy),'Color',[44 162 95]/256, ...
                            'LineWidth',1);
    end
    function update_rho_contour(runs,handle,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2],1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        mask = ((runs.rhosurf(2:end-1,2:end-1,tt) - runs.rbacksurf) < ...
                runs.eddy.drhothresh(1)) .* runs.eddy.vormask(:,:,tt);

        for ii=1:length(handle)
            try
                set(handle(ii),'ZData', ...
                               double(mask(ix,iy)));
            catch ME
                disp(ME.message);
            end
        end
    end


    function [hplot] = plot_eddy_sshcontour(runs,plottype,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2],1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        hold on;
        [~,hplot] = contour(runs.eddy.xr(ix,iy)/1000, ...
                            runs.eddy.yr(ix,iy)/1000, ...
                            runs.eddy.mask(ix,iy,tt), ...
                            'Color','k','LineWidth',1);
    end

    function update_eddy_sshcontour(runs,handle,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2],1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        for ii=1:length(handle)
            try
                set(handle(ii),'ZData', run.eddy.mask(ix,iy,tt));
            catch ME
            end
        end
    end

    function [ht] = set_title(runs,titlestr,tt)
        ht = title([titlestr ' | ' runs.name ' | ' ...
                    num2str(runs.time(tt)/86400)  ' days | t_{nd} = ' ...
                   num2str(runs.time(tt)/runs.eddy.turnover)]);
    end
    function update_title(runs,ht,titlestr,tt)
        set(ht,'String',[titlestr ' | ' runs.name ' | ' ...
                         num2str(runs.time(tt)/86400)  ' days | t_{nd} = ' ...
                          num2str(runs.time(tt)/runs.tscale)]);
    end

    function [hplot] = plot_bathy(runs,plottype,color)
        ix = runs.spng.sx1:runs.spng.sx2;
        iy = runs.spng.sy1:runs.spng.sy2;

        if ~exist('color','var'), color = [1 1 1]*0.85; end
        if strcmpi(plottype,'contour')
            [cc,hplot{1}] = contour(runs.rgrid.xr(ix,iy)/1000,...
                                 runs.rgrid.yr(ix,iy)/1000, ...
                                 runs.rgrid.h(iy,ix)',[200 500 1000 1500 ...
                                2000], 'Color', color);
            clabel(cc, hplot{1}, 'LabelSpacing', 108*2.75, 'Color', color);
            hax = gca;
            sbslcolor = color;
            if runs.bathy.axis == 'y'
                hplot{2} = liney(runs.bathy.xsb/1000,[], sbslcolor);
                hplot{3} = liney(runs.bathy.xsl/1000,[], sbslcolor);
            else
                hplot{2} = linex(runs.bathy.xsb/1000,[], sbslcolor);
                hplot{3} = linex(runs.bathy.xsl/1000,[], sbslcolor);
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