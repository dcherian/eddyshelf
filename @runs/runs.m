classdef runs < handle
properties
    % dir & file names
    name; dir; out_file; ltrans_file; flt_file; givenFile; tracpy_file; ...
        fpos_file;
    % data
    zeta; temp; usurf; vsurf; vorsurf; csdsurf; ubot; vbot; eddsurf; ...
        rhosurf; edcsdyesurf; pvsurf; zdyesurf;
    rbacksurf; % background density at surface
    sgntamp; % sign(runs.eddy.tamp) = -1 if cyclone; 1 otherwise
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
    % volume budget
    volume;
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
    % along-shore supply jet properties
    supply;
    % SSH scale at shelfbreak
    sbssh;
    % shelf baroclinicity
    shelfbc;
    % ubar scale
    ubarscale;
    % along-shore slope jet properties
    jet;
    % threshold values
    eddy_thresh = 0.7;
    % initial params
    params
    % fluxes - cross-shore & along-shore; - energy fluxes
    csflux; asflux;
    % fluxes within a radius of eddy center
    radius;
    % stats of eddy water on shelf
    onshelf;
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
    function [runs] = runs(dir, reset,  do_all, reduced)

        ticstart = tic;

        disp('=======================')
        disp(['runs(' dir ')']);
        disp('=======================')

        read_zeta = 0;

        if ~exist('reset','var'), reset = 0; end
        if ~exist('do_all', 'var'), do_all = 0; end
        if ~exist('reduced', 'var'), reduced = 0; end

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
        runs.rgrid = dc_roms_get_grid(runs.out_file,runs.out_file, zeta0',1);
        runs.rgrid.xr = runs.rgrid.x_rho';
        runs.rgrid.yr = runs.rgrid.y_rho';
        runs.rgrid.dx = mean(1./runs.rgrid.pm(:));
        runs.rgrid.dy = mean(1./runs.rgrid.pn(:));

        % remove needless rgrid matrices
        runs.rgrid.angle = [];

        runs.rbacksurf = ncread(runs.out_file, 'rho', [1 1 runs.rgrid.N 1], [1 1 1 1]);

        % read zeta
        if ~runs.givenFile
            if read_zeta
                runs.zeta = dc_roms_read_data(dir,'zeta',[],{},[],runs.rgrid, ...
                                              'his', 'single');
            end
            runs.time = dc_roms_read_data(dir,'ocean_time',[],{}, ...
                                          [],runs.rgrid, 'his', 'single');
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

        runs.sgntamp = sign(runs.params.eddy.tamp);
        if isnan(runs.params.bg.ubt)
            runs.params.bg.ubt = 0;
        end
        if isnan(runs.params.bg.vbt)
            runs.params.bg.vbt = 0;
        end

        % fill bathy
        [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = ...
                        find_shelfbreak(runs.out_file);
        [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = ...
                        find_shelfbreak(runs.out_file,'slope');
        runs.bathy.h = runs.rgrid.h';
        runs.bathy.betash = runs.params.phys.f0/runs.bathy.hsb * runs.bathy.sl_shelf;
        runs.bathy.betasl = runs.params.phys.f0/runs.bathy.hsb * runs.bathy.sl_slope;

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
        if isempty(runs.spng.sx2)
            runs.spng.sx2 = size(runs.sponge, 1) - 2;
        end
        if isempty(runs.spng.sy2)
            runs.spng.sy2 = size(runs.sponge, 2) - 2;
        end

        % rossby radii
        runs.rrdeep = sqrt(runs.params.phys.N2)*max(runs.bathy.h(:)) ...
                    /abs(mean(runs.rgrid.f(:)))/pi;
        runs.rrshelf = sqrt(runs.params.phys.N2)*max(runs.bathy.hsb) ...
                    /abs(mean(runs.rgrid.f(:)))/pi;

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
        if ~exist([dir '/eddytrack.mat'],'file') || reset == 1
            try
                runs.eddy = track_eddy(dir);
                runs.noeddy = 0;
            catch ME
                disp(ME.message);
                disp('Couldn''t run track_eddy.m');
                runs.noeddy = 1;
            end
        else
            edd = load([dir '/eddytrack.mat'],'eddy');
            runs.eddy = edd.eddy;
            runs.noeddy = 0;
        end


        % extra processing of eddy track
        if ~runs.noeddy
            runs.params.nondim.eddy.Bu = (runs.params.phys.f0 * ...
                                          runs.params.eddy.dia/2 / ...
                                          runs.params.eddy.depth).^2 / ...
                runs.params.phys.N2;

            % scale time by eddy translation
            if runs.bathy.axis == 'y'
                runs.eddy.tscaleind = find_approx(runs.eddy.my, runs.bathy.xsl, 1);
            else
                runs.eddy.tscaleind = find_approx(runs.eddy.mx, runs.bathy.xsl, 1);
            end
            runs.eddy.tscale = runs.eddy.t(runs.eddy.tscaleind) .* 86400;

            % find time when southern edge crosses slopebreak
            runs.eddy.edgtscaleind = find_approx(runs.eddy.vor.se, runs.bathy.xsl, 1);
            runs.eddy.edgtscale = runs.eddy.t(runs.eddy.edgtscaleind) * 86400;

            if ~isfield(runs.eddy, 'tend') | (runs.eddy.tend == 0)
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
            try
                runs.eddy.turnover = (runs.eddy.rhovor.dia(1)/2) ./ runs.eddy.V(1);
            catch ME
                runs.eddy.turnover = (runs.eddy.vor.dia(1)/2) ./ runs.eddy.V(1);
            end

            if ~isfield('imx', runs.eddy)
                runs.eddy.imx = vecfind(runs.rgrid.xr(:,1), runs.eddy.mx);
                runs.eddy.imy = vecfind(runs.rgrid.yr(1,:), runs.eddy.my);
            end

            try
                runs.eddy.fitx.L = addnan(runs.eddy.fitx.L, 1e10);
                runs.eddy.fitx.Lrho = addnan(runs.eddy.fitx.Lrho, 1e10);
                runs.eddy.fity.L = addnan(runs.eddy.fity.L, 1e10);
                runs.eddy.fity.Lrho = addnan(runs.eddy.fity.Lrho, 1e10);
                runs.eddy.fitx.L(runs.eddy.fitx.L < 0) = NaN;
                runs.eddy.fity.L(runs.eddy.fity.L < 0) = NaN;
                runs.eddy.fitx.V0(runs.eddy.fitx.V0 > 1) = NaN;
                runs.eddy.fity.V0(runs.eddy.fity.V0 > 1) = NaN;
            catch
            end

            % non-dimensionalized time
            runs.ndtime = runs.time ./ runs.eddy.turnover;

            % save memory by converting masks to logical
            try
                if ~islogical(runs.eddy.mask)
                    runs.eddy.mask = logical(repnan(runs.eddy.mask, 0));
                    runs.eddy.vor.mask = logical(repnan(runs.eddy.vor.mask, 0));
                    runs.eddy.rhovor.mask = logical(repnan(runs.eddy.rhovor.mask, 0));
                    runs.eddy.rhossh.mask = logical(repnan(runs.eddy.rhossh.mask, 0));
                end
            catch ME
            end

            runs.eddy.xr = single(runs.eddy.xr);
            runs.eddy.yr = single(runs.eddy.yr);

            try
                % drhothresh based on ssh mask if it doesn't exist
                if ~isfield(runs.eddy, 'drhothreshssh')
                    rs = ncread(runs.out_file, 'rho', [1 1 runs.rgrid.N ...
                                        1], [Inf Inf 1 1]);
                    rs = rs(2:end-1,2:end-1) - rs(1,1);
                    runs.eddy.drhothreshssh = ...
                        squeeze(nanmax(nanmax(rs .* ...
                                              fillnan(runs.eddy.mask(:,:,1),0), [], 1), [], 2));
                end
                % drhothresh based on ssh mask if it doesn't exist
                if ~isfield(runs.eddy, 'drhothresh')
                    rs = ncread(runs.out_file, 'rho', [1 1 runs.rgrid.N ...
                                        1], [Inf Inf 1 1]);
                    rs = rs(2:end-1,2:end-1) - rs(1,1);
                    runs.eddy.drhothresh = ...
                        squeeze(nanmax(nanmax(rs .* ...
                                              fillnan(runs.eddy.vor.mask(:,:,1),0), [], 1), [], 2));
                end
            catch ME
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
            % try
            %     runs.eddy.trevind = find(runs.eddy.cvx < 0,1,'first');
            %     runs.eddy.trev = runs.time(runs.eddy.trevind);
            % catch ME
            %     disp('Eddy did not reverse direction');
            %     runs.eddy.trev = nan;
            % end
            % if isempty(runs.eddy.trev), runs.eddy.trev = NaN; end

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

            % remove needless h-matrices
            runs.eddy.h = [];

            if isfield(runs.eddy.vor, 'Ro')
                runs.eddy.Ro = runs.eddy.vor.Ro;
            end

            % calculate βL/f
            try
                dA = 1./runs.rgrid.pm' .* 1./runs.rgrid.pn';
                betal = fillnan(runs.eddy.vor.mask, 0) .* ...
                        abs(bsxfun(@minus, runs.rgrid.f(2:end-1,2:end-1)', ...
                                   permute(runs.eddy.fcen, [3 2 1])))/2;
                betaldA = bsxfun(@times, betal, dA(2:end-1,2:end-1));

                runs.eddy.betahat = squeeze(nansum(nansum(betaldA, 1), 2))' ./ ...
                    runs.eddy.vor.area ./ runs.eddy.fcen';
                runs.eddy.Rh = runs.eddy.Ro ./ runs.eddy.betahat;
            catch ME
            end
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
        if strcmpi(runs.name, 'ew-34')
            try
                runs.roms = floats('roms', runs.dir, runs.rgrid, ...
                                   runs.bathy.xsb, runs.fpos_file);
                runs.filter_cross_sb('roms');
            catch ME
                disp('Reading ROMS floats failed.');
                disp(ME);
            end
        end
        %try
        %    runs.tracpy = floats('tracpy',runs.tracpy_file,runs.rgrid, ...
        %                         runs.bathy.xsb);
        %    runs.filter_cross_sb('tracpy');
        %catch ME
        %    disp('Reading tracpy floats failed.');
        %    disp(ME);
        %end
        %if exist(runs.ltrans_file,'file')
        %    try
        %        runs.ltrans = floats('ltrans',runs.ltrans_file, ...
        %                             runs.rgrid, runs.bathy.xsb);
        %    catch ME
        %        warning('LTRANS data not read in');
        %    end
        %    runs.filter_cross_sb('ltrans');
        %end

        % load streamer data if it exists.
        % if exist([dir '/streamer.mat'], 'file') && reset ~= 1
        %     disp('Loading streamer data');
        %     streamer = load([dir '/streamer.mat'],'streamer');
        %     runs.streamer = streamer.streamer;
        %     clear streamer;
        % end

        % load water mass data
        if exist([dir '/watermass.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading water mass data');
            water = load([dir '/watermass.mat'], 'water');
            runs.water = water.water;
            clear water
            runs.CalcNonShelfWaterCensus;
        end

        % load volume budget data
        if exist([dir '/volume.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading volume budget');
            volume = load([dir '/volume.mat'], 'volume');
            runs.volume = volume.volume;
            clear volume
        end

        % load volume budget data
        if exist([dir '/radius.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading flux diags within radius');
            data = load([dir '/radius.mat'], 'radius');
            runs.radius = data.radius;
            clear data
        end

        % load vorticity budget data
          % load water mass data
        if exist([dir '/vorbudget.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading vorticity budget');
            vorbudget = load([dir '/vorbudget.mat'], 'vorbudget');
            runs.vorbudget = vorbudget.vorbudget;
            clear vorbudget
        end

        % load fluxes if the file exists
        if exist([dir '/fluxes.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading fluxes');
            data = load([dir '/fluxes.mat']);
            if isfield(data, 'csflux'), runs.csflux = data.csflux; end
            if isfield(data, 'asflux'), runs.asflux = data.asflux; end
            clear data

            % find time when flux is increasing up
            try
                trans  = runs.csflux.off.itrans.shelf(:,1);
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

            try
                if ~islogical(runs.csflux.offmask)
                    runs.csflux.offmask = logical(runs.csflux.offmask);
                    runs.csflux.onmask = logical(runs.csflux.onmask);

                    runs.csflux.slopext = single(runs.csflux.slopext);
                    runs.csflux.eddyxt = single(runs.csflux.eddyxt);
                    runs.csflux.off.slopezt = single(runs.csflux.off.slopezt);
                    runs.csflux.off.eddyzt = single(runs.csflux.off.eddyzt);
                    runs.csflux.east.slopezt = single(runs.csflux.east.slopezt);
                    runs.csflux.east.eddyzt = single(runs.csflux.east.eddyzt);
                end
            catch ME
            end
        end

        % load eddy water diagnostics if the file exists
        if exist([dir '/onshelf.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading eddy-on-shelf diagnostics');
            data = load([dir '/onshelf.mat']);
            runs.onshelf = data.onshelf;
            clear data
        end

        % load shelf baroclinicity diagnostics if the file exists
        if exist([dir '/shelfbc.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading shelf baroclinicity diagnostics');
            data = load([dir '/shelfbc.mat']);
            runs.shelfbc = data.shelfbc;
            clear data
        end

        % load jet diagnostics if the file exists
        if exist([dir '/jet.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading jet diagnostics');
            data = load([dir '/jet.mat']);
            runs.jet = data.jet;
            clear data
        end

        % load bottom torque diagnostics if the file exists
        if exist([dir '/bottom.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading bottom torque diagnostics');
            data = load([dir '/bottom.mat']);
            runs.bottom = data.bottom;
            clear data
        end

        % load angular momentum diagnostics if the file exists
        if exist([dir '/angmom.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading angular momentum diagnostics');
            data = load([dir '/angmom.mat']);
            runs.angmom = data.angmom;
            clear data
        end

        % load average streamer section diagnostics if the file exists
        if exist([dir '/avgstreamer.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading streamer diagnostics');
            data = load([dir '/avgstreamer.mat']);
            runs.streamer = data.streamer;
            clear data
        end

        % load average streamer section diagnostics if the file exists
        if exist([dir '/sbssh.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading shelfbreak SSH diagnostics');
            data = load([dir '/sbssh.mat']);
            runs.sbssh = data.sbssh;
            clear data
        end

        if exist([dir '/ubarscale.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading ubar cross isobath scale');
            data = load([dir '/ubarscale.mat']);
            runs.ubarscale = data.ubarscale;
            clear data
        end


        % load average supply jet diagnostics
        if exist([dir '/supplyjet.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading supply diagnostics');
            data = load([dir '/supplyjet.mat']);
            runs.supply = data.supply;
            clear data
        end

        % load trajectory fit
        if exist([dir '/traj.mat'],'file') && reset ~= 1 & ~reduced
            disp('Loading trajectory');
            data = load([dir '/traj.mat']);
            runs.traj = data.traj;
            clear data
        end

        % load eddy energy diagnostics if the file exists
        if exist([dir '/energy.mat'],'file') && reset ~= 1 & ~reduced
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

        % calculate Rhines scale
        try
            runs.bathy.Lbetash = sqrt(runs.eddy.fitx.V0(1) / 2.3 / sqrt(2) / ...
                                      runs.bathy.betash);
        catch ME
            runs.bathy.Lbetash = sqrt(runs.eddy.V(1) / ...
                                      runs.bathy.betash);
        end

        % save some memory
        runs.asflux = [];
        runs.csflux.on = [];
        runs.csflux.off.eddyzt = [];
        runs.csflux.off.slopezt = runs.csflux.off.slopezt(:,:,1,1);

        toc(ticstart);
    end

    function [] = info(runs)
        roms_info(runs.dir);
    end

    function [iloc, xloc] = process_loc(runs, ax, loc)

        if strcmpi(loc, 'sb') & (ax == runs.bathy.axis)
            iloc = runs.bathy.isb;
            xloc = runs.bathy.xsb;
        end

        if ischar(loc)
            xloc = str2double(loc);
            iloc = find_approx(1);
        end
    end

    function [sigma, flux, errflx, V0, L0] = ...
            SlopeFactor(runs, norm_v, norm_L, norm_hsb, norm_Lsh)

        fluxvec = runs.csflux.off.slope(:,1,1);

        if strcmpi(runs.name, 'ew-8040')
            whenstop = 0.90;
        else
            whenstop = 0.85;
        end
        [flux, errflx] = runs.calc_avgflux(fluxvec, 0, 0.05, whenstop);
        [start,stop] = runs.flux_tindices(fluxvec);
        t0 = 1;
        tend = ceil((start+stop)/2);
        [V0, L0, Lz0] = runs.EddyScalesForFlux(t0, tend);

        hsb = runs.bathy.hsb;
        alpha = runs.bathy.sl_shelf;
        Lsh = runs.bathy.L_shelf;
        Lsupp = runs.bathy.Lbetash*1.22;

        if ~exist('norm_v', 'var'), norm_v = V0; end
        if ~exist('norm_L', 'var'), norm_L = L0; end
        if ~exist('norm_hsb', 'var'), norm_hsb = hsb; end
        if ~exist('norm_Lsh', 'var'), norm_Lsh = Lsh; end

        % 1. no Ldef corrections;
        %    assume shelf water hasn't been replaced yet
        %    how much shelf volume can the eddy affect?
        %
        % 2. sometimes I change Lsh when slope becomes steep!
        %    denominator of slfac should have shelf width
        %    that was used in flat-bottom run

        Qsl = (V0 * Lsupp * hsb ...
               * ( (1-alpha*Lsupp/hsb) * (1-exp(-(Lsh/Lsupp))) ...
                   - alpha*Lsh/hsb * exp(-Lsh/Lsupp)));

        sigma = Qsl ...
                / ( norm_v * norm_L * norm_hsb ...
                    * (1 - exp(-(norm_Lsh/norm_L))) );
    end

    function [] = FluxIntegralTimeScale(runs, isobath, factor)
        if ~exist('factor', 'var')
            factor = 2;
        end

        if ~exist('isobath', 'var')
            isobath = 2:4;
        end

        for iso=isobath
            fluxvec = runs.recalculateFlux(factor*runs.bathy.hsb, iso, iso);
            [start,stop] = runs.flux_tindices(fluxvec);
            tvec = double(runs.csflux.time(start:stop));
            tvec = (tvec - tvec(1));
            [tvec, uind] = unique(tvec);
            fluxvec = fluxvec(uind);

            [dof,IT(iso)] = calcdof(fluxvec);
        end

        IT = mean(IT(isobath) * (tvec(2)-tvec(1))/86400)
        turnover = runs.eddy.turnover/86400*2*pi
        IT/turnover
    end

    function [tindex] = process_time(runs, in)

        if ischar(in)
            if strcmpi(in, 'max flux')
                [~,tindex] = runs.calc_maxflux(5);
                return
            end

            if strcmpi(in, 'resistance')
                [~,~,tindex] = runs.locate_resistance;
                return;
            end

            % days
            tindex = find_approx(runs.time, str2double(in)*86400, 1);
        else
            if isinf(in)
                tindex = length(runs.time);
            else
                % actual index
                tindex = in;
            end
        end
    end

    function [varname] = process_varname(runs, in)
        varname = in;

        if strcmpi(in, 'csdye') | strcmpi(in, 'csdsurf')
            varname = runs.csdname;
        end

        if strcmpi(in, 'ssh')
            varname = 'zeta';
        end

        if strcmpi(in, 'eddye')
            varname = runs.eddname;
        end

        if strcmpi(in, 'zdye')
            varname = runs.zdname;
        end

        if strcmpi(in, 'csvel')
            varname = runs.csvelname;
        end

        if strcmpi(in, 'asvel')
            varname = runs.asvelname;
        end
    end

    function [runs] = reload(runs)
        dir = runs.dir;
        keyboard;

        runs.delete;
        runs(dir);
    end

    function [hplt, var, xvec] = plot_profilex(runs, varname, tindex, axname, ix, hax)

        tindex = runs.process_time(tindex);
        varname = runs.process_varname(varname);

        if ~exist('hax')
            figure;
            hax = gca;
            insertAnnotation([runs.name '.plot_profilex']);
        else
            hold on;
        end

        if strcmpi(ix, 'sb')
            ix = runs.bathy.isb;
        end

        if axname == 'x'
            if ~exist('ix', 'var') | isempty(ix) | strcmpi(ix, 'cen')
                ix = runs.eddy.imx(tindex);
            end
        else
            if ~exist('ix', 'var') | isempty(ix) | strcmpi(ix, 'cen')
                ix = runs.eddy.imy(tindex);
            end
        end

        if ~strcmpi(varname, 'zeta') ...
                & ~strcmpi(varname, 'ubar') & ~strcmpi(varname, 'vbar')
            vol = {axname ix ix; 'z' runs.rgrid.N runs.rgrid.N};
        else
            vol = {axname ix ix};
        end

        [var, xax, yax]  = dc_roms_read_data(runs, varname, tindex, vol);

        if axname == 'x'
            % for x = ix, plot against y
            xvec = yax(1,:)/1000;
            mxname = 'my';
        else
            % for y = iy, plot against x
            xvec = xax(:,1)/1000;
            mxname = 'mx';
        end

        if axname ~= runs.bathy.axis
            xvec = xvec - runs.bathy.xsb/1000;
            labx = [upper(axname) ' - ' upper(axname) '_{sb} (km)'];
        else
            eval(['xvec = xvec - runs.eddy.' mxname '(tindex)/1000;']);
            labx = [upper(axname) ' - ' upper(axname) '_{cen} (km)'];
        end

        axes(hax);
        hplt = plot(xvec, var);
        xlabel(labx); ylabel(varname);
    end

    function [hplt] = plotSSHgradient(runs, tindex, hax)

        if ~exist('hax', 'var')
            figure;
            insertAnnotation([runs.name '.plotSSHgradient']);
            hax = gca;
        end

        ix = runs.bathy.isb;

        tind = runs.process_time(tindex);

        if runs.bathy.axis == 'y'
            xvec = runs.rgrid.x_rho(1,:)'; - runs.eddy.mx(tind);
        else
            xvec = runs.rgrid.y_rho(:,1) - runs.eddy.my(tind);
        end

        zeta = dc_roms_read_data(runs, 'zeta', tind, {runs.bathy.axis ix ix});

        dzdx = diff(zeta)./diff(xvec);

        axes(hax);
        hplt = plot(avg1(xvec/1000), dzdx);
    end

    function [] = testAvgSupplyJet(runs)

        [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:, 1, 1));
        tindices = [start stop];

        runs.animate_field('csdye', [], tindices(2), 1);
        linex(runs.supply.x/1000);
    end

    function [] = plot_test1(runs)
        vmagn = hypot(avg1(runs.usurf(:,:,1),2), ...
                      avg1(runs.vsurf(:,:,1),1));
        figure;
        pcolorcen(runs.rgrid.xvor, runs.rgrid.yvor, vmagn);
        clim = caxis;
        hold on
        contour(runs.rgrid.x_v, runs.rgrid.y_v, ...
                  runs.vsurf(:,:,1)');
        contour(runs.rgrid.xvor, runs.rgrid.yvor, ...
                runs.vorsurf(:,:,1), [1 1]*0, 'b');
        caxis(clim);
        linex(runs.eddy.mx(1));
        % this should be zero vorticity contour
        linex([-1 1] * runs.eddy.vor.dia(1)/2 + runs.eddy.mx(1));
        % this should be peak of radial velocity and cartesian v
        linex([-1 1]/sqrt(2) * runs.eddy.vor.dia(1)/2 + ...
              runs.eddy.mx(1));
        center_colorbar;
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

    function [] = checkEddyDiags(runs)

        if ~exist('isobath', 'var'), isobath = 1; end

        [~,~,tind] = runs.locate_resistance;
        for iso=1:length(runs.csflux.x)
            [~,maxloc(iso)] = runs.calc_maxflux(iso);
        end

        figure; hold on
        insertAnnotation([runs.name '.checkEddyDiags']);
        subplot(211); hold on;
        plot(runs.eddy.vor.dia/2000);
        plot(runs.eddy.rhovor.dia/2000);
        plot(smooth(hypot(runs.eddy.fitx.Lrho/1000, ...
                          runs.eddy.fity.Lrho/1000)/sqrt(2), 30));
        plot(smooth(hypot(runs.eddy.fitx.L/1000, ...
                          runs.eddy.fity.L/1000)/sqrt(2), 30));
        legend('L_{vor}','L_{rhovor}','L_{fitrho}','L_{fitv}');
        ylim([0.5 2] * runs.params.eddy.dia/2000);
        linex([tind maxloc]);
        title(runs.name);

        subplot(212); hold on;
        plot(runs.eddy.vor.Vke);
        plot(runs.eddy.rhovor.Vke);
        plot(hypot(runs.eddy.fitx.V0, runs.eddy.fity.V0)/sqrt(2));
        ylim([0.5 2] * runs.eddy.vor.Vke(1));
        legend('V_{vor}','V_{rhovor}','V_{fit}');
        linex([tind maxloc]);
    end

    function [] = plot_velprofiles(runs, vname)

        if ~exist('vname', 'var'), vname = 'u'; end

        [~,~,tind] = runs.locate_resistance();

        it = 1:10:tind;

        if vname == 'u'
            ix = vecfind(runs.rgrid.x_rho(1,:), runs.eddy.mx(it));
            iy = vecfind(runs.rgrid.y_rho(:,1), ...
                         runs.eddy.my(it) - runs.eddy.vor.lmin(it)/3);
        else
            ix = vecfind(runs.rgrid.x_rho(1,:), ...
                         runs.eddy.mx(it) - runs.eddy.vor.lmin(it)/3);
            iy = vecfind(runs.rgrid.y_rho(:,1), runs.eddy.my(it));

        end

        hf = figure; hold all
        lightDarkLines(length(it));
        t0 = 1;
        for tt=1:length(it)
            tt
            u = dc_roms_read_data(runs.dir, vname, it(tt), {'x' ix(tt) ...
                                ix(tt); 'y' iy(tt) iy(tt)}, [], ...
                                  runs.rgrid, 'his');
            figure(hf);
            plot((u./u(end)), ...
                 runs.rgrid.z_u(:,iy(tt),ix(tt)) ./ runs.eddy.Lgauss(1));
        end

        z = [-3:0.01:0];
        prof = (1-erf(abs(z)));
        hplt = plot(prof./prof(end), z, 'b-');
        linex(0); liney(runs.eddy.hcen(tind)./runs.eddy.Lgauss(t0)*-1);
        legend(hplt, '1 - erf(|z|/Lz)', 'Location', 'SouthEast');
        ylabel('z/L_z'); xlabel('U(z)/U(0)'); title(runs.name);

        %runs.animate_field('u', [], it(end), 1);
        %plot(runs.rgrid.x_rho(1,ix(end))/1000, ...
        %     runs.rgrid.y_rho(iy(end),1)/1000, 'kx');
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

    function [V0, L0, Lz0] = EddyScalesForFlux(runs, tstart, tend, nsmth)

        if ~exist('tstart', 'var')
            error('EddyScalesForFlux: provide at least one time instant.');
        end
        if ~exist('tend', 'var'), tend = tstart; end
        if ~exist('nsmth', 'var'), nsmth = 10; end

        % peak velocity in eddy - downscale fit by 2.3
        %V = smooth(hypot(runs.eddy.fitx.V0, runs.eddy.fity.V0), nsmth) / 2.3 / sqrt(2);
        V = smooth(runs.eddy.fitx.V0, nsmth) / sqrt(2*exp(1));
        %V = smooth(runs.eddy.rhovor.Vke, nsmth) / 2.3 / sqrt(2);
        % radius to maximum velocity
        %L = smooth(hypot(runs.eddy.fitx.Lrho, runs.eddy.fity.Lrho), nsmth) / sqrt(2);
        L = smooth(runs.eddy.fitx.L, nsmth);
        %L = smooth(runs.eddy.rhovor.dia, nsmth) / 2;
        % gaussian decay scale in vertical
        Lz = smooth(runs.eddy.Lgauss, nsmth);

        V0 = nanmedian(addnan(V(tstart:tend), 1));
        L0 = nanmedian(addnan(L(tstart:tend), 5e5));
        Lz0 = nanmedian(addnan(Lz(tstart:tend), 1000));

        %if strcmpi(runs.name, 'ew-8352-2') | strcmpi(runs.name, 'ew-8392')
        %    L0 = runs.eddy.fitx.Lrho(tstart);
        %end
    end

    function [fluxscl] = eddyfluxscale(runs)

    % integrate velocity profile from surface to z=-H
        tind = runs.csflux.maxloc;
        syms Vint z x;
        Lz = runs.eddy.Lgauss(1);
        V0 = runs.eddy.V();
        Lx = runs.eddy.vor.dia(1)/2;

        fluxscl = double(V0 * int(exp(-x/Lx)^2, 0, Lx) ...
                  * int(1 - erf(z/Lz), z, -Lz, 0));
    end

    function [] = plot_vertProfile(runs, name, ix, iy, tinds, hax)

        for ii=1:length(tinds)
            [var(ii,:),~,~,zax] = dc_roms_read_data(runs.dir, name, tinds(ii), ...
                                                    {'x' ix ix; 'y' iy iy}, [], runs.rgrid);
        end

        if ~exist('hax', 'var') | isempty(hax)
            figure; hax = gca;
        end
        axes(hax);
        insertAnnotation([runs.name '.plot_vertprofile']);
        lightDarkLines(length(tinds));
        plot(bsxfun(@minus, var, var(1,:)), zax);
        ylabel('Z (m)');
        xlabel('\rho - \rho_{t = 0}');
        title(['(' num2str(ix) ',' num2str(iy) ') | tinds = ' num2str(tinds)]);
    end

    function [] = checkStreamerVertProfile(runs, isobath)
    % check if streamer vertical profile changes much with
    % changing (x,t) limits

    % t limits don't matter much
        xwidth = 40; % in grid points
        tstart = 1;
        tend = length(runs.csflux.time);

        [~,maxloc] = runs.calc_maxflux(isobath);
        tstart = maxloc-50;

        [tstart,tend] = runs.flux_tindices(runs.csflux.off.slope(:, isobath, isobath));

        vertbins = runs.csflux.vertbins(:, isobath);
        slopezt = runs.csflux.off.slopezt(:,:, isobath, isobath);

        fullzprofile = runs.csflux.off.slopewater.vertitrans(:,isobath,isobath);
        newzprofile = trapz(runs.csflux.time(tstart:tend), ...
                            slopezt(:,tstart:tend), 2);

        figure;
        plot(fullzprofile./max(fullzprofile), vertbins);
        hold on
        plot(newzprofile./max(newzprofile), vertbins);
        legend('all time', 'restricted time interval');
    end

    function [] = plot_recalculatedFlux(runs, factor, isobath)
        if ~exist('factor', 'var'), factor = 2; end
        if ~exist('isobath', 'var'), isobath = 1:length(runs.csflux.x); end

        figure;
        insertAnnotation([runs.name '.plot_recalculatedFlux(' factor ')']);
        hold on;
        if length(isobath) > 1, lightDarkLines(length(isobath)); end
        for iso=isobath
            fluxvec = runs.recalculateFlux(runs.bathy.hsb, iso, iso);
            [start, stop] = runs.flux_tindices(fluxvec);
            [maxf, maxl] = runs.calc_maxflux(fluxvec);

            subplot(211); hold on;
            hplt(iso) = plot(runs.csflux.time/86400, fluxvec);
            plot(runs.csflux.time(maxl)/86400, fluxvec(maxl), 'kx');
            plot(runs.csflux.time([start stop])/86400, fluxvec([start stop]), 'kx');

            subplot(212); hold on;
            plot(runs.csflux.time/86400, cumtrapz(runs.csflux.time, fluxvec));
        end
    end

    function [] = CalculateLzUncertainty(runs, debug)

        if ~exist('debug', 'var'), debug = 0; end
        T = runs.eddy.T;
        zT = runs.eddy.zT;

        for tt=1:size(T,1)
            [T0(tt),Lzfit(tt),z0(tt),T1(tt),~,conf(tt,:,:)] = gauss_fit(zT(tt,:)*-1, T(tt,:));
        end

        Lzfit(2,:) = conf(:,1,2);
        Lzfit(3,:) = conf(:,2,2);

        runs.eddy.Lzfit = Lzfit;
        runs.eddy.z0fit = [z0' conf(:,1,3) conf(:,2,3)];
        eddy.fithash = githash([mfilename('fullpath') '.m']);
        eddy = runs.eddy;
        save([runs.dir '/eddytrack.mat'], 'eddy');

        if debug
            figure;
            plot(Lzfit(1,:));
            hold on; plot(runs.eddy.Lgauss);
            plot(conf(:,1,2), 'k--');
            plot(conf(:,2,2), 'k--');
            title(runs.name);
            ylim([0 runs.params.eddy.depth*2]);
        end

    end

    function [cvy] = smoothCenterVelocity(runs, nsmooth, type)
        if ~exist('nsmooth', 'var') | isempty(nsmooth), nsmooth = 10; end
        if ~exist('type', 'var'), type = 'cen'; end

        % number of points to smooth over.
        npts = (nsmooth*runs.eddy.turnover/86400);

        if runs.bathy.axis == 'y'
            cen = smooth(runs.eddy.vor.cy(1:runs.eddy.tend), npts)/1000;
            mcen = runs.eddy.my(1:runs.eddy.tend)'/1000;
            %vel = smooth(runs.eddy.mvy, npts);
            % figure;
            % subplot(211)
            % plot(runs.eddy.my/1000); hold all
            % plot(cy);
            % subplot(212)
            % plot(smooth(runs.eddy.mvy, 15)); hold all
            % plot(vel);
        else
            cen = smooth(runs.eddy.vor.cx(1:runs.eddy.tend), npts)/1000;
            mcen = runs.eddy.mx(1:runs.eddy.tend)'/1000;
        end
        tvec = runs.eddy.t(1:runs.eddy.tend);

        if size(tvec,2) == 1, tvec = tvec'; end

        if strcmpi(type, 'cen')
            % centroid
            cvy = [0; diff(cen)./diff(smooth(tvec', npts))];
        else
            cvy = [0; smooth(diff(mcen)./diff(tvec'), npts)];
        end
    end

    function [meanx, meany, meant, meanh, err] = averageResistance(runs, nsmooth)

        if ~exist('nsmooth', 'var'), nsmooth = []; end

        factors = 1 - [0.55:0.05:0.85];
        N = length(factors);

        xx = nan([1 N]); yy = xx; tt = xx;
        for ff = 1:length(factors)
            try
                [xx(ff), yy(ff), tt(ff), ~] = runs.locate_resistance(nsmooth, factors(ff));
            catch
            end
        end

        xx = cut_nan(xx); yy = cut_nan(yy); tt = cut_nan(tt);

        N = length(xx);
        hh = runs.eddy.hcen(tt);

        meanx = mean(xx); meany = mean(yy); meant = ceil(mean(tt)); meanh = mean(hh);
        err(1) = conft(0.05, N-1) * std(xx)/sqrt(N);
        err(2) = conft(0.05, N-1) * std(yy)/sqrt(N);
        err(3) = conft(0.05, N-1) * std(tt)/sqrt(N);
        err(4) = conft(0.05, N-1) * std(hh)/sqrt(N);

        err = abs(err);
    end

    function [itsl, itse, tsl, tse] = getEddyCenterTimeScales(runs)
    % time indices when eddy center/southern edge cross slopebreak
        itsl = find_approx(runs.eddy.my, runs.bathy.xsl, 1);
        itse = find_approx(runs.eddy.my - runs.eddy.vor.dia(1)/2, ...
                           runs.bathy.xsl, 1);

        tsl = runs.time(itsl);
        tse = runs.time(itse);
    end

    function [] = streamerstruct(runs)
    % calculate along-shelf streamer structure on interpolated grid
        debug = 0;

        error('Deprecated. use avgStreamerVelSection instead.');
        xr = runs.rgrid.x_rho(1,:)';
        cen = runs.eddy.mx;
        % x-grid to interpolate on to
        dx = bsxfun(@minus, xr, cen);
        xmax = max(abs(dx(:)));
        xi = [-1 * xmax: 1000 : xmax]';

        runs.csflux.slopex = [];

        for src = 1:length(runs.csflux.ndloc)
            for isobath=1:length(runs.csflux.ndloc)
                % (x,t)
                matrix = runs.csflux.slopext(:,:,isobath, src);
                nt = size(matrix,2);

                %[~,nt] = runs.calc_maxflux(runs.csflux.off.slope(:,isobath));

                mati = nan([length(xi) nt]);
                for tt = [1:nt]
                    xvec = runs.rgrid.x_rho(1,2:end-1)' - cen(tt);
                    mati(:,tt) = interp1(xvec, matrix(:,tt), xi);
                end

                % integrate in time
                slopex = trapz(double(runs.csflux.time(1:nt)), repnan(mati, 0), 2);

                runs.csflux.slopex.slopexti(:,:,isobath,src) = mati;
                runs.csflux.slopex.flux(:,isobath,src) = slopex;
                % this fit makes no sense without applying offmask.
                % [y0, X, x0] = ...
                %     gauss_fit(xi, runs.csflux.slopex.flux(:,isobath,src), debug);
                % if debug, title(runs.name); end
                % runs.csflux.slopex.Lx(isobath,src) = X;
                % runs.csflux.slopex.y0(isobath,src) = y0;
                % runs.csflux.slopex.x0(isobath,src) = x0;

                % (z,t)
                matrix = runs.csflux.off.slopezt(:,:,isobath, src);
                nt = size(matrix,2);

                % integrate in time
                slopex = trapz(double(runs.csflux.time(1:nt)), matrix, 2);

                runs.csflux.slopez.flux(:,isobath,src) = slopex;
            end
        end

        runs.csflux.slopex.nt = nt;
        runs.csflux.slopex.xi = xi;
        runs.csflux.slopex.comment = ['[y0, Lx] = gauss_fit | flux = ' ...
                            'along-shelf structure(x,isobath,source) | xi = x-vector ' ...
                            'for flux | slopexti = interpolated ' ...
                            'to grid with eddy center as origin | ' ...
                            'nt = number of timesteps this has been calculated ' ...
                            'for.'];
    end

    function [] = plot_checkmaxflux(runs, isobath)

        fluxvec = runs.csflux.off.slope(:,isobath,isobath);
        figure;
        ax = subplot(211);
        [maxflux, maxloc] = runs.calc_maxflux(fluxvec);
        runs.animate_field('v', ax, maxloc, 1);
        liney([runs.bathy.xsb+runs.csflux.R ...
               runs.csflux.x(isobath)]/1000, [], 'r');

        subplot(212);
        plot(runs.csflux.time/86400, fluxvec);
        linex(runs.csflux.time(maxloc)/86400);
    end

    function [] = plot_checkvertprofile(runs)
    % Check whether zpeak in vertical profile changes much with time

        hfig = figure; hold on;
        lightDarkLines(length(runs.csflux.ix));
        for iso = 1:length(runs.csflux.ix)
            % figure; hold on;
            flux = runs.csflux.off.slope(:,iso,iso);
            tvec = runs.csflux.time;
            zvec = runs.csflux.vertbins(:,iso);

            nt = length(flux);
            [~,maxloc] = runs.calc_maxflux(flux);

            indices = round(linspace(maxloc,nt,30));

            clear zmax
            jj = 1;
            for kk=indices
                profile = trapz(tvec(1:kk), ...
                                runs.csflux.off.slopezt(:,1:kk,iso,iso), ...
                                2);
                norm = trapz(zvec, profile);
                % plot(profile./norm, zvec);

                [~,ind] = max(profile);
                zmax(jj) = zvec(ind);
                jj = jj + 1;
            end

            figure(hfig);
            plot(indices, zmax);
            %legend(cellstr(num2str(indices')));
            %title(num2str(iso));
        end
        xlabel('Time index since maxflux');
        ylabel('z-peak in instantaneous flux at isobath');
        legend(cellstr(num2str(runs.csflux.x'/1000, '%.0f')));
    end

    function [profile, zvec] = streamer_ideal_profile(runs, isobath, maxloc)
    % return idealized streamer vertical profile

        if ~exist('maxloc', 'var'), maxloc = 1; end

        hsb = runs.bathy.hsb;
        Lz = runs.eddy.Lgauss(maxloc);
        Ro = runs.eddy.Ro(1);

        a = 2;
        vscalefactor = 1;

        zvec = runs.csflux.vertbins(:,isobath);
        %dh = 2 * Ro * hsb * runs.csflux.ndloc(isobath);
        %zpeak = -(hsb + dh)/2;
        %zscl = -abs(zpeak/ abs(log(0.5))^(1/a));
        %zjoin = -abs(2 * zpeak);

        [zjoin,zpeak] = runs.predict_zpeak(isobath, 'use');
        zpeak = zjoin/3;
        zscl = zjoin;

        %profile = exp(-abs((zvec - zpeak)/zscl).^a) .* (zvec >= zjoin) ...
        %          + exp(-abs((zvec(find(zvec >= zjoin, 1)) - zpeak)/zscl).^a) * ...
        %          (1 - erf(-zvec/Lz/vscalefactor)) .* (zvec < zjoin);

        profile = exp(-abs((zvec - zpeak)/zscl).^a) .* (zvec >= zjoin) ...
                  + exp(-abs((zvec(find(zvec >= zjoin, 1)+1) - zpeak)/zscl).^a) * ...
                  (1 - erf(-zvec/Lz/vscalefactor)) .* (zvec < zjoin);

        %figure; plot(profile, zvec);
    end

    function [] = plot_maxfluxloc(runs)

        hf = figure; hold all;
        insertAnnotation([runs.name '.plot_maxfluxloc']);
        for isobath = 1:length(runs.csflux.x)
            [~,maxloc] = runs.calc_maxflux( ...
                runs.csflux.off.slope(:,isobath,isobath));
            figure(hf);
            plot(runs.csflux.ndloc(isobath), ...
                 (runs.eddy.my(maxloc) - runs.csflux.x(isobath)) ...
                 / runs.eddy.vor.dia(1) * 2, 'kx');
        end

        title(['Distance of eddy center from isobath at time instant ' ...
               'of max flux']);
        ylabel('(X_{eddy} - X_{isobath})/L_{eddy}');
        xlabel('y/R');
        beautify;
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
            tind = [t0 t0];
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

    function [] = checkFluxReference(runs)

        figure;
        plot(runs.eddy.mx/1000, runs.eddy.my/1000);
        if runs.bathy.axis == 'y'
            liney(runs.bathy.xsb/1000 + runs.csflux.R/1000);
        end
        xlabel('X (km)'); ylabel('Y (km)');
    end

    % read eddy-dye at surface and save it
    function [] = read_eddsurf(runs, t0, ntimes)
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [t0 t0];
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

    % read eddy-dye at surface and save it
    function [] = read_pvsurf(runs, t0, ntimes)
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [t0 t0];
        else
            tind = [1 length(runs.time)];
        end

        sz = [size(runs.rgrid.x_rho') length(runs.time)];
        if isempty(runs.pvsurf) | (t0 > size(runs.pvsurf,3))
            runs.pvsurf(:,:,tind(1):tind(2)) = ...
                dc_roms_read_data(runs.dir, 'pv', [tind], ...
                                  {'z' runs.rgrid.N-2 runs.rgrid.N-2}, ...
                                  [], runs.rgrid, 'his', 'single');
        end

        runs.pvsurf = real((runs.pvsurf));
    end

    function [] = read_csdsurf(runs, t0, ntimes)
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        % read zeta
        if ntimes == 1
            tind = [1 1]*t0;
        else
            tind = [1 length(runs.time)];
        end

        try
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
        catch ME
            runs.csdsurf = [];
            runs.read_csdsurf(t0, ntimes);
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

        if ischar(times), nt = 1; else nt = length(times); end
        opt = 'zx';
        velname = opt(1);
        axname = opt(2); % name of axis for plot;

        figure;
        subplot(2,1,1);
        cmap = brighten(cbrewer('seq','Greys',length(times)+3),0);
        cmap = cmap(3:end,:,:); % chuck out lightest colors
        hold all;

        for ii=1:nt
            if strcmpi(times, 'max flux')
                [~,tind] = runs.calc_maxflux;
                times = runs.eddy.t(tind);
            else
                tind = find_approx(runs.time/86400, times(ii), 1);
            end

            if axname == 'y'
                ref = runs.eddy.my(tind); % 0 for x-axis
                loc = runs.eddy.mx(tind);
                locax = 'x';
            else
                ref = runs.eddy.mx(tind);
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
            vvec = vvec - vvec(1);
            % plot
            if velname == 'z'
                hplt(ii) = plot(xvec, vvec, '-', 'Color', cmap(ii,:));
            else
                hplt(ii) = plot(xvec, (vvec), '-', 'Color', cmap(ii,:));
            end

            names{ii} = [num2str(times(ii)) ' | ' ...
                         num2str(vvec(runs.bathy.isb), '%0.2f') ' ' units];
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
        legend(hplt, names, 'Location', 'NorthWest');
        ylabel([velname ' ./ max(' velname ')'])
        xlabel([upper(axname) ' Distance from center / (radius)']);
        title([velname ' | ' runs.name]);
        beautify;

        subplot(2,1,2)
        plot(runs.csflux.time/86400, runs.csflux.off.shelf(:,1));
        limy = ylim; ylim([0 limy(2)]); linex(times);
        ylabel('Flux'); xlabel('Time (days)');
        beautify;
    end

    % read surface velocities for animate_pt & surf vorticity plot
    function [] = read_velsurf(runs, t0, ntimes)
        start = [1 1 runs.rgrid.N 1];
        count = [Inf Inf 1 Inf];
        stride = [1 1 1 1];

        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        if ntimes == 1
            tind = [t0 t0];
        else
            tind = [1 length(runs.time)];
        end

        if isempty(runs.usurf) | (t0 > size(runs.usurf,3)) | ...
                0 ... %any(isnan(fillnan(runs.usurf(:,:,tind(1):tind(2)),0)))
            disp('Reading surface velocity fields...');
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
    function [] = read_velbot(runs, t0, ntimes)
        start = [1 1 1 1];
        count = [Inf Inf 1 Inf];
        stride = [1 1 1 1];

        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = Inf; end

        if ntimes == 1
            tind = [t0 t0];
        else
            tind = [1 length(runs.time)];
        end

        if isempty(runs.ubot) | (t0 > size(runs.ubot,3)) | ...
                0 ... %any(isnan(fillnan(runs.usurf(:,:,tind(1):tind(2)),0)))
            disp('Reading bottom velocity fields...');
            if runs.givenFile
                runs.ubot = double(squeeze(ncread(runs.out_file, ....
                                                   'u',start,count,stride)));
            else
                runs.ubot(:,:,tind(1):tind(2)) = ...
                    dc_roms_read_data(runs.dir,'u', ...
                                      [tind],{'z' 1 1}, ...
                                      [],runs.rgrid, 'his', 'single');
            end
            if runs.givenFile
                runs.vbot = double(squeeze(ncread(runs.out_file, ....
                                               'v',start,count,stride)));
            else
                runs.vbot(:,:,tind(1):tind(2)) = ...
                    dc_roms_read_data(runs.dir,'v', ...
                                      [tind],{'z' 1 1}, ...
                                      [],runs.rgrid, 'his', 'single');
            end
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

    % ubar cross-isobath scale
    function [] = CalcUbarCrossIsobathScale(runs)
        runs.read_velbar;

        ix = runs.spng.sx2-40;
        isb = runs.bathy.isb;

        [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1));

        ubarscale.time = runs.time;
        ubarscale.scale = nan(size(runs.time));
        ubarscale.i0 = nan(size(runs.time));
        ubarscale.ind = nan(size(runs.time));
        ubarscale.minind = nan(size(runs.time));

        profiles = nan([isb+10 length(runs.time)]);

        yvec = runs.rgrid.y_rho(:,1);

        for tt=start:stop
            prof = runs.ubar(ix,1:isb+10,tt);

            [minval,minind] = min(prof);
            i0 = find(prof(1:minind) > 0, 1, 'last');
            ind = i0-1 + find(prof(i0:minind) > 0.33*minval, 1, 'last');
            %ind = find_approx(prof(1:minind), 0.33*minval, 1);

            try
                ubarscale.scale(tt) = runs.bathy.xsb - yvec(ind);
                ubarscale.i0(tt) = i0;
                ubarscale.ind(tt) = ind;
                ubarscale.minind(tt) = minind;
            catch ME
            end
            profiles(:,tt) = prof;
        end

        hash = githash([mfilename('fullpath') '.m']);
        ubarscale.hash = hash;
        ubarscale.ix = ix;
        ubarscale.profiles = profiles;

        runs.ubarscale = ubarscale;

        save([runs.dir '/ubarscale.mat'], 'ubarscale');
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

    function [vscale] = predict_verticalscale(runs, debug)
        if ~exist('debug', 'var'), debug = 0; end

        % Pedlosky pp. 413 (6.15.29)
        phys = runs.params.phys;
        D = runs.bathy.hsb;
        f0 = phys.f0;

        beta_t = -f0/D * runs.bathy.sl_shelf;
        L = runs.eddy.rhovor.dia(1)/2;
        K = 2*pi/L;
        S = phys.N2 * D^2/f0^2/L^2;

        %mu = D * K * sqrt(abs(beta_t * S/(beta_t+phys.beta)));
        mu = (2*pi/L) * sqrt(phys.N2)/f0 * sqrt(abs(beta_t/(phys.beta+beta_t)));
        vscale = 2*pi/mu;

        if debug
            zvec = runs.rgrid.z_r(:,runs.bathy.isb,1);
            plot(cosh(mu * (zvec)/D), zvec);
        end
    end

     function [] = plot_csfluxes(runs, iso)

         error('deprecated!');
         runs.streamerstruct;

         R = runs.eddy.vor.dia(1)/2;
         flux = runs.csflux.slopex.flux(:,iso,iso);

         figure;
         insertAnnotation([runs.name '.plot_csfluxes']);
         hold on;
         plot(runs.csflux.slopex.xi/R, ...
              bsxfun(@rdivide, flux, max(flux, [], 1)));
         legend(cellstr(num2str(runs.csflux.ndloc', '%0.2f')));
         xlabel('(X - X_0)/L_x');
         ylabel('Time integrated Flux / Max');
         linex(0); liney(0);
         beautify;
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
                 runs.csflux.off.shelf(:,1)/1e6);
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
            packfig(2,1,'rows');
            legend('z-centroid','zdye-centroid', ...
                'H_{center}/2','f*dia/N','vertical (Gaussian) scale');
            %end
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
            disp(ME);
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
%             permute(runs.eddy.vor.mask(ixmin:ixmax,iymin:iymax,:),[1 2 4 3])), ...
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
        strmask = reshape(full(runs.streamer.off.mask), runs.streamer.sz4dfull);

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

    function [htext] = add_isobathlabel(runs, isobath)
        htext = text(0.68,0.10, ...
                     ['y/R = ' num2str(runs.csflux.ndloc(isobath), '%.2f')], ...
                     'Units', 'normalized');
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

    function [V0, Lx, x0] = fit_vel(runs, tind)

        if runs.bathy.axis == 'y'
            loc = num2str(runs.eddy.my(tind));
            cen = runs.eddy.mx(tind);
            xvec = runs.rgrid.x_rho(1,:);
            vstr = 'v';
        else
            loc = num2str(runs.eddy.mx(tind));
            cen = runs.eddy.my(tind);
            xvec = run.rgrid.y_rho(:,1);
            vstr = 'u';
        end

        vvec = dc_roms_read_data(runs.dir, vstr, tind, ...
                                 {runs.bathy.axis loc loc; ...
                            'z' runs.rgrid.N runs.rgrid.N}, ...
                                 [], runs.rgrid, 'his', 'double');

        [V0, Lx, x0] = gauss_der_fit(xvec-cen, vvec', 0);
    end

    function [] = ideal_flux_profile(runs)

        index = 2;
        hsb = runs.bathy.hsb;
        xsb = runs.bathy.xsb;
        alpha = runs.bathy.sl_slope;
        L =  runs.eddy.vor.dia(1)/2;
        Lz = runs.eddy.Lgauss(1);

        profile = ...
            runs.csflux.off.slopewater.vertitrans(:,index);
        vertbins = runs.csflux.vertbins(:,index);
        zvec = vertbins./ max(abs(vertbins));

        yt = runs.rgrid.y_rho(vecfind(runs.bathy.h(1,:), - ...
                                      vertbins),1) / L;
        y0 = runs.csflux.x(index)/L;

        plot((erf(y0) - erf(yt)) .* (1 - erf(vertbins/Lz)), vertbins);

    end

    %% animation functions
    function [] = animate_zeta(runs, hax, t0, ntimes)
        if ~exist('hax', 'var'), hax = []; end
        if ~exist('t0', 'var'), t0 = 1; end
        if ~exist('ntimes', 'var'), ntimes = length(runs.time); end

        runs.animate_field('zeta', hax, t0, ntimes);
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

    function [] = plot_streamerwidth(runs)
        cxi = runs.eddy.vor.ee;
        offmask = bsxfun(@times, ...
                          bsxfun(@lt, runs.eddy.xr(:,1), cxi), ...
                          ~runs.sponge(2:end-1,1));
        onmask = bsxfun(@times, 1 - offmask, ...
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

    function [] = animate_floats(runs,type)
    %runs.read_zeta;
        runs.read_csdsurf;
        if strcmpi(type,'ltrans')
            runs.ltrans.animate(runs.rgrid,runs.zeta,runs.eddy);
        end
        if strcmpi(type,'roms')
            runs.roms.animate(runs.rgrid,runs.csdsurf,runs.eddy);
        end
    end


    %% colors and colormaps
    function [colors] = eddyeColormap(runs)
        colors = brighten(cbrewer('seq', 'Reds', 20), 0.4);
        colors = colors(1:end-3,:);
    end
    function [color] = shelfSlopeColor(runs, type)
        if ~exist('type', 'var') | isempty(type)
            type = 'light';
        end
        % color = [1 1 1]*0.55;
        if strcmpi(type, 'light')
            color = [49,130,189]/256;
        else
            color = [8 48 107]/256;
        end
    end
    function [color] = rhoContColor(runs)
    % color = 'k';
        color = [44 162 95]/256;
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
        else
            if strcmpi(plottype,'contourf') || strcmpi(plottype,'contour')
                eval(['[cc,hplot] = ' plottype ...
                      '(runs.rgrid.xr(' range ')/1000,' ...
                      'runs.rgrid.yr(' range ')/1000,'...
                      'double(runs.' varname '(' range ',tt)));']);

                hplot.LevelList = linspace(crange(1),crange(2), 10);
            end
        end

        shading interp;

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

    function [htext] = add_timelabel(runs, tt)
       htext = text(0.05,0.9, ...
                     ['t = ' num2str(runs.time(tt)/86400, '%.0f') ' days'], ...
                     'Units', 'normalized');
   end
   function [] = update_timelabel(runs, htext, tt)
       htext.String = ['t = ' num2str(runs.time(tt)/86400) ' days'];
    end


    function [hplot] = plot_eddy_contour(runs,plottype,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        hold on;
        [~,hplot] = contour(runs.eddy.xr(ix,iy)/1000,runs.eddy.yr(ix,iy)/1000, ...
                            runs.eddy.vor.mask(ix,iy,tt), [1 1], ...
                            'Color','k','LineWidth',1);
    end

    function update_eddy_contour(runs,handle,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        for ii=1:length(handle)
            try
                set(handle(ii),'ZData', ...
                               double(runs.eddy.vor.mask(ix,iy,tt)));
            catch ME
                set(handle(ii),'CData', ...
                               double(runs.eddy.vor.mask(ix,iy,tt)));
            end
        end
    end


    function [hplot] = plot_rho_contour(runs,plottype,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        hold on;
        [~,hplot] = contour(runs.eddy.xr(ix,iy)/1000,runs.eddy.yr(ix,iy)/1000, ...
                            runs.eddy.rhovor.mask(ix,iy,tt), [1 1], ...
                            'Color', runs.rhoContColor,'LineWidth',2);
    end
    function update_rho_contour(runs,handle,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        for ii=1:length(handle)
            try
                set(handle(ii),'ZData', ...
                               double(runs.eddy.rhovor.mask(ix,iy,tt)));
            catch ME
                set(handle(ii),'CData', ...
                               double(runs.eddy.rhovor.mask(ix,iy,tt)));
            end
        end
    end


    function [hplot] = plot_eddy_sshcontour(runs,plottype,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
        iy = max([runs.spng.sy1:runs.spng.sy2]-1,1);

        hold on;
        [~,hplot] = contour(runs.eddy.xr(ix,iy)/1000, ...
                            runs.eddy.yr(ix,iy)/1000, ...
                            runs.eddy.mask(ix,iy,tt), ...
                            [1 1], 'Color','k','LineWidth',1);
    end

    function update_eddy_sshcontour(runs,handle,tt)
        ix = max([runs.spng.sx1:runs.spng.sx2]-1,1);
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
                    num2str(runs.time(tt)/86400)  ' days']);
    end
    function update_title(runs,ht,titlestr,tt)
        set(ht,'String',[titlestr ' | ' runs.name ' | ' ...
                         num2str(runs.time(tt)/86400)  ' days']);
    end

    function [hplot] = plot_bathy(runs,plottype,color)
        ix = runs.spng.sx1:runs.spng.sx2;
        iy = runs.spng.sy1:runs.spng.sy2;

        hold on;
        if ~exist('color','var') | isempty(color), color = [1 1 1]*0.45; end
        if strcmpi(plottype,'contour')
            [hplot{4},hplot{1}] = contour(runs.rgrid.xr(ix,iy)/1000,...
                                 runs.rgrid.yr(ix,iy)/1000, ...
                                 runs.rgrid.h(iy,ix)',[200 500 1000 1500 ...
                                2000], 'Color', color);
            clabel(hplot{4}, hplot{1}, 'LabelSpacing', 108*2.75, 'Color', color);
            hax = gca;
            sbslcolor = color;
            if runs.bathy.axis == 'y'
                hplot{2} = liney(runs.bathy.xsb/1000,[], sbslcolor);
                hplot{3} = liney(runs.bathy.xsl/1000,[], sbslcolor);
            else
                hplot{2} = linex(runs.bathy.xsb/1000,[], sbslcolor);
                hplot{3} = linex(runs.bathy.xsl/1000,[], sbslcolor);
            end

            hplot{2}.LineStyle = '--';
            hplot{2}.LineWidth = 2;
            hplot{3}.LineWidth = 2;
            hplot{2}.Tag = '';
            hplot{3}.Tag = '';
            hplot{3}.LineStyle = '--';
        end
        if strcmpi(plottype,'pcolor')
            hplot{1} = pcolorcen(runs.rgrid.xr(ix,iy)/1000,...
                                 runs.rgrid.yr(ix,iy)/1000, ...
                                 runs.rgrid.h(iy,ix)');
        end
    end

   %% video functions
    function [] = video_init(runs,filename)
        if runs.makeVideo
            runs.makeVideo
            runs.mm_instance = mm_setup('frameDir',['videos/' runs.name '-' filename]);
            runs.mm_instance.pixelSize = [1600 900];
            runs.mm_instance.outputFile = ['videos/' runs.name '-' filename '.mp4'];
            runs.mm_instance.ffmpegArgs = '-q:v 1 -g 5';
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