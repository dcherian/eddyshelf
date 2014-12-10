% Creates initial condition and grid files for ROMS runs
% Modified from Arango's d_initial.m

%% Parameters

% names cannot start with number
[~,machine] = system('hostname');
if strfind(machine,'scylla')
    FOLDER    = '/scylla-a/home/dcherian/ROMS/runs/eddyshelf/topoeddy/run-2/';
    prefix    = 'tes';
    addpath(genpath('/scylla-a/home/dcherian/tools/'));
end
if strfind(machine,'poison')
    FOLDER    = '/home/poison/deepak/ROMS/runs/eddyshelf/';
    prefix    = 'tes';
end
if strfind(machine,'kadal')
    FOLDER = '/media/data/Work/eddyshelf/runs/';
    prefix    = 'tek';
end
if strfind(machine,'login')
    FOLDER = '/mit/dcherian/ROMS/runs/eddyshelf/topoeddy/run-2/';
    prefix    = 'tea';
    addpath(genpath('/mit/dcherian/tools/'));
end

fprintf('\n Writing to %s. ', [FOLDER '/' prefix '_*.nc']);
input('Are you sure?');

GRID_NAME = [prefix '_grd'];
INI_NAME  = [prefix '_ini'];
BRY_NAME  = [prefix '_bry'];
FRC_NAME  = [prefix '_frc'];

% fix file names
GRID_NAME = [FOLDER GRID_NAME '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER GRID_NAME '.nc'];
INI_NAME  = [FOLDER INI_NAME  '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER INI_NAME '.nc'];
BRY_NAME  = [FOLDER BRY_NAME  '.nc'];
FRC_NAME  = [FOLDER FRC_NAME  '.nc'];

% read parameters from file - cheap & easy way
roms_create_params

%% Create Junk IC & Grid Files + set vars to zero

% Set other variables
S.hc = S.Tcline;
S.ncname = INI_NAME;
OA_INTERPOLATE = 0;
INIname = INI_NAME;
GRDname = GRID_NAME;

% Create *new* grid file
c_grid(S.Lm+2,S.Mm+2,GRID_NAME,1);

% Create IC file
[~]=c_initial(S);

%  Set attributes for "ocean_time".

avalue='seconds since 0001-01-01 00:00:00';
[~]=nc_attadd(INIname,'units',avalue,'ocean_time');

avalue='360.0 days in every year';
[~]=nc_attadd(INIname,'calendar',avalue,'ocean_time');


%---------------------------------------------------------------------------
%  Set grid variables.
%---------------------------------------------------------------------------

V=nc_vnames(GRDname);
nvars=length(V.Variables);
%  Horizontal grid variables. Read in for input GRID NetCDF file.

% if (S.spherical),
%   S.lon_rho = nc_read(GRDname, 'lon_rho');
%   S.lat_rho = nc_read(GRDname, 'lat_rho');
%
%   S.lon_u   = nc_read(GRDname, 'lon_u');
%   S.lat_u   = nc_read(GRDname, 'lat_u');
%
%   S.lon_v   = nc_read(GRDname, 'lon_v');
%   S.lat_v   = nc_read(GRDname, 'lat_v');
% else
%   S.x_rho   = nc_read(GRDname, 'x_rho');
%   S.y_rho   = nc_read(GRDname, 'y_rho');
%
%   S.x_u     = nc_read(GRDname, 'x_u');
%   S.y_u     = nc_read(GRDname, 'y_u');
%
%   S.x_v     = nc_read(GRDname, 'x_v');
%   S.y_v     = nc_read(GRDname, 'y_v');
% end

%  Read in Land/Sea mask, if appropriate.

% for n=1:nvars,
%   name=char(V.Variables(n).Name);
%   switch (name),
%     case 'mask_rho'
%       S.mask_rho = nc_read(GRDname, 'mask_rho');
%     case 'mask_u'
%       S.mask_u   = nc_read(GRDname, 'mask_u');
%     case 'mask_v'
%       S.mask_v   = nc_read(GRDname, 'mask_v');
%   end,
% end,


%  Bathymetry.

% S.h = nc_read(GRDname, 'h');

%  Set vertical grid variables.
report = 0;
[S.s_rho, S.Cs_r]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N, ...
			     0, 0);

[S.s_w,   S.Cs_w]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N, ...
			     1, 0);

%---------------------------------------------------------------------------
%  Set zero initial conditions.
%---------------------------------------------------------------------------

Lr = S.Lm+2;   Lu = Lr-1;   Lv = Lr;
Mr = S.Mm+2;   Mu = Mr;     Mv = Mr-1;

S.zeta = zeros([Lr Mr]);
S.ubar = zeros([Lu Mu]);
S.vbar = zeros([Lv Mv]);
S.u    = zeros([Lu Mu S.N]);
S.v    = zeros([Lv Mv S.N]);
S.temp = zeros([Lr Mr S.N]);
S.salt = zeros([Lr Mr S.N]);

%  If Land/Sea masking arrays are not found, initialize them to unity.

if (~isfield(S, 'mask_rho')),  S.mask_rho = ones([Lr Mr]);  end,
if (~isfield(S, 'mask_u'  )),  S.mask_u   = ones([Lu Mu]);  end,
if (~isfield(S, 'mask_v'  )),  S.mask_v   = ones([Lv Mv]);  end,

% salt
S.salt = S0*ones(size(S.salt));

% temperature
S.temp = T0*ones(size(S.temp));

% Fix grids for future use
xmid = ceil(S.Lm/2);
ymid = ceil(S.Mm/2);
zmid = ceil(S.N/2);

if ~flags.telescoping

    S.uniform = 1;

    S.x_rho = repmat([-dx0/2:dx0:X+dx0/2]',[1 S.Mm+2]);
    S.y_rho = repmat(-dy0/2:dy0:Y+dy0/2 ,[S.Lm+2 1]);

    S.x_u = repmat([0:dx0:X]',[1 S.Mm+2]);
    S.y_u = repmat(-dy0/2:dy0:Y+dy0/2,[S.Lm+1 1]);

    S.x_v = repmat((-dx0/2:dx0:X+dx0/2)',[1 S.Mm+1]);
    S.y_v = repmat(0:dy0:Y,[S.Lm+2 1]);

    S.x_psi = repmat([0:dx0:X]',[1 S.Mm+1]);
    S.y_psi = repmat([0:dy0:Y],[S.Lm+1 1]);


    S.pm = 1./dx0 * ones(size(S.x_rho));
    S.pn = 1./dy0 * ones(size(S.x_rho));
else
    S.uniform = 0;

    ddx = grid.dxmax - grid.dxmin;
    ddy = grid.dymax - grid.dymin;

    ixr = 1:S.Lm+2; iyr = 1:S.Mm+2;
    clear dx dy
    dx = grid.dxmin + ( (ixr <= grid.ixp & ixr >= grid.ixn) ...
        + ddx * (ixr > grid.ixp) .* tanh( (ixr-grid.ixp)/grid.xscalep ) ...
        - ddx * (ixr < grid.ixn) .* tanh( (ixr-grid.ixn)/grid.xscalen ));
    dy = grid.dymin + ( (iyr <= grid.iyp & iyr >= grid.iyn) ...
        + ddy * (iyr > grid.iyp) .* tanh( (iyr-grid.iyp)/grid.yscalep ) ...
        - ddy * (iyr < grid.iyn) .* tanh( (iyr-grid.iyn)/grid.yscalen ));

    % pm,pn are dx,dy now. Inverted later
    dx = repmat(dx',[1 S.Mm+2]);
    dy = repmat(dy ,[S.Lm+2 1]);

    S.x_rho = -3*dx(1,1)/2 + cumsum(dx,1);
    S.y_rho = -3*dy(1,1)/2 + cumsum(dy,2);

    S.x_u = avg1(S.x_rho,1); S.y_u = avg1(S.y_rho,1);
    S.x_v = avg1(S.x_rho,2); S.y_v = avg1(S.y_rho,2);
    S.x_psi = avg1(avg1(S.x_rho,1),2); S.y_psi = avg1(avg1(S.y_rho,1),2);
end

S.xl = S.x_v(end,1);
S.el = S.y_u(1,end);

X = S.xl;
Y = S.el;

% Calculate grid metrics
[S.pm,S.pn,S.dndx,S.dmde] = grid_metrics(S,false);
S.angle = zeros(size(S.x_rho));

% matrices for future use
xrmat = repmat(S.x_rho,[1 1 S.N]);
yrmat = repmat(S.y_rho,[1 1 S.N]);
xumat = repmat(S.x_u,[1 1 S.N]);
yumat = repmat(S.y_u,[1 1 S.N]);
xvmat = repmat(S.x_v,[1 1 S.N]);
yvmat = repmat(S.y_v,[1 1 S.N]);

% set lat_* and lon_* just in case
if ~S.spherical
    S.lon_rho = S.x_rho; S.lat_rho = S.y_rho;
    S.lon_u = S.x_u; S.lat_u = S.y_u;
    S.lon_v = S.x_v; S.lat_v = S.y_v;
    S.lon_psi = S.x_psi; S.lat_psi = S.y_psi;
end

% change axes to km if needed
fx = 1; fy = 1; lx = '(m)'; ly = '(m)';
if max(abs(S.x_u(:))) > 3500
    fx = 1000; lx = '(km)';
end
if max(abs(S.x_u(:))) > 3500
    fy = 1000; ly = '(km)';
end

% write parameters to initial conditions .nc file - REST AT THE END
tic;
fprintf('\n Writing params \n');
write_params_to_ini(INI_NAME,flags);
write_params_to_ini(INI_NAME,bathy);
write_params_to_ini(INI_NAME,phys);
write_params_to_ini(INI_NAME,grid);
toc;

fprintf('\n Initialization - %4.1f MB \n\n', monitor_memory_whos);

%% Bathymetry + Coriolis + more grid stuff
% Coriolis with beta. f = f0 @ y=ymid
fnew = f0*ones(size(S.x_rho));
f = fnew + beta * (S.y_rho - S.y_rho(1,ymid));
clear fnew

% make plots to check bathymetry?
bathy_plot = 0

if isnan(bathy.H_shelf)
    bathy.H_shelf = bathy.H_sbreak - bathy.sl_shelf * bathy.L_shelf;
end

if flags.flat_bottom
    S.h = Z * ones(size(S.x_rho));

    bathy.xsb = 0;
    bathy.isb = 0;
    bathy.hsb = Z;
    if bathy.axis == 'x'
        bathy.xsl = X/2;
    else
        bathy.xsl = Y/2;
    end
    bathy.hsl = Z;
else
    % linear bathymetry
    %if flags.linear_bathymetry == 1

    if flags.crooked_bathy
        [S] = bathy2_x(S,bathy,X,Y);
    else
        [S] = bathy_simple(S,bathy,X,Y,bathy.axis);
    end

    switch bathy.axis
        case 'x'
            ax_cs = xrmat(:,1,1);
            ax_as = yrmat(1,:,1)'; % cross-shelf axis
            i_cs = 1; % cross shelf axis
            i_as = 2; % along shelf axis
            hvec = S.h(:,1);
        case 'y'
            ax_cs = yrmat(1,:,1)'; % cross-shelf axis
            ax_as = xrmat(:,1,1); % along shelf axis
            i_cs = 2; % along shelf axis
            i_as = 1; % cross-shelf axis
            hvec = S.h(1,:)';
    end

    % run smoother
    for i=1:bathy.n_passes
        hvec = smooth(hvec,bathy.n_points);
    end
    if bathy_plot
        fbathy = figure;
        subplot(133);
        plot(hvec,'b'); hold on
        plot(hvec,'k');
    end

    % smooth transition to deep water even more
    % find end of slope
    dh2dx2 = diff(hvec,2,1)./avg1(diff(ax_cs,1,1).^2,1);
    [~,isl] = min(dh2dx2(:));
    isl = isl - bathy.n_points;
    % smooth again!
%    for i=1:bathy.n_passes
%       hvec(isl:end) = smooth(hvec(isl:end),bathy.n_points*4);
%    end
%    reconstruct h
    if bathy.axis == 'y'
        S.h = repmat(hvec',[size(S.h,1) 1]);
    else
        S.h = repmat(hvec, [1 size(S.h, 2)]);
    end

    % Calculate Burger numbers
    S_sh = bathy.sl_shelf * sqrt(N2)./f0; % shelf
    S_sl = bathy.sl_slope * sqrt(N2)./f0; % slope

    % Calculate topographic beta
    b_sh = f0 * bathy.sl_shelf / bathy.H_shelf;
    b_sl = f0 * bathy.sl_slope / max(S.h(:));

    % calculate for smoothed bathymetry
    % find shelfbreak
    dh2dx2 = diff(hvec,2,1)./avg1(diff(ax_cs,1,1).^2,1);
    [~,bathy.isb] = max(dh2dx2(:));
    bathy.isb = bathy.isb-1;
    bathy.hsb = hvec(bathy.isb);
    bathy.xsb = ax_cs(bathy.isb);

    % find end of slope
    [~,bathy.isl] = min(dh2dx2(:));
    bathy.isl = bathy.isl;
    bathy.hsl = hvec(bathy.isl);
    bathy.xsl = ax_cs(bathy.isl);
end

if any(S.h(:) < 0), error('h < 0!'); end

% Land - Sea Mask
S.mask_u = ones(size(S.x_u));
S.mask_v = ones(size(S.x_v));
S.mask_rho = ones(size(S.x_rho));
S.mask_psi = ones(size(S.x_psi));

% duplicate in lat, lon fields
S.lat_rho = S.y_rho;
S.lon_rho = S.x_rho;
S.lat_u = S.y_u;
S.lon_u = S.x_u;
S.lat_v = S.y_v;
S.lon_v = S.x_v;
S.lat_psi = S.y_psi;
S.lon_psi = S.x_psi;

% z grids - only z_r is available in ana_grid.h
[z_r]=set_depth(S.Vtransform, S.Vstretching, ...
                S.theta_s, S.theta_b, S.hc, S.N, ...
                1, S.h, S.zeta,0);
%[zrflat]=set_depth(S.Vtransform, S.Vstretching, ...
%    S.theta_s, S.theta_b, S.hc, S.N, ...
%    1, max(S.h(:)).*ones(size(S.h)), S.zeta,0);
[z_u]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 3, S.h, S.zeta,0);
[z_v]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 4, S.h, S.zeta,0);
[z_w]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 5, S.h, S.zeta,0);
zrmat = z_r;
zumat = z_u;
zvmat = z_v;
zwmat = z_w;

clear z_r z_u z_v z_w

Hz = diff(zwmat,1,3); %bsxfun(@rdivide,diff(zwmat,1,3),permute(diff(S.s_w',1,1),[3 2 1]));
Z  = abs(max(S.h(:)));

[rx0,rx1,rx0mat,rx1mat] = stiffness(S.h,zrmat);
bathy_title = sprintf(['\n\n Beckmann & Haidvogel number (r_{x0}) = %f (< 0.2 , max 0.4) \n' ...
            ' Haney number (r_{x1}) = %f (< 9 , maybe 16)'], rx0,rx1);

if bathy_plot
    if ~exist('fbathy','var')
        fbathy = figure;
    else
        figure(fbathy);
    end
    subplot(131)
    if strcmp(bathy.axis,'x')
        ind = ymid;
        pcolor(squeeze(xrmat(:,ind,:))/fx,squeeze(zrmat(:,ind,:)), ...
               ...%        zeros(size(squeeze(zrmat(:,ind,:)))));
        squeeze(rx1mat(:,ind,:)));
        xlabel(['x ' lx]);
        linex([bathy.xsb bathy.xsl]/fx);
    else
        ind = xmid;
        pcolor(squeeze(yrmat(ind,:,:))/fy,squeeze(zrmat(ind,:,:)), ...
               ... %zeros(size(squeeze(rx1mat(ind,:,:)))));
        squeeze(rx1mat(ind,:,:)));
        xlabel(['y ' ly]);
        linex([bathy.xsb bathy.xsl]/fx);
    end
    colorbar;
    title(bathy_title);

    %subplot(132)
    %pcolor(S.x_rho/fx,S.y_rho/fy,zeros(size(S.x_rho)));
    %title(bathy_title); colorbar;
    %xlabel(['x ' lx]); ylabel(['y ' ly]); zlabel('z (m)');
    %beautify;

    %subplot(133)
    %contour(S.x_rho/fx,S.y_rho/fy,-S.h./f,30,'k');
    %xlabel(['x ' lx]); ylabel(['y ' ly]); title('f/h');
    %beautify;

    % calculate rotation angles to estimate effect of S-surface
    % momentum mixing
    % Redi (1982) notation
    % delta^2 = (dz/dx)_s ^2 + (dz/dy)_s ^2
    subplot(132);
    visc = 800;

    dzdx = abs(bsxfun(@rdivide, diff(zrmat, 1, 1), avg1(dx, 1)));
    dzdy = abs(bsxfun(@rdivide, diff(zrmat, 1, 2), avg1(dy, 2)));
    delta2 = avg1(dzdx, 2).^2 + avg1(dzdy, 1).^2;
    factor = fillnan(delta2./(1+delta2), 0);
    dzmat = avg1(avg1(diff(zwmat, 1, 3), 1), 2);
    viscz = visc*factor; %-1*addnan(-visc*factor,0);
    timescale = addnan((dzmat/2).^4 ./ (viscz) /86400, 1000);

    if bathy.axis == 'y'
        contourf(avg1(squeeze(yrmat(xmid,:,:)), 1)/1000, ...
                 avg1(squeeze(zrmat(xmid, :, :)), 1), ...
                 squeeze(factor(xmid, :,:)));
        colormap(flipud(colormap('bone')));
        colorbar;
        %    caxis([0 1]);
    else
        pcolor(squeeze(xrmat(:,ymid,:))/1000, squeeze(zrmat(:,ymid, :)), ...
               squeeze(delta2(:, ymid, :)));
    end
    title('estimated vertical visc factor');
    %spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))
end

% save for later use
bathy.h = S.h;

fprintf('\n Bathy - %4.1f MB \n\n', monitor_memory_whos);

%% Now set initial conditions - first shelfbreak front

% all variables are 0 by default S = S0; T=T0
tic

% reset initial tracers just in case
S.temp = T0*ones(size(S.temp));
S.salt = S0*ones(size(S.salt));

% setup variables for diff_cgrid used later
tgrid.zw = permute(zwmat,[3 2 1]); tgrid.s_w = S.s_w;
tgrid.xmat = xrmat; tgrid.ymat = yrmat; tgrid.zmat = zrmat;
tgrid.s = S.s_rho;

% Create background state (assumes uniform horizontal grid)
% assign initial stratification
bstrat = T0.*ones(size(zrmat));
if strat.z0 > 0 , strat.z0 = strat.z0 * -1; end
% N2 here is phys.N2 = strat.N2
if flags.conststrat
    % constant stratification
    Tz = N2/g/TCOEF * ones(size(zwmat) - [0 0 2]); % at w points except top / bottom face
else
    % non-constant stratification.
    zmat = zwmat(:,:,2:end-1);
    N2mat = N2 .* ...
            ((exp(-(zmat - strat.z0)./strat.Lp) .* (zmat >  strat.z0)) + ...
             (exp( (zmat - strat.z0)./strat.Lm) .* (zmat <= ...
                                                    strat.z0)));
    % clamp min N² to 1e-6
    N2mat(N2mat < 1e-6) = 1e-6;
    Tz = N2mat./g./TCOEF;

    Zscl = 1;
    subplot(1,2,1)
    plot(squeeze(bstrat(end,end,:)), squeeze(zrmat(end,end,:))./Zscl);
    liney(strat.z0./Zscl);
    title('Temp');
    subplot(1,2,2)
    plot(squeeze(N2mat(end,end,:)), squeeze(zmat(end,end,:))./Zscl);
    liney(strat.z0./Zscl);
    title('N^2');
end
for k=size(zrmat,3)-1:-1:1
    bstrat(:,:,k) = bstrat(:,:,k+1) - Tz(:,:,k).*(zrmat(:,:,k+1)-zrmat(:,:,k));
end

% assign background stratification to temp
S.temp = bstrat;

% choose axis + assign appropriate variables
switch bathy.axis
    case 'x'
        ax_cs = xrmat(:,1,1);
        axmat = xrmat;
        ax_as = yrmat(1,:,1)'; % cross-shelf axis
        i_cs = 1; % cross shelf axis
        i_as = 2; % along shelf axis
        zeta_sign = 1; % for free surface height thermal wind
         % uv_sign = 1; % for u/v thermal wind - SAME SIGN?
        flip_flag = 0;

        vel = 'v';

        velzmat = zvmat;
        zetahax = xvmat(:,1,end);
        zetahax2 = yrmat(:,:,end); % interpolation at edges of zeta field

    case 'y'
        ax_cs = yrmat(1,:,1)'; % cross-shelf axis
        axmat = yrmat;
        ax_as = xrmat(:,1,1); % along shelf axis
        i_cs = 2; % along shelf axis
        i_as = 1; % cross-shelf axis
        zeta_sign = -1; % for integration of free surface height
          %uv_sign = -1; % for u/v thermal wind
        flip_flag = 1;

        vel = 'u';
        velzmat = zumat;
        zetahax = yumat(1,:,end);
        zetahax2 = xrmat(:,:,end); % interpolation at edges of zeta field
end

% create front if needed
if flags.front
    % build cosine curve about shelfbreak (do each lat/lon line individually)
    % this gets gradient on a Z level
    % for this to work Tx must be independent of z!
    if strcmpi(front.Tra,'salt')
        coef = -SCOEF;
        S.Tra = S.salt;
    else
        coef = TCOEF;
        S.Tra = S.temp;
    end

    front.dT = front.dRho/abs(coef)/R0; % delta tracer across front
    % divide dT/2 since tanh goes from (-1 to 1)*dT -> 2dT across front
    x0     = bathy.xsb + (zrmat + bathy.hsb)/front.slope;
    S.Tra  = front.dT/2 * (tanh( (axmat-x0)/front.Lx)) .* exp(- (zrmat/front.Lz).^2);
    S.Trax = front.dT/2/front.Lx * sech( (axmat-x0)/front.Lx).^2 .* exp(- (zrmat/front.Lz).^2);
    S.Tra  = S.Tra - min(S.Tra(:));

    clear x x0

    % Make plots to check temperature field
    h_check = figure;
    S = reset_flip(S);
    S = flip_vars(flip_flag,S);
    axt(3) = subplot(243);
    plot(ax_cs/1000,1+S.Trax(:,1)./max(abs(S.Trax(:,1)))); hold on
    plot(ax_cs/1000,1-S.h(:,1)./max(S.h(:)),'r')
    xlabel('cross shelf axis (km)');
    legend('Tracer gradient','bathy')
    S = flip_vars(flip_flag,S);

    dsst = max(max(S.Tra(:,:,end)))-min(min(S.Tra(:,:,end)));
    ssttitle = sprintf('Surface Tra | \\Delta surface Tracer = %.2f ', mean(dsst(:)));

    % check front plots
    if bathy.axis == 'x'
        axt(1) = subplot(241);
        contourf(squeeze(xrmat(:,ymid,:))/1000,squeeze(zrmat(:,ymid,:)),squeeze(S.Tra(:,ymid,:)),20);
        xlabel('x (km)'); ylabel('z (m)'); title('Tracer (Salt?)');  colorbar
        axt(2) = subplot(242);
        contourf(squeeze(xrmat(:,ymid,:))/1000,squeeze(zrmat(:,ymid,:)),squeeze(S.Trax(:,ymid,:)),20);
        xlabel('x (km)'); ylabel('z (m)'); title('Tracer Gradient'); colorbar
    else
        axt(1) = subplot(241);
        contourf(squeeze(yrmat(xmid,:,:))/1000,squeeze(zrmat(xmid,:,:)),squeeze(S.Tra(xmid,:,:)),20);
        xlabel('y (km)'); ylabel('z (m)'); title('Tracer'); colorbar
        axt(2) = subplot(242);
        contourf(squeeze(yrmat(xmid,:,:))/1000,squeeze(zrmat(xmid,:,:)),squeeze(S.Trax(xmid,:,:)),20);
        xlabel('y (km)'); ylabel('z (m)'); title('Tracer Gradient'); colorbar
    end
    axt(4) = subplot(244);
    contourf(xrmat(:,:,end)/1000,yrmat(:,:,end)/1000, S.Tra(:,:,end),40);
    xlabel('x (km)'); ylabel('y (km)'); title(ssttitle); colorbar

    %fig;pcolorcen(squeeze(avg1(xrmat(:,ymid,:),1))/1000,squeeze(avg1(zrmat(:,ymid,:))),squeeze(S.Trax_sig(:,ymid,:)));

    % calculate velocity field
    vel_shear = zeta_sign * g*coef .* bsxfun(@rdivide,avg1(S.Trax,i_as),avg1(f,i_as));
    vmat = zeros(size(velzmat)); % zero vel. at bottom
    for ii=2:size(velzmat,3) % bottom to top
        vmat(:,:,ii) = vmat(:,:,ii-1) + vel_shear(:,:,ii).*(velzmat(:,:,ii)-velzmat(:,:,ii-1));
    end
    eval(['S.' vel '=vmat;']);

    figure(h_check)
    % plot velocity
    if bathy.axis == 'x'
        axt(5) = subplot(245);
        contourf(squeeze(xvmat(:,ymid,:))./fx,squeeze(zvmat(:,ymid,:)),squeeze(S.v(:,ymid,:)));
        colorbar; title('velocity'); xlabel('x (km)'); ylabel('z (m)');
        axt(6) = subplot(246);
        contourf(squeeze(xvmat(:,ymid,:))./fx,squeeze(zvmat(:,ymid,:)),squeeze(vel_shear(:,ymid,:)));
        colorbar; title('velocity shear'); xlabel('x (km)'); ylabel('z (m)');
    else
        axt(5) = subplot(245);
        contourf(squeeze(yumat(xmid,:,:))./fy,squeeze(zumat(xmid,:,:)),squeeze(S.u(xmid,:,:)));
        colorbar; title('velocity'); xlabel('y (km)'); ylabel('z (m)');
        axt(6) = subplot(246);
        contourf(squeeze(yumat(xmid,:,:))./fy,squeeze(zumat(xmid,:,:)),squeeze(vel_shear(xmid,:,:)));
        colorbar; title('velocity shear'); xlabel('y (km)'); ylabel('z (m)');
    end

    % then calculate zeta
    S.zeta = nan([size(S.h,1) size(S.h,2)]);
    tmp = zeta_sign * f0/g * cumtrapz(zetahax,vmat(:,:,end),i_cs);
    if flip_flag,
        S.zeta = S.zeta';
        tmp = tmp';
        zetahax2 = zetahax2';
    end
    S.zeta(:,2:end-1) = (tmp(:,1:end-1) + tmp(:,2:end))/2;
    S.zeta(:,1) = S.zeta(:,2) - (zetahax2(:,2) - zetahax2(:,1)).* ...
                (S.zeta(:,3)-S.zeta(:,2))./(zetahax2(:,3) - zetahax2(:,2));
    S.zeta(:,end) = S.zeta(:,end-1) + ...
                (zetahax2(:,end) - zetahax2(:,end-1)).* ...
                    (S.zeta(:,end-1)-S.zeta(:,end-2))./(zetahax2(:,end-1) - zetahax2(:,end-2));
    if flip_flag, % flip everything just in case
        S.zeta = S.zeta';
        tmp = tmp';
        zetahax2 = zetahax2';
    end

    % plot zeta
    figure(h_check)
    axt(7) = subplot(247);
    contourf(xrmat(:,:,end)./fx,yrmat(:,:,end)./fy,S.zeta);
    %axis image;
    colorbar; title('zeta (front only)');

    % link appropriate axes
    linkaxes([axt(1) axt(2) axt(5) axt(6)],'xy');
    linkaxes([axt(4) axt(7)],'xy');

    % if eddy then line is called after placing last graph
    if ~flags.eddy
        spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))
    end

    % calculate diagnostics
    % define horizontal and vertical scales of jet  as half the core velocity
    % like in Fratantoni et al (2001) & Linder & Gawarkiewicz (1998)
    if bathy.axis == 'x'
        vsurf = vmat(:,ymid,end);
        vvert = squeeze(vmat(:,ymid,:));
    else
        vsurf = vmat(xmid,:,end)';
        vvert = squeeze(vmat(xmid,:,:));
    end

    [vm,imax] = max(abs(vsurf));
    ind = sort(find_approx(abs(vsurf),vm/2,2));
    front.hscale = ax_cs(ind(2))-ax_cs(ind(1));
    vvert = vvert(imax,:);
    ind = find_approx(abs(vvert),vm/2,1);
    if bathy.axis == 'x'
        front.vscale = abs(zrmat(imax,ymid,ind));
    else
        front.vscale = abs(zrmat(xmid,imax,ind));
    end

    if strcmpi(front.Tra,'salt')
        S.salt = S.salt + S.Tra;
    else
        S.temp = S.temp + S.Tra;
    end
    % clear some vars
    clear S.Tra S.Trax S.Traz S.Tz vmat

    % estimate Rossby number of front
    if bathy.axis == 'x'
        vx = diff(S.v, 1, 1)./diff(xvmat, 1, 1);
        nondim.front.Ro = max(abs(vx(:)))./f0;
    else
        uy = diff(S.u, 1, 2)./diff(yumat, 1, 2);
        nondim.front.Ro = max(abs(uy(:)))./f0;
    end

    fprintf('\n Front Parameters: ');
    fprintf(['\n Ro = %.2f | h-scale = %.2f km | v-scale = %.2f m | ' ...
             'vscale/Hsb = %.2f \n\n'], ...
            nondim.front.Ro, front.hscale/1000, front.vscale, front.vscale./bathy.hsb);

    fprintf('\n Writing front params \n');

    write_params_to_ini(INI_NAME,front);
    toc;
end % if shelfbreak front

fprintf('\n Front - %4.1f MB \n\n', monitor_memory_whos);

%% Now create eddy using temperature

if flags.eddy

    % reset in case I'm debugging
    if ~flags.front
        S.temp = bstrat;
        S.u = zeros(size(S.u));
        S.v = zeros(size(S.v));
        S.zeta = zeros(size(S.zeta));
    end

    % Set eddy parameters that depend on something else
    if flags.conststrat
        eddy.Ldef = sqrt(phys.N2)*Z/pi/f0; % deformation radius NH/pi/f
                                           % eddy.dia =
                                           % 2*bathy.L_slope;
    else
        % determine deformation radius in deepest water
        % this works for NS and EW isobaths
        N2vec = squeeze(N2mat(end, end, :));
        zvec  = squeeze(zrmat(end, end, :));

        % calculate vertical modes
        [Vmode, Hmode, c] = vertmode(N2vec, zvec, 1, 0);
        eddy.Ldef = c(1)./f0;

        clear Tz N2mat
    end

    % eddy.Bu is (Ldef / eddy_radius)^2
    if isnan(eddy.dia)
        eddy.dia = 2 * 1./sqrt(eddy.Bu) * eddy.Ldef;
    end

    % set depth according to fL/N
    if (flags.vprof_gaussian | flags.pres_gaussian) & isnan(eddy.depth)
        eddy.depth = eddy.depth_factor * f0 * eddy.dia/2 / ...
                sqrt(N2);
    end

    % check for consistency, just in case
    if ~flags.flat_bottom
        factor = 1;1/sqrt(eddy.Bu) * eddy.dia/2 / bathy.L_slope * pi/ ...
                 bathy.S_sl - sqrt(N2)*bathy.hsb/pi/f0/eddy.Ldef;
    else
        factor = 1;
    end

    if factor > 1.05 || factor < 0.95
        error([' pi/S_sl * Le/Lsl * 1/sqrt(Bu) = ' ...
                num2str(factor)]);
    else
        warning([' pi/S_sl * Le/Lsl * 1/sqrt(Bu) = ' ...
                num2str(factor)]);
    end
    if flags.eddy_zhang
        xtra = (4.3)*eddy.dia/2;
    else
        if eddy.a == 3
            xtra = (2.3)*eddy.dia/2;
        else
            xtra = 2.8 * eddy.dia/2;
        end
    end

    if isnan(eddy.buffer)
        %start eddy 1 deformation radius away
        %from shelfbreak
        eddy.buffer = eddy.Ldef;
    end

    switch bathy.axis % cross-shore axis
        case 'x'
            if isnan(eddy.cx)
                eddy.cx = bathy.xsl+eddy.buffer+xtra;
            end
            if isnan(eddy.cy)
                % add deformation radius buffer away from boundary
                % note there is no sponge at the inflow boundary
                eddy.cy = 2*eddy.Ldef+xtra+eddy.buffer_sp;
            end
            fprintf('Distance from eastern edge = %.2f km \n', ...
                    (X-eddy.cx-xtra)/1000);
            if X-eddy.cx-xtra < eddy.buffer_sp
                error('Eastern edge - eddy edge < buffer');
            end
        case 'y'
            if isnan(eddy.cx)
                if ~flags.ubt_initial
                    % add deformation radius buffer away from boundary
                    % note there is no sponge at the inflow boundary
                    eddy.cx = X-eddy.Ldef-xtra-eddy.buffer_sp; % center of eddy
                else
                    eddy.cx = eddy.buffer_sp + eddy.Ldef+xtra;
                end
            end
            if isnan(eddy.cy)
                if flags.flat_bottom
                    eddy.cy = Y/5 + xtra;
                else
                    base = S.y_rho(1, find_approx(bathy.h(1,:), ...
                                                  1.5*eddy.depth, 1));
                    eddy.cy = base + eddy.buffer+xtra; %597000;
                end
            end
            fprintf('Distance from northern edge = %.2f km \n', ...
                    (Y-eddy.cy-xtra)/1000);
            if Y-eddy.cy-xtra < eddy.buffer_sp
                error('Northern edge - eddy edge < sponge buffer');
            end
    end

    % cylindrical co-ordinates
    r0 = eddy.dia/2; % This is required to account for exponential decay
    [th,r] = cart2pol((S.x_rho-eddy.cx),(S.y_rho-eddy.cy));
    rnorm  = r./r0; % normalized radius
    eddy.ix = find_approx(S.x_rho(:,1),eddy.cx,1);
    eddy.iy = find_approx(S.y_rho(1,:),eddy.cy,1);

    % assume that eddy location is in deep water
    eddy.z = squeeze(zrmat(eddy.ix,eddy.iy,:));
    eddy.xyprof = nan(size(rnorm));
    % eddy temp. field - xy profile - Reference: Katsman et al. (2003)
    if flags.solidbody_katsman % solid body rotation
        if eddy.a <= 2
            error('eddy.a > 2 for solid body profiles! See Katsman et al. (2003)');
        end
        % rnorm = rnorm .* 2;
        % normalized rm = rm_tilde in Katsman et al. (2003)
        % values determined by matching first and second derivatives of
        % temp. profile at rm - see paper for details
        rmnorm = nthroot( (eddy.a-2)/(eddy.a-1) , eddy.a);
        gamma = -2 * (eddy.a-2)./eddy.a ./ (rmnorm)^2;

        exponent = (eddy.a - 1)/eddy.a .* (rnorm.^(eddy.a) - rmnorm.^(eddy.a));
        eddy.xyprof = (gamma/2 * rnorm.^2 + 1) .* (rnorm <= rmnorm) ...
                       + (gamma/2 *rmnorm^2 + 1) .* exp( -1 * exponent ) .* ...
                                                           (rnorm > rmnorm);
    else
        if flags.eddy_zhang
            eddy.xyprof = (1 - rnorm.^2 /2) .* exp(-rnorm.^2/2);
        else
            %  gaussian in xy plane (normalized)
            exponent = (eddy.a - 1)/eddy.a .* (rnorm.^(eddy.a)); % needed for radial calculations later
            eddy.xyprof = exp( -1 * exponent );
        end
    end
    eddy.xyprof = eddy.xyprof./max(eddy.xyprof(:));

    % use half-Gaussian profile & normalize
    if flags.vprof_gaussian
        eddy.zprof = exp(-(eddy.z ./ eddy.depth) .^ 2);
    else
        if flags.pres_gaussian
            eddy.zprof = -eddy.z .* exp(-(eddy.z ./ ...
                                          eddy.depth).^2);
            eddy.zprof = eddy.zprof ./ max(eddy.zprof);
        else % energy in barotropic and BC1 mode
            eddy.zprof = abs((cos((- eddy.z)/Z * pi + eddy.theta0)));
            eddy.depth = Z/2;
        end
    end
    %eddy.zprof = eddy.zprof./trapz(eddy.z,eddy.zprof);

    % add eddy temperature perturbation
    eddy.tz = repmat(permute(eddy.zprof,[3 2 1]),[S.Lm+2 S.Mm+2 1]);
    eddy.temp = eddy.tamp * bsxfun(@times,eddy.xyprof,eddy.tz);
    if ~isnan(eddy.temp)
        S.temp = S.temp + eddy.temp;
    end

    % Radial
    if max(~isnan(eddy.temp(:))) % && flags.use_radial
        % integrated z profile of eddy.temp (HAS to include stratification)
        % needed for zeta calculation
%         int_Tz = nan([size(xrmat,1) size(xrmat,2)]);
%         for i=1:size(eddy.tz,1)
%             for j=1:size(eddy.tz,2)
%                 int_Tz(i,j) = trapz(squeeze(zrflat(i,j,:)),eddy.tz(i,j,:),3);
%             end
%         end

        % SSH calculation is same for gradient wind & geostrophic balance!
        % also same for all profiles
        %S.zeta = S.zeta + TCOEF*eddy.tamp * int_Tz .* eddy.xyprof;
        eddy.zeta = TCOEF * trapz(eddy.z, ...
                 bsxfun(@minus,eddy.temp,eddy.temp(eddy.ix,eddy.iy,:)),3);
        S.zeta = S.zeta + eddy.zeta;
        S.zeta = S.zeta - min(S.zeta(:));

        % Calculate azimuthal velocity shear (r d(theta)/dt)_z using geostrophic balance
        if flags.solidbody_katsman
            dTdr = gamma * rnorm ./ r0 .* (rnorm <= rmnorm) ...
                    + (gamma/2 .* rmnorm^2 + 1) .*  (-(eddy.a-1) ./ r0 .* rnorm.^(eddy.a-1)) ...
                               .* exp(-exponent) .* (rnorm > rmnorm);
        else
            if flags.eddy_zhang
                dTdr = - rnorm./r0 .* exp(-rnorm.^2/2) .* (2 - rnorm.^2/2);
%            else
                % gaussian eddy
%                 S.zeta = S.zeta + -TCOEF * eddy.tamp * int_Tz .* (1-eddy.xyprof);
%                 % Calculate azimuthal velocity shear (r d(theta)/dt)_z using geostrophic balance
%                 rutz = avg1(bsxfun(@times, eddy.temp, ...
%                         g*TCOEF* 1./f .* (-exponent./r *eddy.a)),3);
            end
        end
        % rutz = eddy.tamp *(TCOEF*g) .* bsxfun(@times,eddy.tz,dTdr./f);

        % azimuthal velocity = r d(theta)/dt
        rut = eddy.tamp * (TCOEF*g) .* bsxfun(@times,cumtrapz(eddy.z,eddy.tz,3),dTdr./f);
%         rut = zeros(size(xrmat));
%         for i=2:size(xrmat,3)
%             rut(:,:,i) = rut(:,:,i-1) + rutz(:,:,i-1).*(zrmat(:,:,i)-zrmat(:,:,i-1));
%         end

        % solve quadratic for vel. if gradient wind balance
        vgeo = rut;
        if flags.use_gradient
            rfb2 = r.*f ./ 2;
            sdisc = sqrt(1 + bsxfun(@times,vgeo,2./rfb2));% sqrt(discriminant)
            if isreal(sdisc) % gradient wind doesn't always work with anticyclones
                rut = bsxfun(@times,(-1 + sdisc), rfb2);
                warning('Using gradient wind balance.');
            else
                % cyclostrophic balance doesn't work yet
                error(['gradient wind calculated complex v! - ' ...
                         'Ro > 0.25']);
                %rut = -sqrt(bsxfun(@times,vgeo,-2*rfb2)); % need -r for real solutions
            end
        end

        eddy.u = -1 * bsxfun(@times,rut, sin(th));
        eddy.v =      bsxfun(@times,rut, cos(th));

        % check angular momentum integral
%         % if zero, eddy is isolated
%         L = R0*(1 - TCOEF * eddy.temp) .* ...
%             (bsxfun(@times,eddy.v,(S.x_rho - eddy.cx)) - ...
%              bsxfun(@times,eddy.u,(S.y_rho - eddy.cy)));

        S.u = S.u + avg1(eddy.u,1);
        S.v = S.v + avg1(eddy.v,2);

        % calculate surface vorticity
        vor = avg1(diff(eddy.v(:,:,end),1,1)./diff(xrmat(:,:,end),1,1), 2) - ...
              avg1(diff(eddy.u(:,:,end),1,2)./diff(yrmat(:,:,end),1,2), 1);

        % get vorticity mask
        vormask = vor < 0;
        % extract biggest region
        regions = bwconncomp(vormask, 8);
        nn = [];
        for ll=1:regions.NumObjects
            nn(ll) = length(regions.PixelIdxList{ll});
        end
        [~,ind] = sort(nn, 'descend');
        newmask = zeros(size(vormask));
        newmask(regions.PixelIdxList{ind(1)}) = 1;
        vormask = newmask;

        % area-averaged vor/f
        dA = avg1(avg1( 1./S.pm .* 1./S.pn, 1), 2) .* vormask;

        nondim.eddy.Rovor = abs(sum(sum(vor./avg1(avg1(f,1),2) .* ...
                            vormask .* dA,1),2)) ./ nansum(dA(:));

        % max. azimuthal velocity for future
        [eddy.U,iU] = max(abs(rut(:,eddy.iy,end)));
        eddy.R = r(iU,eddy.iy);

        % calculate nondim parameters
        nondim.eddy.Ro = eddy.U ./mean(f(rnorm < r(iU,eddy.iy)))./eddy.R;
        % The check at this point makes no sense. If I have got here the
        % gradient wind quadratic has real root, so there is no issue. The
        % condition is that Ro using _geostrophic_ velocity is < 0.25
%        if  nondim.eddy.Ro > 0.25, error('Error: Ro > 0.25'); end
        nondim.eddy.Rh = eddy.U/phys.beta/eddy.R^2;
        nondim.eddy.Bu = (eddy.R / eddy.Ldef)^2;;
        nondim.eddy.Ri = N2./(TCOEF*g*eddy.tamp/f0/eddy.R).^2;
        nondim.eddy.Bu_temp = TCOEF *g * Z * eddy.tamp / f0^2 / eddy.R^2;
        nondim.eddy.gamma = bathy.hsb/eddy.depth;
        fprintf('\n Ro (vor/f) = %.2f | Bu = %.2f | Bu_temp = %.2f | Ri = %.2f | Rh = %.2f | Lsl/R = %.2f | H_sb/H_eddy = %.2f\n\n', ....
                nondim.eddy.Rovor,nondim.eddy.Bu,nondim.eddy.Bu_temp, ...
                nondim.eddy.Ri,nondim.eddy.Rh, bathy.L_slope/eddy.R, ...
                nondim.eddy.gamma);


        % calculate Ro using vorticity
        vor = avg1(diff(eddy.v(:,:,end),1,1)./diff(xrmat(:,:,end),1,1),2) - ...
              avg1(diff(eddy.u(:,:,end),1,2)./diff(yrmat(:,:,end),1,2),1);
        Ro1 = vor./avg1(avg1(f,1),2);
        Ro1 = max(abs(Ro1(:)));
        fprintf('\n Max. Ro (vor/f)  = %.2f \n', Ro1);

        Lr = eddy.Ldef;
        vgw1 = sqrt(N2) * Z;
        vr1 = -beta * Lr^2;
        % scaled height
        Hs = vgw1^2/g;
        vs1 = 14.5*vr1 * (vr1/vgw1) * nondim.eddy.Rh*86.4;
        vx = (bg.ubt-beta*Lr^2/2) * 86.4;
        tsl = ceil((eddy.cy-bathy.xsl)/1000/vs1);
        tsb = ceil((eddy.cy-bathy.xsb)/1000/vs1);
        fprintf(['\n Southward vel = %.3f km/day | ' ...
            'center reaches slope at t = %d days @ x=%4d km,' ...
            'shelfbreak at %d days, x=%4d km\n'], ...
            vs1,tsl,ceil(eddy.cx/1000+vx*tsl),tsb,ceil(eddy.cx/ ...
                                                       1000+vx*tsb));

        input('Continue?');
    end

    %% check plots
    xind = ymid; yind = ymid; zind = 20;

    if flags.front
        figure(h_check);
        axt(8) = subplot(248);
        cla
    else
        hfeddy = figure;
    end
    contourf(xrmat(:,:,1)./fx,yrmat(:,:,1)./fy,S.zeta,20); shading flat;
    hcb = colorbar; freezeColors; cbfreeze(hcb)
    hold on
    [C,h] = contour(xrmat(:,:,1)./fx,yrmat(:,:,1)./fy,S.h,...
                floor(linspace(min(S.h(:)),max(S.h(:)),5)),'k');
    clabel(C,h);
    liney((Y-eddy.buffer_sp)/fy,'sponge');
    linex((X-eddy.buffer_sp)/fx,'sponge');
    if flags.telescoping
        linex([S.x_rho(grid.ixn,1) S.x_rho(grid.ixp,1)]/fx,'telescope','w');
        liney([S.y_rho(1,grid.iyn) S.y_rho(1,grid.iyp)]/fy,'telescope', ...
              'w');
        contour(xrmat(:,:,1)./fx, yrmat(:,:,1)./fy, 1./S.pm, [1:0.5: ...
                            grid.dxmax/grid.dxmin] * grid.dxmin, 'w');
        contour(xrmat(:,:,1)./fx, yrmat(:,:,1)./fy, 1./S.pn, [1:0.5: ...
                            grid.dxmax/grid.dymin] * grid.dymin, 'w');
    end
    title('Zeta with eddy');
    if bathy.axis == 'y'
        liney([bathy.xsl bathy.xsb]/fy);
    else
        linex([bathy.xsl bathy.xsb]/fx);
    end
    axis image;
    if flags.front
        linkaxes([axt(4) axt(7) axt(8)],'xy');
    end
%     if flags.front
%         cbfreeze(hcb,'off');
%        spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))
%     end

    xedd = eddy.ix;
    yedd = eddy.iy;

    % check temperature profile
    figure;
    axe(1) = subplot(241);
    contourf(squeeze(xrmat(:,yedd,:))/fx,squeeze(zrmat(:,yedd,:)),squeeze(S.temp(:,yedd,:)),20);
    colorbar;
    title('temp (y=mid)');
    xlabel(['x' lx]); ylabel('z (m)');
    axe(2) = subplot(242);
    contourf(squeeze(yrmat(xedd,:,:))/fy,squeeze(zrmat(xedd,:,:)),squeeze(S.temp(xedd,:,:)),20);
    colorbar;
    title('temp (x=mid)');
    xlabel(['y' ly]); ylabel('z (m)');

    axe(3) = subplot(243);
    contourf(squeeze(xvmat(:,yedd,:))/fx,squeeze(zvmat(:,yedd,:)),squeeze(S.v(:,yedd,:)),20);
    colorbar;
    title('v (y=mid)');
    xlabel(['x' lx]); ylabel('z (m)');
    axe(4) = subplot(244);
    contourf(squeeze(yumat(xedd,:,:))/fy,squeeze(zumat(xedd,:,:)),squeeze(S.u(xedd,:,:)),20);
    colorbar;
    title('u (x=mid)');
    xlabel(['y' ly]); ylabel('z (m)');

    axe(5) = subplot(245);
    contourf(squeeze(xrmat(:,yedd,:))/fx,squeeze(zrmat(:,yedd,:)),squeeze(eddy.temp(:,yedd,:)),20);
    colorbar;
    title('Eddy temp (y=mid)');
    xlabel(['x' lx]); ylabel('z (m)');
    axe(6) = subplot(246);
    contourf(squeeze(yrmat(xedd,:,:))/fy,squeeze(zrmat(xedd,:,:)),squeeze(eddy.temp(xedd,:,:)),20);
    colorbar;
    title('Eddy temp (x=mid)');
    xlabel(['y' ly]); ylabel('z (m)');

    axe(7) = subplot(247);
    contourf(squeeze(xvmat(:,yedd,:))/fx,squeeze(zvmat(:,yedd,:)),squeeze(eddy.v(:,yedd,:)),20);
    colorbar;
    title('eddy v (y=mid)');
    xlabel(['x' lx]); ylabel('z (m)');
    axe(8) = subplot(248);
    contourf(squeeze(yumat(xedd,:,:))/fy,squeeze(zumat(xedd,:,:)),squeeze(eddy.u(xedd,:,:)),20);
    colorbar;
    title('eddy u (x=mid)');
    xlabel(['y' ly]); ylabel('z (m)');

    linkaxes([axe(1) axe(3) axe(5) axe(7)],'xy');
    linkaxes([axe(2) axe(4) axe(6) axe(8)],'xy');

    spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))

    %% estimate speed based on van leeuwin (2007)

    % first, K.E
    %rho_eddy = R0 * (- TCOEF * (eddy.temp));
    %ke = 1/2 * rho_eddy .* (eddy.u.^2 + eddy.v.^2);
    %pe = 1/2 * g * rho_eddy .* zrmat;

    %dV = bsxfun(@times, diff(zwmat, 1, 3), 1./S.pm .* 1./S.pn);

    %MASS = sum(rho_eddy(:) .* dV(:));
    %KE = sum(ke(:) .* dV(:));
    %PE = sum(pe(:) .* -dV(:));
    %gamma = (KE+PE)/MASS/g/Z;

    %h0 = Z/2; % average depth of upper layer
    %Ld = eddy.Ldef; % deformation radius
    %eta = eddy.temp(:,:,S.N/2);
    %num1 = 0;
    %num2 = f0^2/h0 * Ld^2 * trapz(yrmat(1,:,1),trapz(xrmat(:,1,1),eta.^2,1),2);
    %den = (2*Ld^2*trapz(yrmat(1,:,1),trapz(xrmat(:,1,1),eta,1),2));
    %gamma = (num1+num2)/den;

    Vr = beta * eddy.Ldef^2;

    %fprintf(['\n Based on van Leeuwin (2007) : \n ' ...
    %        'Gamma = %.2f, Vr = %.3f m/s, Drift speed = %.3f m/s \n'], ...
    %    gamma, Vr, Vr*(gamma));

    A = max(abs(eddy.zeta(:)));
    Nqg = f0*eddy.Ldef/g * Vr;
    gamma = Nqg/A - 1;

    % estimated westward translation speed
    % (3) in Early et al. (2011) is for quasi stable state expected at
    % t ~ 20/(beta * Ld) ~ 150-200 days. Based on fig. 8 in paper,
    % These are actual estimates. The reduction factor is applied later.
    Vest_zon = Vr * (Nqg/A - 1);
    Vest_mer = Vr/2 * Nqg/A * 1/2;
    fprintf(['\n Based on Early et al. (2011) : \n ' ...
            'Vr = %.3f m/s, Estimated (zonal, meridional) speed = (%.3f, %.3f) m/s\n'], ...
            Vr, Vest_zon, Vest_mer);
    if Vest_zon > 0, error('Estimated zonal velocity is EASTWARD'); ...
            end
    %%

    % clear variables to save space
%    clear rut rutz dTdr strat r rnorm rfb2 sdisc
    %eddy.tz = []; eddy.temp = []; eddy.u = []; eddy.v = [];

    fprintf('\n Eddy - %4.1f MB \n\n', monitor_memory_whos);

    toc;
end

%% add barotropic velocity for advection (OBC initial condition also)

% modify velocity and free surface fields
if flags.ubt_initial == 1
    if isnan(bg.ubt)
        bg.ubt = 1/2*abs(Vest_zon) + eddy.U/eddy.nl
    else
        if isnan(eddy.nl)
            eddy.nl = eddy.U/bg.ubt;
        end
    end
%     if flags.localize_jet
%         % find appropriate indices
%         idz1 = find_approx(squeeze(yrmat(1,:,1)),eddy.cy - eddy.dia/2 - buffer,1);
%         idz2 = find_approx(squeeze(yrmat(1,:,1)),eddy.cy + eddy.dia/2 + buffer,1);
%         idu1 = find_approx(squeeze(yumat(1,:,1)),eddy.cy - eddy.dia/2 - buffer,1);
%         idu2 = find_approx(squeeze(yumat(1,:,1)),eddy.cy + eddy.dia/2 + buffer,1);
%
%         S.u(:,idu1:idu2,:) = S.u(:,idu1:idu2,:) + ubt;
%
%         z1 = cumtrapz(squeeze(yrmat(1,:,1)),-f./g * ubt,2);
%         z1 = bsxfun(@minus,z1,z1(:,idz1));
%         S.zeta(:,idz1:idz2) = S.zeta(:,idz1:idz2) + z1(:,idz1:idz2,:);
%         S.zeta(:,idz2+1:end) = repmat(S.zeta(:,idz2),[1 size(S.zeta(:,idz2+1:end),2)]);
%         plot(S.zeta(30,:));
%
%     else
        if bg.ubt ~=0 && bg.vbt ~=0
            error('Adding both x and y barotropic velocity in initial condition. WILL NOT WORK!');
        end

        % modify ubt / vbt based on shear
        % bg.ubt , bg.vbt are scalars = value of background vel if no shear
        % if shear, then value of bg.ubt/vbt does not matter. one must be
        % zero to choose where to apply shear - u or v?
        if flags.bg_shear == 1
            if bg.vbt == 0 && bg.ubt == 0
                warning('No background velocity added');
            else
                bg.shear = bg.shear_fac * max(abs(vor(:)));
                if bg.vbt == 0
                    % need uy shear
                    bg.shear = bg.shear * -1;
                    uu  = bg.shear * (yrmat(1,:,1));
                    uu  = uu - mean(uu(:));
                    ubt = repmat(uu,[size(yrmat,1) 1]);
                else
                    if bg.ubt == 0
                        vv  = bg.shear * (xrmat(:,1,1));
                        vv  = vv - mean(vv(:));
                        vbt = repmat(vv,[1 size(xrmat,2)]);
                    end
                end
            end
        else
            ubt = bg.ubt * ones(size(f));
            vbt = bg.vbt * ones(size(f));
        end

        if bg.ubt ~= 0
            S.u = bsxfun(@plus,S.u,avg1(ubt,1));
            if flags.fplanezeta
                S.zeta = S.zeta + cumtrapz(squeeze(yrmat(1,:,1)),-f0./g * ubt,2);
            else
                S.zeta = S.zeta + cumtrapz(squeeze(yrmat(1,:,1)),-f./g .* ubt,2);
            end
        end
        if bg.vbt ~=0
            S.v = bsxfun(@plus,S.v,avg1(vbt,2));
            if ~flags.fplanezeta
               S.zeta = S.zeta + cumtrapz(squeeze(xrmat(:,1,1)),f./g .* vbt,1);
            else
               S.zeta = S.zeta + cumtrapz(squeeze(xrmat(:,1,1)),f0./g * vbt,1);
            end
        end
end

% write eddy params now (late) because eddy.nl might be calculated
% in this cell
fprintf('\n Writing eddy, strat, bg params');
write_params_to_ini(INI_NAME,eddy);
write_params_to_ini(INI_NAME,strat);
write_params_to_ini(INI_NAME,bg);

fprintf('\n BT vel - %4.1f MB \n\n', monitor_memory_whos);

%% non-dimensional parameters

if exist('S_sh','var');nondim.S_sh = S_sh; end
if exist('S_sl','var');nondim.S_sl = S_sl; end
nondim.comment = ['Ro = Rossby number | S_sh = shelf Burger number |' ...
                  'S_sl = slope Burger number | Bu_temp = Bu based on eddy temp perturbation' ...
                  ' | Bu = traditional eddy burger number = NH/fR | Ri = Richardson number | ' ...
                  ' Rh = Rhines number = V/beta/R^2'];

write_params_to_ini(INI_NAME,nondim);
fprintf('\n Parameters written \n');
toc;

%% Misc calculations (ubar,vbar,rho,pv) - shouldn't require changes

% Add random perturbation
zeta0 = S.zeta; % save for thermal wind check later
if flags.perturb_zeta == 1
    rng('shuffle');
    perturb = randn(size(S.zeta));
    perturb = perturb./max(abs(perturb(:))) * 10^(-4)*3;
    perturb = perturb - mean(perturb(:));
    S.zeta = S.zeta + perturb; % .* ~(xrmat < x1 | xrmat > x2);
end

% remove mean zeta
S.zeta = S.zeta - nanmean(S.zeta(:));

% ubar & vbar - use arango's uv_barotropic - works faster
if ~flags.spinup
    % recalculate to account for zeta stretching grid
    [z_w]=set_depth(S.Vtransform, S.Vstretching, ...
                     S.theta_s, S.theta_b, S.hc, S.N, ...
                     5, S.h, S.zeta,0);
    [S.ubar,S.vbar] = uv_barotropic(S.u,S.v,Hz);
end

S.rho = R0*(1 - TCOEF * (S.temp-T0) + SCOEF * (S.salt-S0));

toc;

% setup for pv calculation
if calc_pv
    grid1.xu = xumat;
    grid1.yu = yumat;
    grid1.zu = zumat;

    grid1.xv = xvmat;
    grid1.yv = yvmat;
    grid1.zv = zvmat;

    grid1.xr = xrmat;
    grid1.yr = yrmat;
    grid1.zr = zrmat;
    grid1.zw = zwmat;

    grid1.s_w = S.s_w; grid1.s_rho = S.s_rho;

    [pv,xpv,ypv,zpv] = pv_cgrid(grid1,S.u,S.v,S.rho,f,R0);

    pvmin = min(pv(:));
    pvmid = pv(xmid,ymid,zmid);

    % figure;
    % contourf(xpv,zpv,squeeze(pv(:,ymid,:))',40);
    % colorbar;
    % title(['PV | PV_{min}/PV_{mid} = ' num2str(pvmin/pvmid)]);
    % xlabel('x'); ylabel('z (m)');
    clear grid1
end

fprintf('\n ubar,vbar,pv- %4.1f MB \n\n', monitor_memory_whos);

%% Sanity Checks + check thermal wind

check_thermalwind = 0;


VZ  = diff(S.v,1,3) ./ diff(zvmat,1,3);
UZ  = diff(S.u,1,3) ./ diff(zumat,1,3);

if check_thermalwind
    % calculate % error in thermal wind balance
    xind = 30; yind = 15; zind = 30;

    % these are needed later too - calculated for geostrophic component of
    % % velocity
    % VZ  = diff(     avg1(bsxfun(@times,vgeo, cos(th)),2),1,3) ./ diff(zvmat,1,3);
    % UZ  = diff(-1 * avg1(bsxfun(@times,vgeo, sin(th)),1),1,3) ./ diff(zumat,1,3);

    dxarr = squeeze(diff(xrmat(:,:,end),1,1));
    dyarr = squeeze(diff(yrmat(:,:,end),1,2));

    dRdx = bsxfun(@rdivide, diff_cgrid(tgrid,S.rho,1) * g/R0, avg1(f,1));
    dRdy = bsxfun(@rdivide, diff_cgrid(tgrid,S.rho,2) * g/R0, avg1(f,2));

    dzdx = bsxfun(@rdivide, diff(zeta0,1,1), dxarr);
    dzdy = bsxfun(@rdivide, diff(zeta0,1,2), dyarr);

    % PERCENTAGE error in thermal wind balance
    diff_u = (avg1(UZ,2) - avg1(avg1(dRdy,1),3))./max(abs(UZ(:))) * 100;
    diff_v = (avg1(VZ,1) + avg1(avg1(dRdx,2),3))./max(abs(VZ(:))) * 100;
    %
    % figure;
    % subplot(211)
    % pcolorcen(squeeze(diff_v(:,yind,:))'); colorbar
    % title('% error (v_z - dT/dx * \alpha g/f)');
    % subplot(212)
    % pcolorcen(squeeze(diff_u(xind,:,:))'); colorbar
    % title('% error (-u_z - dT/dy * \alpha g/f)/max(u_z)');

    figure;
    subplot(221)
    pcolorcen(diff_v(:,:,zmid)'); colorbar
    title('% error (v_z - d\rho/dx * g/f/R0)|_{z=mid}');
    subplot(222)
    pcolorcen(squeeze(diff_u(:,:,zmid))'); colorbar
    title('% error (-u_z - d\rho/dy * g/f/R0)|_{z=mid}');
    subplot(223)
    pcolorcen((((avg1((dzdx ./ avg1(f,1)* g),2)-avg1(S.v(:,:,end),1)))./ ...
                    max(abs(S.v(:)))*100)'); colorbar
    title('% error (g/f d\zeta /dx - v_{z=0})');
    subplot(224)
    pcolorcen(((avg1(dzdy ./ avg1(f,2)* g,1) + avg1(S.u(:,:,end),2)) ./ ...
                max(abs(S.u(:)))*100)'); colorbar
    title(' % error (g/f d\zeta /dy + u_{z=0})');
end

% sanity checks
if min(S.temp(:)) < 3
    error('Temperature less than or close to 0.');
end
if max(S.zeta(:)) > 1
    error('Zeta > 1m.');
end

if min(S.salt(:)) == 0 && (S.NT-S.NPT) > 1
    error('Salt set to zero');
end

if max(abs(S.v(:))) > 1.0 || max(abs(S.u(:))) > 1.0
    fprintf('\n max(u) = %.3f m/s | max(v) = %.3f m/s \n',max((abs(S.v(:)))), max((abs(S.u(:)))));
    input('Really high velocities. Are you sure?');
end

Ri = abs(fillnan(N2./(avg1(VZ.^2,1) + avg1(UZ.^2,2)),Inf));
if min(Ri(:)) <= 0.3
    figure; imagesc(addnan(Ri(:,:,end)',10)); title('Ri < 10'); colorbar;
    error('Ri <= 0.3.');
else
    fprintf('\n Min. Ri = %.3f\n\n', nanmin(Ri(:)));
end

Ri = nanmin(Ri(:));

clear dRdx dRdy dzdx dzdy

fprintf('\n Sanity checks - %4.1f MB \n\n', monitor_memory_whos);

%% Passive Tracers - Initial condition + writing

if S.NPT > 0
    % create variables first
    names = {'cross shelf dye'; 'z dye'; 'eddy dye'; 'along shelf dye'};
    try
        dc_roms_passive_tracer(S,names,1);
    catch ME
        warning('Passive tracer variables already created?');
    end

    % assign initial condition
    % important, we do not want grid scale discontinuity.
    % taper everything with gaussian or similar

    % set dye_01 = cross-shelf label & dye_03 = eddy dye
    if bathy.axis == 'y'
        if ~flags.front
            % cross-shelf dye
            dye_01 = yrmat;
        end
        % along-shelf dye
        dye_04 = xrmat;
    else
        if ~flags.front
            % cross shelf dye
            dye_01 = xrmat;
        end
        % along shelf dye
        dye_04 = yrmat;
    end
    % set cross-shelf dye (dye_01) same as frontal structure
    if flags.front
        dye_01 = S.Tra;
    end

    % second dye is always z-label
    dye_02 = zrmat;

    % fourth dye tags eddy
     % set eddy dye
     if flags.eddy
        dye_03 = zeros(size(yrmat));
        %1e-3 is good threshold for tamp=0.4
        dye_03(eddy.temp > (1e-3/0.4*eddy.tamp)) = 1;
        % now smooth out edges
        nfilter=5;
        for kk=1:S.N
            dye_03(:,:,kk) = filter2(ones(nfilter,nfilter)/nfilter.^2, ...
                                dye_03(:,:,kk));
        end
      end
%     dye_01 = zeros(size(xrmat)); dye_02 = dye_01;
%
%     buffer = 6;
%     % decrease from 1 to 0
%     buffer_val_dec = repmat(cos( pi/2*(0:buffer-1)/(buffer-1) ), ...
%                                 [size(dye_01,1) 1 size(dye_01,3)]);
%     % increase from 0 to 1
%     buffer_val_inc = repmat(fliplr(cos( pi/2*(0:buffer-1)/(buffer-1) )), ...
%                                 [size(dye_01,1) 1 size(dye_01,3)]);
%
%     dye_01(:,1:bathy.isb,:) = 1;
%     dye_01(:,bathy.isb+1:bathy.isb+buffer,:) = buffer_val_dec;
%
%     dye_02(:,bathy.isb-buffer:bathy.isb-1,:) = buffer_val_inc;
%     dye_02(:,bathy.isb:bathy.isl,:) = 1;
%     dye_02(:,bathy.isl+1:bathy.isl+buffer,:) = buffer_val_dec;
%
%     figure;
%     subplot(211)
%     pcolorcen(dye_01(:,:,end)');
%     title('dye_{01}');
%     subplot(212)
%     pcolorcen(dye_02(:,:,end)');
%     title('dye_{02}');

    % write to file
    for ii=1:S.NPT
        vname = sprintf('dye_%02d',ii);
        eval(['ncwrite(S.ncname,''' vname ''',' vname ');']);
    end

    fprintf('\n Passive Tracer - %4.1f MB \n\n', monitor_memory_whos);
end

%% Floats - figure out seeding locations

if flags.floats
    if ~flags.front
       figure(hfeddy);
    end
    str = 'select rectangle for float deployment';
    title(str);
    disp(str);
    %[floatx,floaty] = select_rect();
    fprintf('\n Floats- %4.1f MB \n\n', monitor_memory_whos);
end

%% Open Boundary Conditions - Initialize and write

% modified from arango/d_obc_roms2roms.m
if flags.OBC == 1
    % copy from initial conditions structure
    Sbr = S;

    Sbr.ncname = BRY_NAME;
    %  Set switches for boundary segments to process
    Sbr.boundary(1) = OBC.west;
    Sbr.boundary(2) = OBC.east;
    Sbr.boundary(3) = OBC.south;
    Sbr.boundary(4) = OBC.north;

    VarGrd = {'spherical',                                                ...
          'Vtransform', 'Vstretching',                                ...
          'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
          's_rho', 'Cs_r', 's_w', 'Cs_w'};

    % create junk file with proper attributes
    [status] = c_boundary(Sbr);

    % write OBC data, grid and passive tracer info to OBC file
    dc_roms_create_bry_file(Sbr);

    bry_time = 0;

    boundaries = fieldnames(OBC);
    range   = {'(1,:,:)','(end,:,:)','(:,1,:)','(:,end,:)'};
    VarList = {'zeta'; 'u'; 'ubar';'v'; 'vbar';'temp';'salt'}; % variables to write

    % set boundary condition data
    if flags.OBC_from_initial
        % if background flow is in initial condition, copy that to boundary
        % condition

        for mm = 1:size(boundaries,1);
            if OBC.(char(boundaries(mm,:))) % if open boundary
                for jj = 1:size(VarList,1)
                    varname = sprintf('%s_%s',char(VarList{jj}),char(boundaries(mm,:)));
                    eval([varname '= squeeze(S.' char(VarList{jj}) char(range{mm}) ');']);
                end
            end
        end

    else % analytically calculate boundary condition
        % Set boundary condition fields
        % if flags.ubt_deep
        %     idy1_x = find(S.h(end,:) == max(S.h(end,:)));
        %     idy2_x = idy1_x(end);
        %     idy1_x = idy1_x(1);
        %
        %     % ADD for vbar if needed
        % else
        %    idy1_x = 1;
        %    idy2_x = S.Mm+2;
        %
        %    idy1_y = 1;
        %    idy2_y = S.Mm+1;
        % end
        if OBC.east
                ubar_east = zeros(size(S.u(end,:,1)));
                ubar_east(1,idy1_x:idy2_x) = ubt;

                vbar_east = zeros(size(S.v(end,:,1)));

                zeta_east = zeros(size(S.zeta(end,:)));
                zeta_east(end,:) = zeta_east(1,:) + cumtrapz(squeeze(yrmat(1,:)), ...
                                                    -f(end,:)./g * ubt,2);
                zeta_east = zeta_east - nanmean(zeta_east);
        end

        if OBC.north
                ubar_north = zeros(size(S.u(:,end,1)));

                vbar_north = zeros(size(S.v(:,end,1)));
                vbar_north(1,:) = vbt;

                zeta_north = zeros(size(S.zeta(:,end)));
                zeta_north(:,end) = zeta_north(:,1) + cumtrapz(squeeze(xrmat(:,end,1)), ...
                                                    f(:,end)./g * vbt,1);
                zeta_north = zeta_north - nanmean(zeta_north);
        end
    end

    % Write boundary conditions to file
    for mm = 1:size(boundaries,1);
        if OBC.(char(boundaries(mm,:))) % if open boundary
            for jj = 1:size(VarList,1)
                varname = sprintf('%s_%s',char(VarList{jj}),char(boundaries(mm,:)));
                disp(['Writing ' varname]);
                eval(['nc_write(Sbr.ncname,''' varname ''',' varname ',1);']);
            end
        end
    end
    nc_write(Sbr.ncname,'bry_time',bry_time,1);

    % set and write passive tracer data boundary conditions
    if S.NPT > 0
       for mm=1:size(boundaries,1)
           if OBC.(char(boundaries(mm,:)))
               for ii=1:S.NPT
                   varname = sprintf('dye_%s_%02d',char(boundaries(mm,:)),ii);
                   disp(['Writing ', varname]);
                   eval([varname ' = ' sprintf('dye_%02d',ii) char(range{mm}) ';']);
                   eval(['nc_write(Sbr.ncname,''' varname ''',' varname ',1);']);
               end
           end
       end
    end
    fprintf('\n OBC - %4.1f MB \n\n', monitor_memory_whos);
end

%% Wind Forcing - Initialize and write

if flags.wind
    % first find width of deep water region
    indh = find(S.h(:,1) == max(S.h(:)),1);
    iX = size(xrmat,1);
    wind.Lx = xrmat(end,1,end) - xrmat(indh,1,end);

    % beta * V (depth integrated) = curl(tau)/rho0
    wind.tau0 = rho0 * phys.beta * wind.Lx * wind.v * Z;

    sms_time = [0 wind.ramp 300] * 86400;

    sustr = zeros([size(S.h) length(sms_time)]);
    svstr = zeros([size(S.h) length(sms_time)]);

    svstr(indh:end,:,2) = (0 + (wind.tau0 - 0) *  (xrmat(indh:end,:,1)-xrmat(indh,1,1)) ...
                        /wind.Lx)/rho0;
    svstr(:,:,3) = svstr(:,:,2);

    % create file
    dc_roms_wind_forcing(S,FRC_NAME);

    % write to file
    ncwrite(FRC_NAME,'sms_time',sms_time);
    ncwrite(FRC_NAME,'sustr',sustr);
    ncwrite(FRC_NAME,'svstr',svstr);

    fprintf('\n Wind - %4.1f MB \n\n', monitor_memory_whos);

    write_params_to_ini(INI_NAME,wind);
    toc;
end

%% Check plots

make_plot = 0;

if make_plot

    limx = [min(xrmat(:,1,1)) max(xrmat(:,1,1))]./fx;
    limy = [min(yrmat(1,:,1)) max(yrmat(1,:,1))]./fy;

    if flags.eddy
        xind = xedd;
        yind = yedd;
    else
        xind = xmid;
        yind = ymid;
    end

    %% Plot all fields
    figure;
    ax(1) = subplot(241);
    contourf(squeeze(yumat(xind,:,:))./fx,squeeze(zumat(xind,:,:)),squeeze(S.u(xind,:,:)),20);
    colorbar;
    title('u');
    xlabel(['y ' lx]); ylabel('z (m)');

    ax(2) = subplot(242);
    contourf(squeeze(yvmat(xind,:,:))./fx,squeeze(zvmat(xind,:,:)),squeeze(S.v(xind,:,:)),20);
    colorbar;
    title('v');
    xlabel(['y ' lx]); ylabel('z (m)');

    ax(3) = subplot(243);
    contourf(squeeze(yrmat(xind,:,:))./fy,squeeze(zrmat(xind,:,:)),squeeze(S.rho(xind,:,:)),20);
    colorbar;
    title('rho');
    xlabel(['y ' lx]); ylabel('z (m)');

    % fix pv script so that it returns proper co-ordinates
    if calc_pv
        ax(4) = subplot(244);
        contourf(repmat(xpv(:,1),1,size(pv,3))./fx,squeeze(zpv(:,ymid,:)),squeeze(pv(:,ymid,:)),20);
        colorbar;
        title('PV');
        xlabel(['x ' lx]); ylabel('z (m)');
    end

    % need to interpolate to constant z - surface
    ax(5) = subplot(245);
    contourf(xumat(:,:,end)./fx,yumat(:,:,end)./fy,squeeze(S.u(:,:,end)));
    colorbar;
    title('surface u');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    linex((eddy.cx-eddy.dia/2)/fx);
    linex((eddy.cx+eddy.dia/2)/fx);
    liney((eddy.cy-eddy.dia/2)/fy);
    liney((eddy.cy+eddy.dia/2)/fy);
    axis square;

    ax(6) = subplot(246);
    contourf(xvmat(:,:,end)./fx,yvmat(:,:,end)./fy,squeeze(S.v(:,:,end)));
    colorbar;
    title('surface v');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    linex((eddy.cx-eddy.dia/2)/fx);
    linex((eddy.cx+eddy.dia/2)/fx);
    liney((eddy.cy-eddy.dia/2)/fy);
    liney((eddy.cy+eddy.dia/2)/fy);
    axis square;

    ax(7) = subplot(247);
    contourf(xrmat(:,:,end)./fx,yrmat(:,:,end)./fy,squeeze(S.rho(:,:,end)));
    colorbar;
    linex((eddy.cx-eddy.dia/2)/fx);
    linex((eddy.cx+eddy.dia/2)/fx);
    liney((eddy.cy-eddy.dia/2)/fy);
    liney((eddy.cy+eddy.dia/2)/fy);
    title('SSRho');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    ax(8) = subplot(248);
    contourf(S.x_rho(:,1)./fx,squeeze(S.y_rho(1,:))./fy,squeeze(S.zeta(:,:,1))');
    colorbar;
    linex((eddy.cx-eddy.dia/2)/fx);
    linex((eddy.cx+eddy.dia/2)/fx);
    liney((eddy.cy-eddy.dia/2)/fy);
    liney((eddy.cy+eddy.dia/2)/fy);
    title('SSH (zeta)');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    linkaxes([ax(1:4)],'xy');
    linkaxes([ax(5:8)],'xy');

    spaceplots(gcf,0.03*ones([1 4]),0.05*ones([1 4]))

end

%% clear some vars

clear pv

%% Write to Grid & IC file

toc;
fprintf('\n Started writing files...\n');

% write git hash!
[~, hash] = system('TERM=xterm-256color git log -n 1 --pretty=format:''%H''');
% remove bash escape characters
hash = hash(9:48)

ncwriteatt(GRID_NAME, '/', 'git_hash', hash);
ncwriteatt(INI_NAME, '/',  'git_hash', hash);
ncwriteatt(BRY_NAME, '/',  'git_hash', hash);

% grid file
ncwrite(GRID_NAME,'xl',S.xl);
ncwrite(GRID_NAME,'el',S.el);
ncwrite(GRID_NAME,'f',f);
ncwrite(GRID_NAME,'h',S.h);
ncwrite(GRID_NAME, 'mask_u',    S.mask_u);
ncwrite(GRID_NAME, 'mask_v',    S.mask_v);
ncwrite(GRID_NAME, 'mask_rho',  S.mask_rho);
ncwrite(GRID_NAME, 'mask_psi',  S.mask_psi);
ncwrite(GRID_NAME, 'spherical', S.spherical);
ncwrite(GRID_NAME, 'pm',        S.pm);
ncwrite(GRID_NAME, 'pn',        S.pn);
ncwrite(GRID_NAME, 'dndx',      S.dndx);
ncwrite(GRID_NAME, 'dmde',      S.dmde);
ncwrite(GRID_NAME, 'angle',     S.angle);

ncwrite(INIname,   'spherical',   S.spherical);
ncwrite(INIname,   'Vtransform',  S.Vtransform);
ncwrite(INIname,   'Vstretching', S.Vstretching);
ncwrite(INIname,   'theta_s',     S.theta_s);
ncwrite(INIname,   'theta_b',     S.theta_b);
ncwrite(INIname,   'Tcline',      S.Tcline);
ncwrite(INIname,   'hc',          S.hc);
ncwrite(INIname,   's_rho',       S.s_rho);
ncwrite(INIname,   's_w',         S.s_w);
ncwrite(INIname,   'Cs_r',        S.Cs_r);
ncwrite(INIname,   'Cs_w',        S.Cs_w);

% IC file

IniRec = 1;                               % NetCDF time record

S.ocean_time = 0;                % initial conditions time (s)

ncwrite(INIname, 'ocean_time', S.ocean_time);
ncwrite(INIname,   'h',           S.h);

if (S.spherical),
    ncwrite(INIname, 'lon_rho',     S.lon_rho);
    ncwrite(INIname, 'lat_rho',     S.lat_rho);
    ncwrite(INIname, 'lon_u',       S.lon_u);
    ncwrite(INIname, 'lat_u',       S.lat_u);
    ncwrite(INIname, 'lon_v',       S.lon_v);
    ncwrite(INIname, 'lat_v',       S.lat_v);

    ncwrite(GRID_NAME, 'lon_rho',     S.lon_rho);
    ncwrite(GRID_NAME, 'lat_rho',     S.lat_rho);
    ncwrite(GRID_NAME, 'lon_u',       S.lon_u);
    ncwrite(GRID_NAME, 'lat_u',       S.lat_u);
    ncwrite(GRID_NAME, 'lon_v',       S.lon_v);
    ncwrite(GRID_NAME, 'lat_v',       S.lat_v);
else
    ncwrite(INIname, 'x_rho',       S.x_rho);
    ncwrite(INIname, 'y_rho',       S.y_rho);
    ncwrite(INIname, 'x_u',         S.x_u);
    ncwrite(INIname, 'y_u',         S.y_u);
    ncwrite(INIname, 'x_v',         S.x_v);
    ncwrite(INIname, 'y_v',         S.y_v);

    ncwrite(GRID_NAME, 'x_rho',       S.x_rho);
    ncwrite(GRID_NAME, 'y_rho',       S.y_rho);
    ncwrite(GRID_NAME, 'x_u',         S.x_u);
    ncwrite(GRID_NAME, 'y_u',         S.y_u);
    ncwrite(GRID_NAME, 'x_v',         S.x_v);
    ncwrite(GRID_NAME, 'y_v',         S.y_v);
    ncwrite(GRID_NAME, 'x_psi',       S.x_psi);
    ncwrite(GRID_NAME, 'y_psi',       S.y_psi);
end,
fprintf('\n Grid stuff written \n');
toc;

ncwrite(INIname, 'zeta', S.zeta);
ncwrite(INIname, 'ubar', S.ubar);
ncwrite(INIname, 'vbar', S.vbar);
ncwrite(INIname, 'u',    S.u);
ncwrite(INIname, 'v',    S.v);
ncwrite(INIname, 'temp', S.temp);
ncwrite(INIname, 'salt', S.salt);
fprintf('\n Hydrodynamic vars written \n');
toc;

fprintf('\n\n Files %s | %s ', GRID_NAME,INI_NAME);
if flags.OBC, fprintf('| %s ',BRY_NAME); end
if flags.wind, fprintf('| %s ',FRC_NAME); end
fprintf('written.\n');
toc;

%% Grid and time step information

dx = min(1./S.pm(:)); dy = min(1./S.pn(:));
DX = sqrt(min(dx)^2 + min(dy)^2);
fprintf('\n\n=============================  SUMMARY  ========================================');
fprintf(['\n\n (nx,ny,nz) = (%d,%d,%d) | (X,Y,Z) = (%.2f , %.2f , %.2f) m |' ...
        '(dx,dy) = ( %.2f , %.2f) m'], ...
        S.Lm,S.Mm,S.N,max(xrmat(:)), max(yrmat(:)), max(abs(zrmat(:))),dx, dy);
fprintf('\n\n Beckmann & Haidvogel number = %f (< 0.2 , max 0.4) \n \t\t\t\tHaney number = %f (< 9 , maybe 16)', rx0,rx1);

% From Utility/metrics.F
% Barotropic courant number
dt = 60;
ndtfast = 17;
Cbt = sqrt(g*max(S.h(:))) * dt/ndtfast * sqrt(1/dx^2 + 1/dy^2);
Cbc = sqrt(N2)*min(S.h(:))/pi * dt * sqrt(1/dx^2 + 1/dy^2);
Cbc7 = 7 * dt * sqrt(1/dx^2 + 1/dy^2);
fprintf('\n\n (dt)_bt < %.2f s | (dt)_bc < %.2f s\n\n', DX/(sqrt(g*min(S.h(:)))), DX/(sqrt(N2)*min(S.h(:))/pi))
fprintf('\n\n Assuming dt = %.2f, ndtfast = %d, C_bt = %.3f | C_bc = %.3f  | C_bc7 = %.3f\n\n Min. Ri = %.2f', dt,ndtfast,Cbt,Cbc,Cbc7,min(Ri(:)));

fprintf('\n Bathy Parameters');
if ~flags.flat_bottom
    fprintf('\n\t\t S_shelf = %.2f | S_slope = %.2f | Beta_shelf = %1.2e | Beta_slope = %1.2e \n', ...
    S_sh,S_sl,b_sh,b_sl);
end
if flags.wind, cprintf('Red',sprintf('Wind tau0 = %.2e \n\n',wind.tau0)); end
if flags.eddy,
    fprintf('\n Eddy Parameters: ');
    fprintf('\n Ro (vor/f) = %.2f | max. Ro (vor/f)= %.2f | Bu = %.2f | Bu_temp = %.2f | Ri = %.2f | Rh = %.2f | Lsl/R = %.2f | H_sb/H_eddy = %.2f\n\n', ....
            nondim.eddy.Rovor, Ro1,nondim.eddy.Bu,nondim.eddy.Bu_temp, ...
            nondim.eddy.Ri,nondim.eddy.Rh, bathy.L_slope/eddy.R, ...
            nondim.eddy.gamma);

    fprintf('\n Deploy float in center of eddy = (%d,%d) \n\n',eddy.ix,eddy.iy);
end
if flags.floats
   xlo = find_approx(xrmat(:,1,1),min(floatx(:))*1000,1);
   xhi = find_approx(xrmat(:,1,1),max(floatx(:))*1000,1);
   ylo = find_approx(yrmat(1,:,1),min(floaty(:))*1000,1);
   yhi = find_approx(yrmat(1,:,1),max(floaty(:))*1000,1);
   fprintf('\n Float deployment locations : (%d:%d , %d:%d)', xlo,xhi,ylo,yhi);
end
if beta == 0, warning('f-plane!!!'); end
fprintf('\n\n');

% figure out appropriate viscosity
% I need min sqrt(visc4) = 42 :) for stability with dt=300 and good
% results - see runew-05
grdscl = sqrt(1./S.pm .* 1./S.pn);
grdmax = max(grdscl(:));
factor = (grdscl/grdmax).^3;
fprintf('Max visc4 = %e  | max diff4 = %e\n', (sqrt(5e6)/min(factor(:)))^2, ...
        (18/min(factor(:)))^2);

%% Old sbfront code

%     [S.Trax,mask_z] = hor_grad_tracer(axmat,ax_as,ax_cs,zrmat,i_cs,i_as,front,bathy);
% %     S.Tz = mask_z .* exp(-(zrmat./front.LTz).^2);
% %     S.Tz = S.Tz ./ nanmax(S.Tz(:));
% %     S.Tz = N2/g/TCOEF .* (repnan(S.Tz,1));
%     S.Traz = 0;%N2/g/coef .* ones(size(zrmat));
%
%     % use chain rule to get gradient on SIGMA LEVEL
%     dzdx_s = diff(zrmat,1,i_cs)./diff(axmat,1,i_cs);
%     S.Trax_sig = avg1(S.Trax,i_cs) + dzdx_s .* avg1(S.Traz,i_cs);
%
%     % then integrate to make T front
%     % top to bottom integration done earlier
%     % integrate in horizontal dirn.from left to right using gradients on SIGMA LEVEL
%     [S,axmat] = reset_flip(S,axmat);
%     [S,axmat] = flip_vars(flip_flag,S,axmat);
%     if bathy.loc == 'h'
%         for i=size(S.Tra,1)-1:-1:1
%             S.Tra(i,:,:) = S.Tra(i+1,:,:) + S.Trax_sig(i,:,:)  .*(axmat(i+1,:,:)-axmat(i,:,:));
%         end
%     else
%         for i=2:size(S.Tra,1)
%             S.Tra(i,:,:) = S.Tra(i-1,:,:) + S.Trax_sig(i-1,:,:).*(axmat(i,:,:)-axmat(i-1,:,:));
%         end
%     end
%     S = flip_vars(flip_flag,S);

%% old bathy code

%         if flags.old_bathy
%             bathy.sl_slope2 = bathy.sl_slope;
%
%             % main shelf
%             [hx,hy,hdeep] = bathy_crooked(S.x_rho,S.y_rho,bathy,X,Y);
%
%             % second section
%             B2 = bathy;
%             B2.L_entry = 0 *1000;
%             B2.L_shelf = 50*1000;
%             [hx2,hy2,~] = bathy_crooked(S.x_rho,S.y_rho,B2,X,Y);
%
%             h1 = hx.*hy .*~(S.x_rho > (X-bathy.L_entry-bathy.L_slope));
%             h2 = hx2 .* (S.y_rho > Y-B2.L_shelf) .*(S.x_rho > (X-bathy.L_entry-bathy.L_slope));
%
%         %     subplot(221); imagescnan(hx'); set(gca,'ydir','normal')
%         %     subplot(222); imagescnan(hy'); set(gca,'ydir','normal')
%         %     subplot(223); imagescnan(hx2');  set(gca,'ydir','normal')
%         %     subplot(224); imagescnan(hy2');  set(gca,'ydir','normal')
%
%             % final bathymetry
%             S.h = h1+h2;
%             S.h(S.h == 0) = hdeep;
%         end

%
%     if flags.tanh_bathymetry == 1
%         scale = 20000;
%         bathy.L_deep = Y - bathy.L_shelf;
%         bathy.H_deep  = Z;
%         S.hflat = Z*ones(size(S.h)); % constant depth
%         % y-z profile
%         %hy = H_deep + (H_shelf-H_deep)*(1+tanh( ((S.y_rho-L_deep)/20000) ))/2;
%         %hx = H_shelf - (H_shelf-H_deep)*(1+tanh( ((S.x_rho-L_entry)/20000) ))/2;
%
%         hx = (1-tanh( ((S.x_rho-L_entry)/scale) ))/2;
%         hy = (1+tanh( ((S.y_rho-L_deep)/scale) ))/2;
%         S.h = H_deep + (H_shelf-H_deep) * (hx.*hy);
%     end

%% using Hernan's balance operator code (with prsgd31 algo)

% K.Lm = S.Lm;
% K.Mm = S.Mm;
% K.N = S.N;
% K.LNM_depth = 0; % integrate from bottom for FS
% K.g = 9.81;
% K.rho0 = 1025;
% K.alpha = TCOEF * ones(size(S.temp,1),size(S.temp,2));
% K.beta  = SCOEF * ones(size(S.temp,1),size(S.temp,2));
%
% K.f = f;
% K.pm = S.pm;
% K.pn = S.pn;
% K.pmon_u = avg1(S.pm ./ S.pn , 1);
% K.pnom_v = avg1(S.pn ./ S.pm , 2);
% K.rmask = S.mask_rho;
% K.umask = S.mask_u;
% K.vmask = S.mask_v;
% K.Hz = Hz;
% K.Zr = zrmat;
% K.Zw = zwmat;
%
% drho = rho_balance(K,S.temp-T0,S.salt-S0);
% [u,v,zeta_rhs] = uv_balance(K,drho);
% [zeta,~] = zeta_balance(K,0,drho);

%% old calculate velocity field
% first re-calculate temperature (density) gradient
% tgrid.xmat = xrmat; tgrid.ymat = yrmat; tgrid.zmat = zrmat; tgrid.s = S.s_rho;
% rgrid.zw = permute(zwmat,[3 2 1]); rgrid.s_w = S.s_w;
%
% Tgrad1 = avg1(avg1(horgrad_cgrid(rgrid,tgrid,S.temp,i_cs),i_cs),i_as);
% Tgrad = nan(size(velzmat));
% if flip_flag
%     Tgrad1 = permute(Tgrad1, [2 1 3]);
%     Tgrad = permute(Tgrad, [2 1 3]);
% end
% % % add values at end
% Tgrad1(size(Tgrad1,1)+1,:,:) = Tgrad1(size(Tgrad1,1),:,:);
% Tgrad(1,:,:) = Tgrad1(2,:,:);
% Tgrad(2:end,:,:) = Tgrad1;
% clear Tgrad1
% if flip_flag, Tgrad = permute(Tgrad, [2 1 3]); end

%% check gradient calculation
%
% tgrid.xmat = xrmat; tgrid.ymat = yrmat; tgrid.zmat = zrmat; tgrid.s = S.s_rho;
% rgrid.zw = permute(zwmat,[3 2 1]); rgrid.s_w = S.s_w;

% if flags.front
%     Tgrad1 = avg1(avg1(diff_cgrid(tgrid,S.temp,i_cs),i_cs),i_as);
%     Tgrad = nan(size(velzmat));
%     if flip_flag
%         Tgrad1 = permute(Tgrad1, [2 1 3]);
%         Tgrad = permute(Tgrad, [2 1 3]);
%     end
%     % % add values at end
%     Tgrad1(size(Tgrad1,1)+1,:,:) = Tgrad1(size(Tgrad1,1),:,:);
%     Tgrad(1,:,:) = Tgrad1(2,:,:);
%     Tgrad(2:end,:,:) = Tgrad1;
%     clear Tgrad1
%     if flip_flag, Tgrad = permute(Tgrad, [2 1 3]); end
%
%
%     dT = abs(Tgrad-avg1(S.Trax,i_as))./max(abs(S.Trax(:))) * 100;
%     yind = S.Mm/2;
%     figure
%     ax(1) = subplot(131);
%     contourf(squeeze(xrmat(:,yind,:)), squeeze(zrmat(:,yind,:)), squeeze(S.Trax(:,yind,:)));
%     title('imposed T_x');
%     colorbar; clim = caxis;
%     ax(2) = subplot(132);
%     contourf(squeeze(xrmat(:,yind,:)), squeeze(zrmat(:,yind,:)), squeeze(Tgrad(:,yind,:)));
%     title('Calculated T_x');
%     caxis(clim); colorbar
%     ax(3) = subplot(133);
%     contourf(squeeze(xrmat(:,yind,:)), squeeze(zrmat(:,yind,:)), squeeze(dT(:,yind,:)));
%     title('% error');
%     colorbar
%     linkaxes(ax,'xy');
% end

%% check eddy profile (normalized)
% subplot(131)
% pcolorcen(S.x_rho/1000,S.y_rho/1000,eddy.xyprof);
% axis('square');colorbar; xlabel('x (km)'); ylabel('y (km)');
% subplot(132)
% plot(S.x_rho(:,1)/1000,eddy.xyprof(:,ymid),'*-');xlabel('x (km)');
% subplot(133)
% plot(eddy.zprof,squeeze(zrmat(1,1,:)),'*-'); ylabel('z (m)');

%% cosine eddy profile
    % temp. field -  z-profile (normalized)
%    % cos profile
%     ind = find_approx(deep_z, -1 * eddy.depth,1);
%     eddy.zprof = [zeros(ind-eddy.Ncos,1); (1-cos(pi * [0:eddy.Ncos]'/eddy.Ncos))/2; ones(S.N-ind-1,1)];

%% old ubar,vbar integration
%     for ii=1:S.Lm+2
%         for jj=1:S.Mm+2
%             if ii~=S.Lm+2
%                 ubar(ii,jj) = trapz(squeeze(zumat(ii,jj,:)),S.u(ii,jj,:),3)./hrmat(ii,jj);
%             end
%             if jj~=S.Mm+2
%                 vbar(ii,jj) = trapz(squeeze(zvmat(ii,jj,:)),S.v(ii,jj,:),3)./hrmat(ii,jj);
%             end
%         end
%     end

%% zeta from dynamic height - does not work
%     S.zeta = nan(size(zrmat(:,:,1)));
%     dha = nan(size(zrmat));
%     for mm = 1:size(zrmat,1)
%         % sw_gpan requires matrices as (z,x)
%         s0 = squeeze(S.salt(mm,:,:))';
%         t0 = squeeze(S.temp(mm,:,:))';
%         p0 = abs(squeeze( zrmat(mm,:,:)))';
%         dha(mm,:,:) = sw_gpan(s0,t0,p0)'./g;
%     end
%
%     % relative to bottom
%     dha = bsxfun(@minus,dha(:,:,1),dha);
%
%     S.zeta = dha(:,:,end);
%     S.zeta = S.zeta - min(S.zeta(:));
%     pcolorcen(S.zeta); colorbar; axis square;
   % clear dha

%% old stuff for eddyshelf



%     % use cylindrical co-rodinates
%     [ut1,ur1] = cart2pol(avg1(S.u,2),avg1(S.v,1));

%     fvdx = cumtrapz( xvmat(:,1,1), avg1(f,2)/g .* S.v(:,:,end),1);
%     fudy = cumtrapz( yumat(1,:,1), avg1(f,1)/g .* S.u(:,:,end),2);
%
%     S.zeta(2:end-1,2:end-1) = avg1(avg1( avg1(fvdx,1) - avg1(fudy,2) ,1),2);
%     S.zeta(:,1) = S.zeta(:,2) - ...
%                         (yrmat(:,2,end) - yrmat(:,1,end)).*(S.zeta(:,3)-S.zeta(:,2)) ...
%                             ./(yrmat(:,3,end) - yrmat(:,2,end));
%     S.zeta(:,end) = S.zeta(:,end-1) + ...
%                         (yrmat(:,end,end) - yrmat(:,end-1,end)).*(S.zeta(:,end-1)-S.zeta(:,end-2)) ...
%                                ./(yrmat(:,end-1,end) - yrmat(:,end-2,end));
%     S.zeta(1,:) = S.zeta(2,:) - ...
%                     (xrmat(2,:,end) - xrmat(1,:,end)).*(S.zeta(3,:)-S.zeta(2,:)) ...
%                             ./(xrmat(3,:,end) - xrmat(2,:,end));
%     S.zeta(end,:) = S.zeta(end-1,:) + ...
%                         (xrmat(end,:,end) - xrmat(end-1,:,end)).*(S.zeta(end-1,:)-S.zeta(end-2,:)) ...
%                                ./(xrmat(end-1,:,end) - xrmat(end-2,:,end));
%     clear fvdx fudy

%% geopotential anomaly thermal wind stuff

% % Thermal wind balance. Calculate dynamic height and then integrate.
% % calculate w.r.t bottom
% hp = TCOEF * cumtrapz(squeeze(zrmat(1,1,:)),(S.temp - T0),3);
% vg = g/f0*bsxfun(@rdivide,(hp(2:end,:,:) - hp(1:end-1,:,:)),xrmat(2:end,1,1)-xrmat(1:end-1,1,1));
% ug = -g/f0*bsxfun(@rdivide,(hp(:,2:end,:) - hp(:,1:end-1,:)),squeeze(yrmat(:,2:end,1)-xrmat(:,1:end-1,1)));
%
% % Average to C-grid
% vavgx = (vg(1:end-1,:,:) + vg(2:end,:,:))/2;
% vavgxy = (vavgx(:,1:end-1,:) + vavgx(:,2:end,:))/2;
% % insert extra rows
% S.v = nan(size(xvmat));
% S.v(2:end-1,:,:) = vavgxy;
% S.v(1,:,:) = S.v(2,:,:) - (xvmat(2,:,:)-xvmat(1,:,:)).*(S.v(3,:,:)-S.v(2,:,:))./(xvmat(3,:,:)-xvmat(2,:,:));
% S.v(end,:,:) = S.v(end-1,:,:) + (xvmat(end,:,:)-xvmat(end-1,:,:)).*(S.v(end-1,:,:)-S.v(end-2,:,:))./(xvmat(end-1,:,:)-xvmat(end-2,:,:));
%
% % zeta
% S.zeta = nan([size(xrmat,1) size(xrmat,2)]);
% tmp = f0/g * cumtrapz(xvmat(:,1,1),S.v(:,:,end),1);
% S.zeta(:,2:end-1) = (tmp(:,1:end-1) + tmp(:,2:end))/2;
% S.zeta(:,1) = S.zeta(:,2) - (yrmat(:,2,end) - yrmat(:,1,end)).*(S.zeta(:,3)-S.zeta(:,2))./(yrmat(:,3,end) - yrmat(:,2,end));
% S.zeta(:,end) = S.zeta(:,end-1) + (yrmat(:,end,end) - yrmat(:,end-1,end)).*(S.zeta(:,end-1)-S.zeta(:,end-2))./(yrmat(:,end-1,end) - yrmat(:,end-2,end));
