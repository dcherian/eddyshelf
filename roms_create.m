% Creates initial condition and grid files for ROMS runs
% Modified from Arango's d_initial.m

%% Parameters
% names cannot start with number
FOLDER    = 'runs\';
prefix    = 'te';
GRID_NAME = [prefix '_grd'];
INI_NAME  = [prefix '_ini'];
BRY_NAME  = [prefix '_bry'];
FRC_NAME  = [prefix '_frc'];

% fix file names
GRID_NAME = [FOLDER GRID_NAME '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER GRID_NAME '.nc'];
INI_NAME  = [FOLDER INI_NAME  '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER INI_NAME '.nc'];
BRY_NAME  = [FOLDER BRY_NAME  '.nc'];
FRC_NAME  = [FOLDER FRC_NAME  '.nc'];

% Grid Parameters
S.spherical = 0; % 0 - Cartesian, 1 - Spherical

% WikiROMS - Note that there are Lm by Mm computational points. 
% If you want to create a grid that's neatly divisible by powers of 2, 
% make sure Lm and Mm have those factors.
S.Lm = 400;
S.Mm = 320;
S.N  = 40;

% Domain Extent (in m)
X = S.Lm * 1000;120;
Y = S.Mm * 1000;100;
Z = 1500;

% tracers
S.NPT = 3; % number of passive tracers
S.NT = 2+S.NPT; % total number of tracers

% vertical stretching
S.Vtransform = 2;
S.Vstretching = 4;
S.theta_s = 3.0;     %  S-coordinate surface control parameter.
S.theta_b = 3.0;     %  S-coordinate bottom  control parameter.
S.Tcline  = 100.0;    %  S-coordinate surface/bottom stretching width (m)

% coriolis parameters
lat_ref = 45;
f0    = 2 * (2*pi/86400) * sind(lat_ref);
beta  = 2e-11;

% Physical Parameters
N2    = 1e-5;
T0    = 20;
S0    = 28;
R0    = 1027; % only for EOS purposes
TCOEF = 1.7e-4;
SCOEF = 7.6e-4;
g     = 9.81;
rho0  = 1025; % Boussinesq

%% save physical parameters to structure

phys.rho0 = rho0;
phys.T0 = T0; phys.S0 = S0; phys.R0 = R0;
phys.SCOEF = SCOEF; phys.TCOEF = TCOEF;
phys.N2 = N2; phys.g = g;
phys.lat_ref = lat_ref; phys.f0 = f0; phys.beta = beta;
phys.comment = ['rho0 = Boussinessq approx. | (T0,S0,R0,TCOEF,SCOEF) = linear EOS | ' ...
                '(lat_ref,f0,beta) = coriolis beta plane parameters | ' ...
                ' g = 9.81 m/s2 | N2 = Brunt Vaisalla frequency'];
            
%% Options

calc_pv = 0;

flags.tanh_bathymetry = 0;
flags.linear_bathymetry = 1;
flags.flat_bottom = 0; % set depth in Z above
flags.old_bathy = 0;
flags.crooked_bathy = 0;

flags.perturb_zeta = 0; % add random perturbation to zeta
flags.spinup = 0; % if spinup, do not initialize ubar/vbar fields.

flags.front = 0; % create shelfbreak front
flags.eddy  = 1; % create eddy
flags.wind  = 0; % create wind forcing file
flags.ubt_initial = 0; % add barotropic velocity to initial condition?
flags.fplanezeta = 1; % f-plane solution for zeta (BT vel)
flags.floats = 0; % need to figure out float seeding locations?
flags.bg_shear = 0;

% eddy momentum balance options
flags.use_cartesian = 0; % use cartesian forumlation
flags.use_radial    = 1; % use radial
flags.use_gradient  = 1; % use gradient wind balance instead of geostrophic
flags.solidbody_katsman  = 0; % solid body core profile?
flags.eddy_zhang = ~flags.solidbody_katsman;
flags.vprof_gaussian = ~flags.eddy_zhang; % eddy is gaussian in vertical?

% OBC + barotropic velocity options
flags.OBC = 0;  % create OBC file and set open boundaries
flags.OBC_from_initial = 1; % copy OBC data from initial condition?

flags.comment = ['solidbody_katsman = solid body core profile for eddy (Katsman et al. 2003) | ' ...
    'eddy_zhang = use Zhang et al. (2013) profile | ' ...
    'OBC_from_initial = copy OBC data from IC? | use_gradient = gradient ' ...
    'wind balance instead of geostrophic | use_radial = use expression in' ...
    ' radial instead of cartesian co-ordinates | perturb_zeta = add random' ...
    ' perturbation to initial free surface field | ubt_initial = add background' ...
    ' barotropic velocity field | fplanezeta = f-plane solution for zeta (BT vel) |' ...
    ' vprof_gaussian = if 1, then eddy is Gaussian in vertical'];

% DO NOT CHANGE THIS ORDER
if flags.OBC
    OBC.west  = true;           % process western  boundary segment
    OBC.east  = false;           % process eastern  boundary segment
    OBC.south = false;           % process southern boundary segment
    OBC.north = false;            % process northern boundary segment
end

% flags.ubt_deep = 0; % nudge to ubt only in deep water - NOT WORKING
%flags.localize_jet = 0;% buffer around eddy where velocity should exist - NOT NEEDED

% Barotropic background flow parameters
bg.ubt = 1;0.04; % m/s barotropic velocity
bg.vbt = 0;-0.04; % m/s barotropic velocity
bg.shear_fac = 0.2; 
bg.shear = NaN; % set later as bg.shear_fac * max(eddy vorticity)
bg.comment = ['shear = shear_fac * max(eddy vorticity) | ', ...
              'ubt,vbt = whichever is non-zero gets assigned shear ', ...
              'if flags.bg_shear = 0, then ubt/vbt is added (again non-zero)'];

% Bathymetry parameters - all measurements in m
bathy.H_shelf  = 100;
bathy.L_shelf  = 20 * 1000;
bathy.L_slope  =  50 * 1000;
bathy.axis = 'y'; % CROSS SHELF AXIS
bathy.loc  = 'l'; % h - high end of axis; l - low end
bathy.sl_shelf = 0.0005;
bathy.sl_slope = 0.04;
% bathymetry smoothing options
bathy.n_points = 4;
bathy.n_passes = 6;

% curved bathymetry
bathy.L_shelf2 =  30 * 1000;
bathy.L_entry  = 200* 1000; % deep water to initialize eddy in
bathy.L_tilt   = 130 * 1000;

bathy.comment = ['H_shelf = depth at coast | L_shelf = shelf width | ' ...
                 'L_slope = slope width | axis = cross-shelf axis for bathy ' ...
                 ' | loc = High/Low end of axis | sl_* = slope of shelf/slope' ...
                 ' L_shelf2 = width for smaller shelf (crooked isobaths) | ' ...
                 ' L_tile = length over which shelf width changes | ' ...
                 ' L_entry = length of smaller shelf | n_points = number ' ...
                 ' of points to smooth over | n_passes = number of smoothing' ...
                 ' passes | isb,xsb,hsb = index,axis loc, depth at shelfbreak' ...
                 ' isl,xsl,hsl = index, axis loc, depth at end of continental slope'];

% Eddy parameters - all distances in m
eddy.dia    = 50*1000;
eddy.R      = NaN; % radius of max. vel - determined later
eddy.depth  = NaN; % depth below which flow is 'compensated' = Z/2 - determined later
eddy.tamp   = 0.40; % controls gradient
eddy.a      = 3;  % ? in Katsman et al. (2003) - NOT FOR ZHANG PROFILE
eddy.buffer = 10*1000; % distance from start of deep water (domain edge) to 4.4*r0
                       % should account for sponges
eddy.cx     = NaN; % determined using buffer later
eddy.cy     = NaN; %              "
eddy.theta0 = 0; % surface phase anomaly from Zhang et al. (2013)
%eddy.Ncos  = 10; % no. of points over which the cosine modulates to zero 
eddy.comment = ['dia = diameter | depth = vertical scale | tamp = amplitude' ...
                ' of temp. perturbation | a = alpha in Katsman et al. (2003)' ...
                ' | (cx,cy) = (x,y) location of center | (ix,iy) = indices of center | ' ...
                'U = max. azimuthal velocity | R = radius of max. vel (U)' ...
                '| theta0 = surface phase anomaly | buffer = distance from domain edge ' ...
                '/start of deep water to 4.4*r0'];

% Shelfbreak front parameters
front.LTleft  = 12.5 * 1000; % length scale for temperature (m) - onshore
front.LTright = 8*1000; % length scale - offshore
front.LTz     = 100; % Gaussian decay scale in the vertical for temperature
front.slope   = 100/4000; % non-dimensional
front.Tx0     = 1.2/60000/SCOEF/R0; % max. magnitude of temperature gradient
front.Tra     = 'salt';
front.comment = ['LTleft = onshore length scale | LTright = offshore length scale' ...
                 ' LTz = vertical scale | slope = frontal slope | Tx0 = amplitude' ...
                 ' of gradient'];
             
% wind stress parameters
wind.tau0 = 0; % set later 
wind.ramp = 2; % days
wind.v    = 0.05; % m/s
             
% write parameters to initial conditions .nc file - AT THE END
    
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

if (S.spherical),
  S.lon_rho = nc_read(GRDname, 'lon_rho');
  S.lat_rho = nc_read(GRDname, 'lat_rho');
  
  S.lon_u   = nc_read(GRDname, 'lon_u');
  S.lat_u   = nc_read(GRDname, 'lat_u');
  
  S.lon_v   = nc_read(GRDname, 'lon_v');
  S.lat_v   = nc_read(GRDname, 'lat_v');
else  
  S.x_rho   = nc_read(GRDname, 'x_rho');
  S.y_rho   = nc_read(GRDname, 'y_rho');
  
  S.x_u     = nc_read(GRDname, 'x_u');
  S.y_u     = nc_read(GRDname, 'y_u');
  
  S.x_v     = nc_read(GRDname, 'x_v');
  S.y_v     = nc_read(GRDname, 'y_v');  
end

%  Read in Land/Sea mask, if appropriate.

for n=1:nvars,
  name=char(V.Variables(n).Name);
  switch (name),
    case 'mask_rho'
      S.mask_rho = nc_read(GRDname, 'mask_rho');
    case 'mask_u'
      S.mask_u   = nc_read(GRDname, 'mask_u');
    case 'mask_v'
      S.mask_v   = nc_read(GRDname, 'mask_v');
  end,
end,


%  Bathymetry.

S.h = nc_read(GRDname, 'h');

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

%---------------------------------------------------------------------------
%  Write out grid variables.
%---------------------------------------------------------------------------
			 
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

% salt
S.salt = S0*ones(size(S.salt));

% temperature
S.temp = T0*ones(size(S.temp));

%Fix grids

% for future use
xmid = ceil(S.Lm/2);
ymid = ceil(S.Mm/2);
zmid = ceil(S.N/2);

% x,y grids
dx = X/S.Lm;
dy = Y/S.Mm;

S.x_rho = repmat([-dx/2:dx:X+dx/2]',[1 S.Mm+2]);
S.y_rho = repmat(-dy/2:dy:Y+dy/2 ,[S.Lm+2 1]);

S.x_u = repmat([0:dx:X]',[1 S.Mm+2]);
S.y_u = repmat(-dy/2:dy:Y+dy/2,[S.Lm+1 1]);

S.x_v = repmat((-dx/2:dx:X+dx/2)',[1 S.Mm+1]);
S.y_v = repmat(0:dy:Y,[S.Lm+2 1]);

S.x_psi = repmat([0:dx:X]',[1 S.Mm+1]);
S.y_psi = repmat([0:dy:Y],[S.Lm+1 1]);

xrmat = repmat(S.x_rho,[1 1 S.N]);
yrmat = repmat(S.y_rho,[1 1 S.N]);
xumat = repmat(S.x_u,[1 1 S.N]);
yumat = repmat(S.y_u,[1 1 S.N]);
xvmat = repmat(S.x_v,[1 1 S.N]);
yvmat = repmat(S.y_v,[1 1 S.N]);

% change axes to km if needed
fx = 1; fy = 1; lx = '(m)'; ly = '(m)';
if max(abs(S.x_u(:))) > 3500
    fx = 1000; lx = '(km)';
end
if max(abs(S.x_u(:))) > 3500
    fy = 1000; ly = '(km)';
end

fprintf('\n Initialization - %4.1f MB \n\n', monitor_memory_whos);

%% Bathymetry + Coriolis + more grid stuff

% Coriolis with beta. f = f0 @ y=ymid
fnew = f0*ones(size(S.x_rho));
f = fnew + beta * (S.y_rho - S.y_rho(1,ymid));
clear fnew

if flags.flat_bottom
    S.h = Z * ones(size(S.h));
else
    % linear bathymetry
    if flags.linear_bathymetry == 1

        if flags.crooked_bathy
            [S] = bathy2_x(S,bathy,X,Y);
        else
            [S] = bathy_simple(S,bathy,X,Y,bathy.axis);
        end

        if flags.old_bathy
            bathy.sl_slope2 = bathy.sl_slope;

            % main shelf
            [hx,hy,hdeep] = bathy_crooked(S.x_rho,S.y_rho,bathy,X,Y);

            % second section
            B2 = bathy;
            B2.L_entry = 0 *1000;
            B2.L_shelf = 50*1000;
            [hx2,hy2,~] = bathy_crooked(S.x_rho,S.y_rho,B2,X,Y);

            h1 = hx.*hy .*~(S.x_rho > (X-bathy.L_entry-bathy.L_slope));
            h2 = hx2 .* (S.y_rho > Y-B2.L_shelf) .*(S.x_rho > (X-bathy.L_entry-bathy.L_slope));

        %     subplot(221); imagescnan(hx'); set(gca,'ydir','normal')
        %     subplot(222); imagescnan(hy'); set(gca,'ydir','normal')
        %     subplot(223); imagescnan(hx2');  set(gca,'ydir','normal')
        %     subplot(224); imagescnan(hy2');  set(gca,'ydir','normal')

            % final bathymetry
            S.h = h1+h2;
            S.h(S.h == 0) = hdeep;
        end

        % run smoother   
    %    kernel = [1 2 1];

        for i=1:bathy.n_passes
            for mm = 1:size(S.h,1);
                S.h(mm,:) = smooth(S.h(mm,:),bathy.n_points);
                %S.h(mm,:) = filter(kernel,1,S.h(mm,:));
            end
            for mm = 1:size(S.h,2);
                S.h(:,mm) = smooth(S.h(:,mm),bathy.n_points);
            end
        end  

        % Calculate Burger numbers
        S_sh = bathy.sl_shelf * sqrt(N2)./min(f(:)); % shelf
        S_sl = bathy.sl_slope * sqrt(N2)./min(f(:)); % slope

        % Calculate topographic beta
        b_sh = f0 * bathy.sl_shelf / bathy.H_shelf;
        b_sl = f0 * bathy.sl_slope / max(S.h(:));

    end

    if flags.tanh_bathymetry == 1
        scale = 20000;  
        bathy.L_deep = Y - bathy.L_shelf;
        bathy.H_deep  = Z; 
        S.hflat = Z*ones(size(S.h)); % constant depth
        % y-z profile
        %hy = H_deep + (H_shelf-H_deep)*(1+tanh( ((S.y_rho-L_deep)/20000) ))/2;
        %hx = H_shelf - (H_shelf-H_deep)*(1+tanh( ((S.x_rho-L_entry)/20000) ))/2;

        hx = (1-tanh( ((S.x_rho-L_entry)/scale) ))/2;
        hy = (1+tanh( ((S.y_rho-L_deep)/scale) ))/2;
        S.h = H_deep + (H_shelf-H_deep) * (hx.*hy);
    end
end

% Calculate weird stuff
[S.pm,S.pn,S.dndx,S.dmde] = grid_metrics(S,false);

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

Hz = diff(zwmat,1,3); bsxfun(@rdivide,diff(zwmat,1,3),permute(diff(S.s_w',1,1),[3 2 1]));
Z  = abs(max(S.h(:)));

[rx0,rx1] = stiffness(S.h,zrmat);
bathy_title = sprintf(['\n\n Beckmann & Haidvogel number (r_{x0}) = %f (< 0.2 , max 0.4) \n' ...
            '\t\t\t\t Haney number (r_{x1}) = %f (< 9 , maybe 16)'], rx0,rx1);

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

% find shelfbreak
if ~flags.flat_bottom
    dh2dx2 = diff(hvec,2,1)./avg1(diff(ax_cs,1,1).^2,1);
    [~,bathy.isb] = max(dh2dx2(:));
    bathy.isb = bathy.isb-1;
    bathy.hsb = hvec(bathy.isb);
    bathy.xsb = ax_cs(bathy.isb);

    % find end of slope
    [~,bathy.isl] = min(dh2dx2(:));
    bathy.isl = bathy.isl-1;
    bathy.hsl = hvec(bathy.isl);
    bathy.xsl = ax_cs(bathy.isl);
else
    bathy.xsb = 0;
    bathy.isb = 0;
    bathy.hsb = 0;
end

figure;
subplot(131)
if strcmp(bathy.axis,'x')
    ind = ymid;
    pcolorcen(squeeze(xrmat(:,ind,:))/fx,squeeze(zrmat(:,ind,:)),squeeze(zeros(size(xrmat(:,ind,:)))));
    xlabel(['x ' lx]);
    linex(bathy.xsb/fx,['sb = ' num2str(floor(bathy.hsb)) 'm'],'k');
else 
    ind = xmid;
    pcolorcen(squeeze(yrmat(ind,:,:))/fy,squeeze(zrmat(ind,:,:)),squeeze(zeros(size(xrmat(ind,:,:)))));
    xlabel(['y ' ly]);
    linex(bathy.xsb/fy,['sb = ' num2str(floor(bathy.hsb)) 'm'],'k');
end

title(['Z = ' num2str(Z) ' m']);  ylabel(['z (m)']);
shading faceted; beautify;

subplot(132)
surf(S.x_rho/fx,S.y_rho/fy,-S.h); shading interp
title(bathy_title); colorbar;
xlabel(['x ' lx]); ylabel(['y ' ly]); zlabel('z (m)');
beautify;

subplot(133)
contour(S.x_rho/fx,S.y_rho/fy,-S.h./f,30,'k');
xlabel(['x ' lx]); ylabel(['y ' ly]); title('f/h');
beautify;

spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))

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
Tz = N2/g/TCOEF * ones(size(zwmat) - [0 0 2]); % at w points except top / bottom face
strat = T0.*ones(size(zrmat));
strat(1,:,:) = T0;
%strat_flat = strat;
for k=size(zrmat,3)-1:-1:1
    strat(:,:,k) = strat(:,:,k+1) - Tz(:,:,k).*(zrmat(:,:,k+1)-zrmat(:,:,k));
%    strat_flat(:,:,k) = strat_flat(:,:,k+1) - Tz(:,:,k).*(zrflat(:,:,k+1)-zrflat(:,:,k));
end

S.temp = strat;

clear S.Tz

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
    
    [S.Trax,mask_z] = hor_grad_tracer(axmat,ax_as,ax_cs,zrmat,i_cs,i_as,front,bathy);
%     S.Tz = mask_z .* exp(-(zrmat./front.LTz).^2);
%     S.Tz = S.Tz ./ nanmax(S.Tz(:));
%     S.Tz = N2/g/TCOEF .* (repnan(S.Tz,1));
    S.Traz = 0;%N2/g/coef .* ones(size(zrmat));

    % use chain rule to get gradient on SIGMA LEVEL
    dzdx_s = diff(zrmat,1,i_cs)./diff(axmat,1,i_cs);
    S.Trax_sig = avg1(S.Trax,i_cs) + dzdx_s .* avg1(S.Traz,i_cs);

    % then integrate to make T front
    % top to bottom integration done earlier
    % integrate in horizontal dirn.from left to right using gradients on SIGMA LEVEL
    [S,axmat] = reset_flip(S,axmat);
    [S,axmat] = flip_vars(flip_flag,S,axmat);
    if bathy.loc == 'h'
        for i=size(S.Tra,1)-1:-1:1
            S.Tra(i,:,:) = S.Tra(i+1,:,:) + S.Trax_sig(i,:,:)  .*(axmat(i+1,:,:)-axmat(i,:,:));
        end
    else
        for i=2:size(S.Tra,1)
            S.Tra(i,:,:) = S.Tra(i-1,:,:) + S.Trax_sig(i-1,:,:).*(axmat(i,:,:)-axmat(i-1,:,:));
        end
    end
    S = flip_vars(flip_flag,S);

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
    ssttitle = sprintf('Surface Tra | \\Delta surface Tracer = %.2f C', mean(dsst(:)));

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
    colorbar; title('zeta (front only)');

    % link appropriate axes
    linkaxes(axt,'x');
    linkaxes([axt(1) axt(2) axt(5) axt(6)],'xy');
    linkaxes([axt(4) axt(7)],'xy');
    
    spaceplots(0.03*ones([1 4]),0.05*ones([1 2]))
    
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
        S.salt = S.Tra;
    else
        S.temp = S.Tra;
    end
    % clear some vars
    clear S.Tra S.Trax S.Traz S.Tz
end % if shelfbreak front

fprintf('\n Front - %4.1f MB \n\n', monitor_memory_whos);

%% Now create eddy using temperature

if flags.eddy
    if flags.use_cartesian == 1 && flags.use_radial == 1
        error('Both cartesian & radial formulae specified. Correct it!');
    end
    
    % reset in case I'm debugging
    if ~flags.front
        S.temp = strat; 
        S.u = zeros(size(S.u));
        S.v = zeros(size(S.v));
        S.zeta = zeros(size(S.zeta));               
    end
    
    switch bathy.axis
        case 'x'
            eddy.cx = bathy.xsl+eddy.buffer+4.4*eddy.dia/2;
            eddy.cy = Y-eddy.buffer-4.4*eddy.dia/2;
            if X-eddy.cx-4.4*eddy.dia/2 < eddy.buffer
                warning('Eastern edge - eddy center < buffer');
            else
                fprintf('Distance from eastern edge = %.2f km \n', ...
                    (X-eddy.cx-4.4*eddy.dia/2)/1000);
            end
        case 'y'
            eddy.cx = X-eddy.buffer-4.4*eddy.dia/2; % center of eddy
            eddy.cy = bathy.xsl+eddy.buffer+4.4*eddy.dia/2; %597000; %    "
            if Y-eddy.cy-4.4*eddy.dia/2 < eddy.buffer
                warning('Northern edge - eddy center < buffer');
            else
                fprintf('Distance from northern edge = %.2f km \n', ...
                    (Y-eddy.cy-4.4*eddy.dia/2)/1000);
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

    % temp. field -  z-profile (normalized)
%    % cos profile
%     ind = find_approx(deep_z, -1 * eddy.depth,1);
%     eddy.zprof = [zeros(ind-eddy.Ncos,1); (1-cos(pi * [0:eddy.Ncos]'/eddy.Ncos))/2; ones(S.N-ind-1,1)];
    
    % use half-Gaussian profile & normalize
    if flags.vprof_gaussian
        eddy.zprof = exp(-(eddy.z ./ eddy.depth) .^ 2);
    else % energy in barotropic and BC1 mode  
        eddy.zprof = (1+sin((Z/2 - eddy.z)/Z * pi + eddy.theta0))/2;
        eddy.depth = Z/2;
    end
    %eddy.zprof = eddy.zprof./trapz(eddy.z,eddy.zprof);

    % add eddy temperature perturbation
    eddy.tz = repmat(permute(eddy.zprof,[3 2 1]),[S.Lm+2 S.Mm+2 1]);
    eddy.temp = eddy.tamp * bsxfun(@times,eddy.xyprof,eddy.tz);
    if ~isnan(eddy.temp)
        S.temp = S.temp + eddy.temp;
    end

    % Radial
    if flags.use_radial && max(~isnan(eddy.temp(:)))
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
        S.zeta = S.zeta + TCOEF * trapz(eddy.z, ...
                 bsxfun(@minus,eddy.temp,eddy.temp(eddy.ix,eddy.iy,:)),3);
        S.zeta = S.zeta - min(S.zeta(:));
        
        % Calculate azimuthal velocity shear (r d(theta)/dt)_z using geostrophic balance
        if flags.solidbody_katsman
            dTdr = gamma * rnorm ./ r0 .* (rnorm <= rmnorm) ...
                    + (gamma/2 .* rmnorm^2 + 1) .*  (-(eddy.a-1) ./ r0 .* rnorm.^(eddy.a-1)) ...
                               .* exp(-exponent) .* (rnorm > rmnorm);
        else
            if flags.eddy_zhang
                dTdr = - rnorm./r0 .* exp(-rnorm.^2/2) .* (2- rnorm.^2/2);
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
        S.u = S.u + avg1(eddy.u,1);
        S.v = S.v + avg1(eddy.v,2);
        
        % save max. azimuthal velocity for future   
        [eddy.U,iU] = max(abs(rut(:,eddy.iy,end)));
        eddy.R = r(iU,eddy.iy);
        
        % calculate nondim parameters
        nondim.eddy.Ro = eddy.U ./mean(f(rnorm < r(iU,eddy.iy)))./eddy.R;
        if  nondim.eddy.Ro > 0.25, error('Error: Ro > 0.25'); end
        nondim.eddy.Rh = eddy.U/phys.beta/eddy.R^2;
        nondim.eddy.Bu = f0^2  *eddy.R^2 / N2/Z^2;
        nondim.eddy.Ri = N2./(TCOEF*g*eddy.tamp/f0/eddy.R).^2;
        nondim.eddy.Bu_temp = TCOEF *g * Z * eddy.tamp / f0^2 / eddy.R^2;
        fprintf('\n max. Ro = %.2f | Bu = %.2f | Bu_temp = %.2f | Ri = %.2f | Rh = %.2f\n\n', ....
                nondim.eddy.Ro,nondim.eddy.Bu,nondim.eddy.Bu_temp,nondim.eddy.Ri,nondim.eddy.Rh);
        
        
        % calculate Ro using vorticity
        vor = avg1(diff(eddy.v(:,:,end),1,1)./diff(xrmat(:,:,end),1,1),2) - ...
              avg1(diff(eddy.u(:,:,end),1,2)./diff(yrmat(:,:,end),1,2),1);
        Ro1 = vor./avg1(avg1(f,1),2);
        Ro1 = max(abs(Ro1(:)));
        fprintf('\n Max. Ro (vor/f)  = %.2f \n', Ro1);

    end
% 
%     if flags.use_cartesian % FIX FOR BACKGROUND STATE
%         % calculate Tx at v points and Ty and u points
%         Txv1 = avg1(avg1(diff(S.temp,1,1)./diff(xrmat,1,1),2),1);
%         Txv = [Txv1(1,:,:);Txv1;Txv1(end,:,:)];
%         clear Txv1
% 
%         Tyu1 = avg1(avg1(diff(S.temp,1,2)./diff(yrmat,1,2),2),1);
%         Tyu = [Tyu1(:,1,:) Tyu1 Tyu1(:,end,:)];
%         clear Tyu1;
%         
%         % v field
%         vz = avg1(g*TCOEF * bsxfun(@times,avg1(1./f,2),Txv),3);
%         S.v = zeros(size(xvmat));
%         for i=2:size(xvmat,3)
%             S.v(:,:,i) = S.v(:,:,i-1) + vz(:,:,i-1).*(zvmat(:,:,i)-zvmat(:,:,i-1));
%         end
% 
%         % u field
%         uz = avg1(-g*TCOEF * bsxfun(@times,avg1(1./f,1),Tyu),3);
%         S.u = zeros(size(xumat));
%         for i=2:size(xumat,3)
%             S.u(:,:,i) = S.u(:,:,i-1) + uz(:,:,i-1).*(zumat(:,:,i)-zumat(:,:,i-1));
%         end 
%     end

    % check plots
    xind = ymid; yind = ymid; zind = 20;

    if flags.front
        figure(h_check);
        axt(8) = subplot(248);
    else
        hfeddy = figure;
    end
    contourf(xrmat(:,:,1)./fx,yrmat(:,:,1)./fy,S.zeta,20); shading flat;
    colorbar
    freezeColors;cbfreeze
    hold on
    [C,h] = contour(xrmat(:,:,1)./fx,yrmat(:,:,1)./fy,S.h,...
                floor(linspace(min(S.h(:)),max(S.h(:)),5)),'k');
    clabel(C,h);
    title('Zeta with eddy');
    if flags.front
        linkaxes([axt(4) axt(7) axt(8)],'xy');
    end
    if bathy.axis == 'y'
        liney([bathy.xsl bathy.xsb]/fy);
    else
        linex([bathy.xsl bathy.xsb]/fx);
    end
    axis image;
    
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
    
    % clear variables to save space
%    clear rut rutz dTdr strat r rnorm rfb2 sdisc
    eddy.tz = []; eddy.temp = []; eddy.u = []; eddy.v = [];

    fprintf('\n Eddy - %4.1f MB \n\n', monitor_memory_whos);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add barotropic velocity for advection (OBC initial condition also)

% modify velocity and free surface fields
if flags.ubt_initial == 1
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
                error('fplanezeta not enabled.');
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

fprintf('\n BT vel - %4.1f MB \n\n', monitor_memory_whos);

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

% calculate % error in thermal wind balance
xind = 30; yind = 15; zind = 30;

% these are needed later too - calculated for geostrophic component of
% % velocity
% VZ  = diff(     avg1(bsxfun(@times,vgeo, cos(th)),2),1,3) ./ diff(zvmat,1,3);
% UZ  = diff(-1 * avg1(bsxfun(@times,vgeo, sin(th)),1),1,3) ./ diff(zumat,1,3);

VZ  = diff(S.v,1,3) ./ diff(zvmat,1,3);
UZ  = diff(S.u,1,3) ./ diff(zumat,1,3);
        
dx = squeeze(diff(xrmat(:,:,end),1,1));
dy = squeeze(diff(yrmat(:,:,end),1,2));

dRdx = bsxfun(@rdivide, diff_cgrid(tgrid,S.rho,1) * g/R0, avg1(f,1));
dRdy = bsxfun(@rdivide, diff_cgrid(tgrid,S.rho,2) * g/R0, avg1(f,2));

dzdx = bsxfun(@rdivide, diff(zeta0,1,1), dx);
dzdy = bsxfun(@rdivide, diff(zeta0,1,2), dy);

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
    names = {'shelf water','slope water'};
    try
        dc_roms_passive_tracer(S,names);
    catch ME
        warning('Passive tracer variables already created?');
    end
    
    % assign initial condition
    % important, we do not want grid scale discontinuity.
    % taper everything with gaussian or similar
    
    dye_01 = xrmat;
    dye_02 = yrmat;
    dye_03 = zrmat;
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
        eval(['nc_write(S.ncname,''' vname ''',' vname ',1);']);
    end
    clear dye_01 dye_02 dye_03

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
                   eval([varname ' = ' sprintf('dye_%02d',ii) char(range{mm}) ';']);
                   eval(['nc_write(Sbr.ncname,''' varname ''',' varname ',1);']);
               end
           end
       end
       nc_write(Sbr.ncname,'dye_time',bry_time,1);
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
end

%% non-dimensional parameters

if exist('S_sh','var');nondim.S_sh = S_sh; end
if exist('S_sl','var');nondim.S_sl = S_sl; end
nondim.comment = ['Ro = Rossby number | S_sh = shelf Burger number |' ...
                  'S_sl = slope Burger number | Bu_temp = Bu based on eddy temp perturbation' ...
                  ' | Bu = traditional eddy burger number = NH/fR | Ri = Richardson number | ' ...
                  ' Rh = Rhines number = V/beta/R^2'];

%% Check plots

make_plot = 1;

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

fprintf('\n Started writing files...\n');
tic;
% grid file
ncwrite(GRID_NAME,'spherical',S.spherical);
ncwrite(GRID_NAME,'xl',X);
ncwrite(GRID_NAME,'el',Y);
ncwrite(GRID_NAME,'f',f);
ncwrite(GRID_NAME,'h',S.h);
ncwrite(GRID_NAME, 'x_rho',       S.x_rho);
ncwrite(GRID_NAME, 'y_rho',       S.y_rho);
ncwrite(GRID_NAME, 'x_u',         S.x_u);
ncwrite(GRID_NAME, 'y_u',         S.y_u);
ncwrite(GRID_NAME, 'x_v',         S.x_v);
ncwrite(GRID_NAME, 'y_v',         S.y_v);
ncwrite(GRID_NAME, 'x_psi',       S.x_psi);
ncwrite(GRID_NAME, 'y_psi',       S.y_psi);
ncwrite(GRID_NAME, 'mask_u',    S.mask_u);
ncwrite(GRID_NAME, 'mask_v',    S.mask_v);
ncwrite(GRID_NAME, 'mask_rho',  S.mask_rho);
ncwrite(GRID_NAME, 'mask_psi',  S.mask_psi);

ncwrite(GRID_NAME, 'lon_rho',     S.lon_rho);
ncwrite(GRID_NAME, 'lat_rho',     S.lat_rho);
ncwrite(GRID_NAME, 'lon_u',       S.lon_u);
ncwrite(GRID_NAME, 'lat_u',       S.lat_u);
ncwrite(GRID_NAME, 'lon_v',       S.lon_v);
ncwrite(GRID_NAME, 'lat_v',       S.lat_v);

ncwrite(GRID_NAME, 'pm',       S.pm);
ncwrite(GRID_NAME, 'pn',       S.pn);
ncwrite(GRID_NAME, 'dndx',       S.dndx);
ncwrite(GRID_NAME, 'dmde',       S.dmde);
ncwrite(GRID_NAME, 'angle', zeros(size(S.h)));

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
else
  ncwrite(INIname, 'x_rho',       S.x_rho);
  ncwrite(INIname, 'y_rho',       S.y_rho);
  ncwrite(INIname, 'x_u',         S.x_u);
  ncwrite(INIname, 'y_u',         S.y_u);
  ncwrite(INIname, 'x_v',         S.x_v);
  ncwrite(INIname, 'y_v',         S.y_v);
end,

ncwrite(INIname, 'zeta', S.zeta);
ncwrite(INIname, 'ubar', S.ubar);
ncwrite(INIname, 'vbar', S.vbar);
ncwrite(INIname, 'u',    S.u);
ncwrite(INIname, 'v',    S.v);
ncwrite(INIname, 'temp', S.temp);
ncwrite(INIname, 'salt', S.salt);

% write parameters to IC file
write_params_to_ini(INI_NAME,flags);
write_params_to_ini(INI_NAME,bathy);
write_params_to_ini(INI_NAME,phys);
write_params_to_ini(INI_NAME,bg);
% write_params_to_ini(INI_NAME,ubt,'u_background_barotropic');
% write_params_to_ini(INI_NAME,vbt,'v_background_barotropic');
if flags.front, write_params_to_ini(INI_NAME,front); end
if flags.eddy,  write_params_to_ini(INI_NAME,eddy); end
if flags.wind,  write_params_to_ini(INI_NAME,wind); end
if exist('nondim','var'),  write_params_to_ini(INI_NAME,nondim); end

fprintf('\n\n Files %s | %s ', GRID_NAME,INI_NAME);
if flags.OBC, fprintf('| %s ',BRY_NAME); end
if flags.wind, fprintf('| %s ',FRC_NAME); end
fprintf('written.\n');
toc;

%% Grid and time step information

dx = min(dx(:)); dy = min(dy(:));
fprintf('\n\n=============================  SUMMARY  ========================================');
fprintf(['\n\n (nx,ny,nz) = (%d,%d,%d) | (X,Y,Z) = (%.2f , %.2f , %.2f) m |' ...
        '(dx,dy) = ( %.2f , %.2f) m'], ...
        S.Lm,S.Mm,S.N,max(xrmat(:)), max(yrmat(:)), max(abs(zrmat(:))),dx, dy);

DX = sqrt(min(dx)^2 + min(dy)^2);

% From Utility/metrics.F
% Barotropic courant number
dt = 120;
ndtfast = 35;
Cbt = sqrt(g*max(S.h(:))) * dt/ndtfast * sqrt(1/dx^2 + 1/dy^2);
Cbc = sqrt(N2)*min(S.h(:))/pi * dt * sqrt(1/dx^2 + 1/dy^2);
Cbc7 = 7 * dt * sqrt(1/dx^2 + 1/dy^2);

% print to screen
fprintf('\n\n Beckmann & Haidvogel number = %f (< 0.2 , max 0.4) \n \t\t\t\tHaney number = %f (< 9 , maybe 16)', rx0,rx1);
fprintf('\n\n Assuming dt = %.2f, ndtfast = %d, C_bt = %.3f | C_bc = %.3f  | C_bc7 = %.3f\n\n Min. Ri = %.2f', dt,ndtfast,Cbt,Cbc,Cbc7,min(Ri(:)));

fprintf('\n Bathy Parameters');
fprintf('\n\t\t S_shelf = %.2f | S_slope = %.2f | Beta_shelf = %1.2e | Beta_slope = %1.2e \n', ...
    S_sh,S_sl,b_sh,b_sl);
if flags.wind, cprintf('Red',sprintf('Wind tau0 = %.2e \n\n',wind.tau0)); end
if flags.eddy, 
    fprintf('\n Eddy Parameters: ');
    fprintf('\n\t\t max. Ro = %.2f | Bu = %.2f | Bu_temp = %.2f | Ri = %.2f | Rh = %.2f\n\n', ....
                nondim.eddy.Ro,nondim.eddy.Bu,nondim.eddy.Bu_temp,nondim.eddy.Ri,nondim.eddy.Rh);
    fprintf('\n Deploy float in center of eddy = (%d,%d) \n\n',eddy.ix,eddy.iy);
end
if flags.floats
   xlo = find_approx(xrmat(:,1,1),min(floatx(:))*1000,1); 
   xhi = find_approx(xrmat(:,1,1),max(floatx(:))*1000,1); 
   ylo = find_approx(yrmat(1,:,1),min(floaty(:))*1000,1); 
   yhi = find_approx(yrmat(1,:,1),max(floaty(:))*1000,1); 
   fprintf('\n Float deployment locations : (%d:%d , %d:%d)', xlo,xhi,ylo,yhi);
end
fprintf('\n\n');
%fprintf('\n\n (dt)_bt < %.2f s | (dt)_bc < %.2f s\n\n', DX/(sqrt(g*min(S.h(:)))), DX/(sqrt(N2)*min(S.h(:))/pi))

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

%% old solid body profile - zeta
%             expo = exp( -1 * exponent );
%             z1 = (gamma/2 * rnorm.^2 + 1) .* (rnorm <= rmnorm);
%             z2 = (gamma/2 *rmnorm.^2 + 1) .* expo .* (rnorm > rmnorm);
            
            % hack in matching at boundary - THIS SHOULD NOT BE NEEDED
            %z1 = (z1 - nanmin(fillnan(z1(:),0)) .* (z1 ~= 0));
            % S.zeta = S.zeta + TCOEF*eddy.tamp * int_Tz .* (z1 + z2);      
%             z12 = z1+z2;
%             idx = 250;
%             plot(z1(:,idx)); hold on; 
%             plot(z2(:,idx),'r'); 
%             plot(z12(:,idx),'c'); 
%             plot(S.zeta(:,idx),'k');
%             legend('z1','z2','z12','zeta');

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
