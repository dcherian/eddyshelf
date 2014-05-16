% Grid Parameters
S.spherical = 0; % 0 - Cartesian, 1 - Spherical

% WikiROMS - Note that there are Lm by Mm computational points.
% If you want to create a grid that's neatly divisible by powers of 2,
% make sure Lm and Mm have those factors.
S.Lm = 412;
S.Mm = 252;
S.N  = 72;

%set value of dx,dy for uniform grid
% also, min dx,dy for telescoped grid.
dx0 = 1000;
dy0 = dx0;

% Domain Extent (in m) - rewritten later
X = S.Lm * dx0;120;
Y = S.Mm * dy0;100;
Z = 2200;

% tracers
S.NPT = 3; % number of passive tracers
S.NT = 2+S.NPT; % total number of tracers

% vertical stretching
S.Vtransform = 2;
S.Vstretching = 4;
S.theta_s = 1.0;     %  S-coordinate surface control parameter.
S.theta_b = 1.0;     %  S-coordinate bottom  control parameter.
S.Tcline  = 500;    %  S-coordinate surface/bottom stretching width (m)

% coriolis parameters
lat_ref = 45;
f0    = 5e-5; 2 * (2*pi/86400) * sind(lat_ref);
beta  = 6e-11;

% Physical Parameters
N2    = 1e-5;
T0    = 60;
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

flags.perturb_zeta = 0; % add random perturbation to zeta
flags.spinup = 0; % if spinup, do not initialize ubar/vbar fields.

flags.front = 0; % create shelfbreak front
flags.eddy  = 1; % create eddy
flags.wind  = 0; % create wind forcing file
flags.floats = 0; % need to figure out float seeding locations?
flags.ubt_initial = 1; % add barotropic velocity to initial condition?
flags.OBC = 1;  % create OBC file and set open boundaries
flags.OBC_from_initial = 1; % copy OBC data from initial condition?
% flags.ubt_deep = 0; % nudge to ubt only in deep water - NOT WORKING
%flags.localize_jet = 0;% buffer around eddy where velocity should exist - NOT NEEDED

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
    OBC.east  = true;           % process eastern  boundary segment
    OBC.south = false;           % process southern boundary segment
    OBC.north = true;            % process northern boundary segment
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Barotropic background flow parameters
flags.fplanezeta = 1; % f-plane solution for zeta (BT vel)
flags.bg_shear = 0;

bg.ubt = NaN; % m/s barotropic velocity
              % if NaN; eddy.nl is used to determine it later
bg.vbt = 0;-0.04; % m/s barotropic velocity
bg.shear_fac = 0.2;
bg.shear = NaN; % set later as bg.shear_fac * max(eddy vorticity)
bg.comment = ['shear = shear_fac * max(eddy vorticity) | ', ...
              'ubt,vbt = whichever is non-zero gets assigned shear ', ...
              'if flags.bg_shear = 0, then ubt/vbt is added (again non-zero)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BATHY
% Bathymetry parameters - all measurements in m
%flags.tanh_bathymetry = 0;
%flags.linear_bathymetry = 1;
%flags.old_bathy = 0;
flags.flat_bottom = 0; % set depth in Z above
flags.crooked_bathy = 0;

bathy.S_sh = 0; % Slope Burger number for shelf
bathy.S_sl = 1.5; % slope Burger number for slope

bathy.H_shelf  = 50;
bathy.L_shelf  = 40 * 1000;
bathy.L_slope  = 50 * 1000;
bathy.axis = 'y'; % CROSS SHELF AXIS
bathy.loc  = 'l'; % h - high end of axis; l - low end
bathy.sl_shelf = bathy.S_sh * f0/sqrt(N2);
bathy.sl_slope = bathy.S_sl * f0/sqrt(N2);
bathy.sl_deep = 0;1/8 * f0/sqrt(N2);

% bathymetry smoothing options
bathy.n_points = 4;
bathy.n_passes = 6;

% curved bathymetry
if flags.crooked_bathy
    bathy.L_shelf2 =  30 * 1000;
    bathy.L_entry  = 200* 1000; % deep water to initialize eddy in
    bathy.L_tilt   = 130 * 1000;
end

bathy.comment = ['H_shelf = depth at coast | L_shelf = shelf width | ' ...
                 'L_slope = slope width | axis = cross-shelf axis for bathy ' ...
                 ' | loc = High/Low end of axis | sl_* = slope of shelf/slope/deep bottom' ...
                 ' L_shelf2 = width for smaller shelf (crooked isobaths) | ' ...
                 ' L_tile = length over which shelf width changes | ' ...
                 ' L_entry = length of smaller shelf | n_points = number ' ...
                 ' of points to smooth over | n_passes = number of smoothing' ...
                 ' passes | isb,xsb,hsb = index,axis loc, depth at shelfbreak' ...
                 ' isl,xsl,hsl = index, axis loc, depth at end of ' ...
                 'continental slope | S_sh, S_sl = slope burger ' ...
                 'numbers for shelf and slope. These are calculated ' ...
                 'later too in nondim.'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRID TELESCOPING
flags.telescoping = 1; % telescope dx,dy

grid.dxmin = dx0;
grid.dymin = dy0;
   
% p - positive side : axis > center
% n - negative side : axis < center
grid.xscalep = 75;
grid.xscalen = 75;
grid.yscalep = 75;
grid.yscalen = 75;

if bathy.axis == 'y'
    grid.dxmax = 3*dx0;
    grid.dymax = 2*dy0;

    % telescope for ix > ixp & ix < ixn
    % similarly for iy > iyp & iy < iyn
    grid.ixp = ceil(0.9*S.Lm);
    grid.ixn = floor(0.1*S.Lm);
    grid.iyp = ceil((bathy.L_shelf + bathy.L_slope)*1.50/dx0);
    grid.iyn = 1;
else
    grid.dxmax = 2*dx0;
    grid.dymax = 3*dy0;

    % telescope for ix > ixp & ix < ixn
    % similarly for iy > iyp & iy < iyn
    grid.ixp = ceil((bathy.L_shelf + bathy.L_slope)*1.50/dx0);
    grid.ixn = 1;
    grid.iyp = ceil(0.9*S.Mm);
    grid.iyn = floor(0.1*S.Mm);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDDY
% eddy momentum balance options
%flags.use_radial    = 1; % use radial
flags.use_gradient  = 1; % use gradient wind balance instead of geostrophic
flags.solidbody_katsman  = 1; % Katsman et al. (2003) solid body core profile?
flags.eddy_zhang = ~flags.solidbody_katsman;
flags.vprof_gaussian = 0;%~flags.eddy_zhang; % eddy is gaussian in vertical?

% Eddy parameters - all distances in m
eddy.Bu     = 1; % ratio of (eddy radius to deformation radius)^2
eddy.nl     = 12.15; % eddy velocity scale / eddy translation velocity
                 % parameter
                 % if bg.ubt = NaN; this is used to determine it later

eddy.dia    = NaN; % sqrt(Bu) * NH/pi/f0 - determined later
eddy.R      = NaN; % radius of max. vel - determined later
eddy.depth  = NaN; % depth below which flow is 'compensated' = Z/2 - determined later
eddy.tamp   = 0.15; % controls gradient
eddy.buffer_sp = 50*1000; % distance from  4.3 (2.3) *r0 to sponge edge
eddy.buffer = NaN;7.5*1000; % distance from start of deep water to 4.3 (2.3) * dia
eddy.cx     = NaN; % if NaN, determined using buffer later
eddy.cy     = NaN;Y/2; %              "
eddy.theta0 = pi/2; % surface phase anomaly from Zhang et al. (2013)
                    % 7/16 * pi for WCR
eddy.comment = ['dia = diameter | depth = vertical scale | tamp = amplitude' ...
                ' of temp. perturbation | a = alpha in Katsman et al. (2003)' ...
                ' | (cx,cy) = (x,y) location of center | (ix,iy) = indices of center | ' ...
                'U = max. azimuthal velocity | R = radius of max. vel (U)' ...
                '| theta0 = surface phase anomaly | buffer = distance from domain edge ' ...
                '/start of deep water to 4.2*r0 (zhang) or 2.2 r0 ' ...
                '(katsman) | nl = U/c in chelton terminology | Bu = ' ...
                'burger number = (eddy.R /  deformation radius)^2 | ' ...
                'Ldef = deformation radius NH/(pi*f0)'];

if flags.solidbody_katsman
    eddy.a      = 3;  % ? in Katsman et al. (2003) - NOT FOR ZHANG PROFILE
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shelfbreak front parameters
% front.LTleft  = 12.5 * 1000; % length scale for temperature (m) - onshore
% front.LTright = 8*1000; % length scale - offshore
% front.LTz     = 100; % Gaussian decay scale in the vertical for temperature
% front.Tx0     = 1.2/60000/SCOEF/R0; % max. magnitude of temperature gradient
% front.comment = ['LTleft = onshore length scale | LTright = offshore length scale' ...
%                  ' LTz = vertical scale | slope = frontal slope | Tx0 = amplitude' ...
%                  ' of gradient'];

front.dRho    = 0.6; % delta Rho across front
front.Lx      = 15 * 1000; % m - horizontal scale
front.Lz      = 80; % m - vertical scale
front.slope   = 100/4000; % non-dimensional - frontal slope
front.Tra     = 'salt';
front.comment = ['Lx = horizontal scale | Tra = tracer var for front | ' ...
                 'Lz = vertical scale | slope = frontal slope | dT = change in' ...
                 'tracer value across front | dRho = change in density across front'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wind stress parameters
wind.tau0 = 0; % set later
wind.ramp = 2; % days
wind.v    = 0.05; % m/s
