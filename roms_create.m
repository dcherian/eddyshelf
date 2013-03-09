% Creates initial condition and grid files for ROMS runs
% Modified from Arango's d_initial.m

%% Parameters
% names cannot start with number
FOLDER    = 'runs\';
GRID_NAME = 'topo_grd_01';
INI_NAME  = 'topo_ini_01';
BRY_NAME  = 'topo_bry_01';

% Grid Parameters
S.spherical = 0; % 0 - Cartesian, 1 - Spherical

% grid information
S.Lm = 102;
S.Mm = 100;
S.N  = 40;

% Domain Extent (in m)
X = 150000;
Y = 200000;
Z = 2000;

% tracers
S.NPT = 0; % number of passive tracers
S.NT = 2+S.NPT; % total number of tracers

% vertical stretching
S.Vtransform = 2;
S.Vstretching = 4;
S.theta_s = 3.0;     %  S-coordinate surface control parameter.
S.theta_b = 2.0;     %  S-coordinate bottom  control parameter.
S.Tcline  = 100.0;    %  S-coordinate surface/bottom stretching width (m)

% coriolis parameters
lat_ref = 45;
f0    = 2 * (2*pi/86400) * sind(lat_ref);
beta  = 2e-11;

% Physical Parameters
N2    = 1e-5;
T0    = 20;
S0    = 35;
R0    = 1027; % only for EOS purposes
TCOEF = 1.7e-4;
SCOEF = 7.6e-4;
g     = 9.81;
rho0  = 1025; % Boussinesq

% fix file names
GRID_NAME = [FOLDER GRID_NAME '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER GRID_NAME '.nc'];
INI_NAME  = [FOLDER INI_NAME '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER INI_NAME '.nc'];
BRY_NAME  = [FOLDER BRY_NAME '.nc'];

%% Options

flags.tanh_bathymetry = 0;
flags.linear_bathymetry = 1;
flags.old_bathy = 0;

flags.perturb_zeta = 0; % add random perturbation to zeta
flags.spinup = 0; % if spinup, do not initialize ubar/vbar fields.

flags.front = 0; % create shelfbreak front
flags.eddy  = 0; % create eddy
flags.ubt_initial = 0; % add barotropic velocity to initial condition?

% momentum balance options
flags.use_cartesian = 0; % cartesian
flags.use_radial    = 1; % use radial
flags.use_gradient  = 1; % use gradient wind balance

% OBC + barotropic velocity options
flags.OBC = 1;  % create OBC file and set open boundaries
flags.OBC_from_initial = 1; % copy OBC data from initial condition?

if flags.OBC
    OBC.west  = false;           % process western  boundary segment
    OBC.east  = false;           % process eastern  boundary segment
    OBC.south = false;           % process southern boundary segment
    OBC.north = true;            % process northern boundary segment
end

% flags.ubt_deep = 0; % nudge to ubt only in deep water - NOT WORKING
%flags.localize_jet = 0;% buffer around eddy where velocity should exist - NOT NEEDED

% Barotropic background flow parameters
ubt = 0;0.05; % m/s barotropic velocity
vbt = -0.05; % m/s barotropic velocity

% Bathymetry parameters - all measurements in m
bathy.H_shelf  = 100;
bathy.L_shelf  = 20 * 1000;
bathy.L_slope  =  20 * 1000;
bathy.axis = 'x'; % CROSS SHELF AXIS
bathy.loc  = 'l'; % h - high end of axis; l - low end
bathy.sl_shelf = 0.0005;
bathy.sl_slope = 0.05;

% curved bathymetry
bathy.L_shelf2 =  0 * 1000;
bathy.L_entry  = 0; %200* 1000; % deep water to initialize eddy in
bathy.L_tilt   = 0; %130 * 1000;

% Eddy parameters - all distances in m
eddy.dia   = 60*1000/2; % in m
eddy.depth = 500; % depth below which flow is 'compensated'
%eddy.Ncos  = 10; % no. of points over which the cosine modulates to zero
eddy.tamp  = 20; % controls gradient
eddy.a     = 2;  % ? in Katsman et al. (2003)
eddy.cx    = X-eddy.dia-20000; % center of eddy
eddy.cy    = Y-eddy.dia-25000; %        " 

% Shelfbreak front parameters
front.LTleft  = 12.5 * 1000; % length scale for temperature (m) - onshore
front.LTright = 8*1000; % length scale - offshore
front.slope   = 500/4000; % non-dimensional
front.Tx0     = -1e-4; % max. magnitude of temperature gradient

%%%%%%%%%%%%%%%%%
    
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


%---------------------------------------------------------------------------
%  Interpolate OA of temperature and salinity from standard levels to
%  model depths.
%---------------------------------------------------------------------------
% 
% if (OA_INTERPOLATE),
% 
%   disp(' ')
%   disp([ 'Interpolating from OA fields, please wait ...']);
%   
%   InpRec = 1;
% 
%   Zoa=nc_read(OAname, 'zout');
% 
%   oa_temp=nc_read(OAname, 'temp', InpRec);
%   oa_salt=nc_read(OAname, 'salt', InpRec);
% 
%   for j=1:Mr,
%     for i=1:Lr,
%       Zroms = squeeze(z_r(i,j,:));
%       Toa   = squeeze(oa_temp(i,j,:));
%       Soa   = squeeze(oa_salt(i,j,:));
% 
%       S.temp(i,j,:) = interp1(Zoa, Toa, Zroms, method);
%       S.salt(i,j,:) = interp1(Zoa, Soa, Zroms, method);
%     end,
%   end,
%   
% end,

% salt
S.salt = S0*ones(size(S.salt));

% temperature
S.temp = T0*ones(size(S.temp));

%% Fix grids

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

%% Bathymetry + Coriolis + more grid stuff

% linear bathymetry
if flags.linear_bathymetry == 1
    
    [S] = bathy_simple(S,bathy,X,Y,bathy.axis);
%   [S] = bathy2_x(S,B,X,Y);
%     ix1 = find_approx(S.x_rho(:,1),X-bathy.L_entry-bathy.L_tilt,1);
%     ix2 = find_approx(S.x_rho(:,1),X-bathy.L_entry,1);
%     iy1 = find_approx(S.y_rho(1,:),Y-bathy.L_shelf,1);
%     iy2 = find_approx(S.y_rho(1,:),Y-bathy.L_shelf2,1);
        
    if flags.old_bathy
        bathy.sl_slope2 = bathy.sl_slope;

        % main shelf
        [hx,hy,hdeep] = bathy(S.x_rho,S.y_rho,B,X,Y);

        % second section
        B2 = B;
        B2.L_entry = 0 *1000;
        B2.L_shelf = 50*1000;
        [hx2,hy2,~] = bathy(S.x_rho,S.y_rho,B2,X,Y);

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

bathy.smoother.n_points = 4;
bathy.smoother.n_passes = 2;
    
    for i=1:bathy.smoother.n_passes
        for mm = 1:size(S.h,1);
            S.h(mm,:) = smooth(S.h(mm,:),bathy.smoother.n_points);
            %S.h(mm,:) = filter(kernel,1,S.h(mm,:));
        end
        for mm = 1:size(S.h,2);
            S.h(:,mm) = smooth(S.h(:,mm),bathy.smoother.n_points);
        end
        % smooth corners a little more
%         for mm=-1*n_points:1*n_points
%             S.h(ix1+mm,:) = smooth(S.h(ix1+mm,:),n_points);
%             S.h(ix2+mm,:) = smooth(S.h(ix2+mm,:),n_points);
%         end
    end  
    
%     subplot(131)
%     plot(-S.h');
%     subplot(132)
%     contour(S.x_rho/1000,S.y_rho/1000,-S.h,40,'k');
%     xlabel('x'); ylabel('y'); zlabel('z')
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

% Coriolis with beta. f = f0 @ y=ymid
fnew = f0*ones(size(S.x_rho));
f = fnew + beta * (S.y_rho - S.y_rho(1,ymid));

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
[zrflat]=set_depth(S.Vtransform, S.Vstretching, ...
    S.theta_s, S.theta_b, S.hc, S.N, ...
    1, max(S.h(:)).*ones(size(S.h)), S.zeta,0);
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

hrmat = repmat(S.h,[1 1 S.N]);
Hz = diff(zwmat,1,3); bsxfun(@rdivide,diff(zwmat,1,3),permute(diff(S.s_w',1,1),[3 2 1]));
Z  = abs(max(S.h(:)));

[rx0,rx1] = stiffness(S.h,zrmat);
bathy_title = sprintf(['\n\n Beckmann & Haidvogel number (r_{x0}) = %f (< 0.2 , max 0.4) \n' ...
            '\t\t\t\t Haney number (r_{x1}) = %f (< 9 , maybe 16)'], rx0,rx1);

figure;
subplot(121)
ind = S.Mm/2;
pcolorcen(squeeze(xrmat(:,ind,:))/fx,squeeze(zrmat(:,ind,:)),squeeze(zeros(size(xrmat(:,ind,:)))));
title(['Z = ' num2str(Z) ' m']); xlabel(['x ' lx]); ylabel(['z (m)']);
shading faceted; beautify;

subplot(122)
surf(S.x_rho/fx,S.y_rho/fy,-S.h); shading interp
title(bathy_title); colorbar;
xlabel(['x ' lx]); ylabel(['y ' ly]); zlabel('z (m)');
beautify;

% save for later use
bathy.h = S.h;

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
strat_flat = strat;
for k=size(zrmat,3)-1:-1:1
    strat(:,:,k) = strat(:,:,k+1) - Tz(:,:,k).*(zrmat(:,:,k+1)-zrmat(:,:,k));
    strat_flat(:,:,k) = strat_flat(:,:,k+1) - Tz(:,:,k).*(zrflat(:,:,k+1)-zrflat(:,:,k));
end

S.temp = strat;

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
    S.Tx = hor_grad_tracer(axmat,ax_as,ax_cs,zrmat,i_cs,i_as,front,bathy);
    S.Tz = N2/g/TCOEF * ones(size(zwmat) - [0 0 1]);

    % use chain rule to get gradient on SIGMA LEVEL
    dzdx_s = diff(zrmat,1,i_cs)./diff(axmat,1,i_cs);
    S.Tx_sig = avg1(S.Tx,i_cs) + dzdx_s .* avg1(S.Tz,i_cs);

    % then integrate to make T front
    % top to bottom integration done earlier
    % integrate in horizontal dirn.from left to right using gradients on SIGMA LEVEL
    [S,axmat] = reset_flip(S,axmat);
    [S,axmat] = flip_vars(flip_flag,S,axmat);
    if bathy.loc == 'h'
        for i=size(S.temp,1)-1:-1:1
            S.temp(i,:,:) = S.temp(i+1,:,:) + S.Tx_sig(i,:,:)  .*(axmat(i+1,:,:)-axmat(i,:,:));
        end
    else
        for i=2:size(S.temp,1)
            S.temp(i,:,:) = S.temp(i-1,:,:) + S.Tx_sig(i-1,:,:).*(axmat(i,:,:)-axmat(i-1,:,:));
        end
    end
    S = flip_vars(flip_flag,S);

    % Make plots to check temperature field
    h_check = figure;
    S = reset_flip(S);
    S = flip_vars(flip_flag,S);
    axt(3) = subplot(243);
    plot(ax_cs/1000,1+S.Tx(:,1)./max(abs(S.Tx(:,1)))); hold on
    plot(ax_cs/1000,1-S.h(:,1)./max(S.h(:)),'r')
    xlabel('cross shelf axis (km)'); 
    legend('Temp gradient','bathy')
    S = flip_vars(flip_flag,S);

    dsst = max(S.temp(:,:,end))-min(S.temp(:,:,end));
    ssttitle = sprintf('SST | \\Delta SST = %.2f C', mean(dsst(:)));

    % check front plots
    if bathy.axis == 'x'
        axt(1) = subplot(241);
        contourf(squeeze(xrmat(:,ymid,:))/1000,squeeze(zrmat(:,ymid,:)),squeeze(S.temp(:,ymid,:)),40);
        xlabel('x (km)'); ylabel('z (m)'); title('Temperature');  colorbar
        axt(2) = subplot(242);
        contourf(squeeze(xrmat(:,ymid,:))/1000,squeeze(zrmat(:,ymid,:)),squeeze(S.Tx(:,ymid,:)),40);
        xlabel('x (km)'); ylabel('z (m)'); title('Temperature Gradient'); colorbar
    else
        axt(1) = subplot(241);
        contourf(squeeze(yrmat(xmid,:,:))/1000,squeeze(zrmat(xmid,:,:)),squeeze(S.temp(xmid,:,:)),40);
        xlabel('y (km)'); ylabel('z (m)'); title('Temperature'); colorbar
        axt(2) = subplot(242);
        contourf(squeeze(yrmat(xmid,:,:))/1000,squeeze(zrmat(xmid,:,:)),squeeze(S.temp(xmid,:,:)),40);
        xlabel('y (km)'); ylabel('z (m)'); title('Temperature'); colorbar
    end
    axt(4) = subplot(244);
    contourf(xrmat(:,:,end)/1000,yrmat(:,:,end)/1000, S.temp(:,:,end),40);
    xlabel('x (km)'); ylabel('y (km)'); title(ssttitle); colorbar

    % calculate velocity field
    vel_shear = zeta_sign * g*TCOEF .* bsxfun(@rdivide,avg1(S.Tx,i_as),avg1(f,i_as));
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
        contourf(squeeze(yvmat(xmid,:,:))./fy,squeeze(zvmat(xmid,:,:)),squeeze(S.v(xmid,:,:)));
        colorbar; title('velocity'); xlabel('y (km)'); ylabel('z (m)');
        axt(6) = subplot(246);
        contourf(squeeze(yvmat(xmid,:,:))./fy,squeeze(zvmat(xmid,:,:)),squeeze(vel_shear(xmid,:,:)));
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
end % if shelfbreak front

%% Now create eddy

if flags.eddy
    if flags.use_cartesian == 1 && flags.use_radial == 1
        error('Both cartesian & radial formulae specified. Correct it!');
    end

    % cylindrical co-ordinates
    r0 = eddy.dia/2;
    [th,r] = cart2pol((S.x_rho-eddy.cx),(S.y_rho-eddy.cy));
    rnorm = r./r0; % normalized radius
    eddy.ix = find_approx(S.x_rho(:,1),eddy.cx,1);
    eddy.iy = find_approx(S.y_rho(1,:),eddy.cy,1);

    % assume that eddy location is in deep water
    eddy.z = squeeze(zrmat(eddy.ix,eddy.iy,:));

    % temp. field - in xy plane (normalized)
    exponent = (eddy.a - 1)/eddy.a .* (rnorm.^(eddy.a)); % needed for radial calculations later
    eddy.xyprof = exp( -1 * exponent ); 
    eddy.xyprof = eddy.xyprof./max(eddy.xyprof(:));

    % temp. field -  z-profile (normalized)
%    % cos profile
%     ind = find_approx(deep_z, -1 * eddy.depth,1);
%     eddy.zprof = [zeros(ind-eddy.Ncos,1); (1-cos(pi * [0:eddy.Ncos]'/eddy.Ncos))/2; ones(S.N-ind-1,1)];
    % use half-Gaussian profile
    eddy.zprof = exp(-(eddy.z ./ eddy.depth) .^ 2);    
    int_zprof = trapz(deep_z,eddy.zprof); % for normalization
    eddy.zprof = eddy.zprof./int_zprof;

    % add eddy temperature perturbation
    eddy.tz = strat_flat .* repmat(permute(eddy.zprof,[3 2 1]),[S.Lm+2 S.Mm+2 1]);
    eddy.temp = eddy.tamp * bsxfun(@times,eddy.xyprof,eddy.tz);
    if ~isnan(eddy.temp)
        S.temp = S.temp + eddy.temp;
    end

    % Radial
    if flags.use_radial && max(~isnan(eddy.temp(:)))
        % integrated z profile of eddy.temp (HAS to include stratification)
        int_Tz = nan([size(xrmat,1) size(xrmat,2)]);
        for i=1:size(eddy.tz,1)
            for j=1:size(eddy.tz,2)
                int_Tz(i,j) = trapz(squeeze(zrflat(i,j,:)),eddy.tz(i,j,:),3);
            end
        end
        
        % SSH calculation is same for gradient wind & geostrophic balance?
        S.zeta = S.zeta + -TCOEF * eddy.tamp * int_Tz .* (1-exp(-exponent));
        S.zeta = S.zeta - min(S.zeta(:));

        % CHECK UNITS &  WHY IS PROFILE WEIRD NEAR CENTER. 
        % ANS =  r d(theta)/dt NOT d(theta)/dt

        % azimuthal velocity = r d(theta)/dt
        rutz = avg1(bsxfun(@times, eddy.temp, ...
                    g*TCOEF* 1./f .* (-exponent./r *eddy.a)),3);
        rut = zeros(size(xrmat));
        for i=2:size(xrmat,3)
            rut(:,:,i) = rut(:,:,i-1) + rutz(:,:,i-1).*(zrmat(:,:,i)-zrmat(:,:,i-1));
        end

        % solve quadratic for vel. if gradient wind balance
        if flags.use_gradient
            rhs = rut;
            rfb2 = r.*f ./ 2;            
            sdisc = sqrt(1 + bsxfun(@times,rhs,2./rfb2));% sqrt(discriminant)
            rut = bsxfun(@times,(-1 + sdisc), rfb2);
            if ~isreal(rut)
                warning(['gradient wind calculated complex v! - ' ...
                         'shifting to cyclostrophic balance']);
                rut = -sqrt(bsxfun(@times,rhs,-2*rfb2)); % need -r for real solutions
                fprintf('\n max. Ro = %.2f \n', max(abs(rut(:)))./f0./r0);
            end
        end
        
        eddy.u = -1 * bsxfun(@times,rut, sin(th));
        eddy.v =      bsxfun(@times,rut, cos(th));

        S.u = S.u + avg1(eddy.u,1);
        S.v = S.v + avg1(eddy.v,2);
    end

    if flags.use_cartesian % FIX FOR BACKGROUND STATE
        % calculate Tx at v points and Ty and u points
        Txv1 = avg1(avg1(diff(S.temp,1,1)./diff(xrmat,1,1),2),1);
        Txv = [Txv1(1,:,:);Txv1;Txv1(end,:,:)];
        clear Txv1

        Tyu1 = avg1(avg1(diff(S.temp,1,2)./diff(yrmat,1,2),2),1);
        Tyu = [Tyu1(:,1,:) Tyu1 Tyu1(:,end,:)];
        clear Tyu1;
        
        % v field
        vz = avg1(g*TCOEF * bsxfun(@times,avg1(1./f,2),Txv),3);
        S.v = zeros(size(xvmat));
        for i=2:size(xvmat,3)
            S.v(:,:,i) = S.v(:,:,i-1) + vz(:,:,i-1).*(zvmat(:,:,i)-zvmat(:,:,i-1));
        end

        % u field
        uz = avg1(-g*TCOEF * bsxfun(@times,avg1(1./f,1),Tyu),3);
        S.u = zeros(size(xumat));
        for i=2:size(xumat,3)
            S.u(:,:,i) = S.u(:,:,i-1) + uz(:,:,i-1).*(zumat(:,:,i)-zumat(:,:,i-1));
        end 
    end

    % check plots
    xind = ymid; yind = ymid; zind = 20;

    if flags.front
        figure(h_check);
        axt(8) = subplot(248);
    else
        figure
    end
    contourf(xrmat(:,:,1)./fx,yrmat(:,:,1)./fy,S.zeta,20); shading flat;
    colorbar
    %freezeColors;cbfreeze
    %hold on
    %[C,h] = contour(xrmat(:,:,1)./fx,yrmat(:,:,1)./fy,S.h,5,'k');
    %clabel(C,h);
    title('Zeta with eddy');
    axis image;
    if flags.front
        linkaxes([axt(4) axt(7) axt(8)],'xy');
    end
    
    xedd = find_approx(xrmat(:,1,1),eddy.cx,1);
    yedd = find_approx(yrmat(1,:,1),eddy.cy,1);
    
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
        if ubt ~=0 && vbt ~=0 
            error('Adding both x and y barotropic velocity in initial condition. WILL NOT WORK!');
        end
        
        if vbt == 0
            S.u = S.u + ubt;
            S.zeta = S.zeta + cumtrapz(squeeze(yrmat(1,:,1)),-f./g * ubt,2);
        else
            S.v = S.v + vbt; 
            S.zeta = S.zeta + cumtrapz(squeeze(xrmat(:,1,1)),f./g * vbt,1);
        end
%    end
end

%% Misc calculations (ubar,vbar,pv) - shouldn't require changes

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
toc;

% setup for pv calculation
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

rho = R0 - TCOEF * (S.temp-T0);

[pv,xpv,ypv,zpv] = pv_cgrid(grid1,S.u,S.v,rho,f,R0);

pvmin = min(pv(:));
pvmid = pv(xmid,ymid,zmid);

% figure;
% contourf(xpv,zpv,squeeze(pv(:,ymid,:))',40);
% colorbar;
% title(['PV | PV_{min}/PV_{mid} = ' num2str(pvmin/pvmid)]);
% xlabel('x'); ylabel('z (m)');

%% Sanity Checks + check thermal wind

% calculate % error in thermal wind balance
xind = 30; yind = 15; zind = 30;

% these are needed later too
VZ  = diff(S.v,1,3) ./ diff(zvmat,1,3);
UZ  = diff(S.u,1,3) ./ diff(zumat,1,3);
        
dx = squeeze(diff(xrmat,1,1));
dy = squeeze(diff(yrmat,1,2));

dTdx = bsxfun(@rdivide, diff_cgrid(tgrid,S.temp,1) *  TCOEF*g, avg1(f,1));
dTdy = bsxfun(@rdivide, diff_cgrid(tgrid,S.temp,2) * -TCOEF*g, avg1(f,2));

dzdx = bsxfun(@rdivide, diff(zeta0,1,1), dx(:,:,end));
dzdy = bsxfun(@rdivide, diff(zeta0,1,2), dy(:,:,end));

% PERCENTAGE error in thermal wind balance
diff_u = (avg1(UZ,2) - avg1(avg1(dTdy,1),3))./max(abs(UZ(:))) * 100;
diff_v = (avg1(VZ,1) - avg1(avg1(dTdx,2),3))./max(abs(VZ(:))) * 100;
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
pcolorcen(diff_v(:,:,S.N/2)'); colorbar
title('% error (v_z - dT/dx * \alpha g/f)|_{z=mid}');
subplot(222)
pcolorcen(squeeze(diff_u(:,:,S.N/2))'); colorbar
title('% error (-u_z - dT/dy * \alpha g/f)|_{z=mid}');
subplot(223)
pcolorcen(((avg1((dzdx ./ avg1(f,1)* g),2)-avg1(S.v(:,:,end),1)))./ ...
                max(abs(S.v(:)))*100); colorbar
title('% error (g/f d\zeta /dx - v_{z=0})');
subplot(224)
pcolorcen((avg1(dzdy ./ avg1(f,2)* g,1) + avg1(S.u(:,:,end),2)) ./ ...
            max(abs(S.u(:)))*100); colorbar
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
    figure; imagesc(addnan(Ri',10)); title('Ri < 10'); colorbar;
    error('Ri <= 0.3.');
else
    fprintf('\n Min. Ri = %.3f\n\n', nanmin(Ri(:)));
end 

%% Passive Tracers - Initial condition + writing

if S.NPT > 0
    % create variables first
    names = {'initial x';'initial y';'initial z'};
    try
        dc_roms_passive_tracer(S,names);
    catch ME
        warning('Passive tracer variables already created?');
    end
    
    % assign initial condition
    dye_01 = xrmat;
    dye_02 = yrmat;
    dye_03 = zrmat;    
    
    % write to file
    for ii=1:S.NPT
        vname = sprintf('dye_%02d',ii);
        eval(['nc_write(S.ncname,''' vname ''',' vname ',1);']);
    end
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
    VarList = {'zeta'; 'ubar'; 'vbar';'temp';'salt'}; % variables to write
    
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
       for ii=1:S.NPT
           varname = sprintf('dye_%s_%02d',char(boundaries(mm,:)),ii);
           eval([varname ' = ' sprintf('dye_%02d',ii) char(range{mm}) ';']);
           eval(['nc_write(Sbr.ncname,''' varname ''',' varname ',1);']);
       end
       nc_write(Sbr.ncname,'dye_time',bry_time,1);
    end
end

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
    contourf(squeeze(yumat(xmid,:,:))./fx,squeeze(zumat(xmid,:,:)),squeeze(S.u(xmid,:,:)),20);
    colorbar;
    title('u');
    xlabel(['y ' lx]); ylabel('z (m)');
    
    ax(2) = subplot(242);
    contourf(squeeze(yvmat(xmid,:,:))./fx,squeeze(zvmat(xmid,:,:)),squeeze(S.v(xmid,:,:)),20);
    colorbar;
    title('v');
    xlabel(['y ' lx]); ylabel('z (m)');
    
    ax(3) = subplot(243);
    contourf(squeeze(yrmat(xmid,:,:))./fx,squeeze(zrmat(xmid,:,:)),squeeze(S.temp(xmid,:,:)),20);
    colorbar;
    title('temp');
    xlabel(['y ' lx]); ylabel('z (m)');
    
    % fix pv script so that it returns proper co-ordinates
    ax(4) = subplot(244);
    contourf(repmat(xpv(:,1),1,size(pv,3))./fx,squeeze(zpv(:,ymid,:)),squeeze(pv(:,ymid,:)),20);
    colorbar;
    title('PV');
    xlabel(['x ' lx]); ylabel('z (m)');
    
    % need to interpolate to constant z - surface
    ax(5) = subplot(245);
    contourf(xumat(:,:,end)./fx,yumat(:,:,end)./fy,squeeze(S.u(:,:,end)));
    colorbar;
    title('surface u');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    ax(6) = subplot(246);
    contourf(xvmat(:,:,end)./fx,yvmat(:,:,end)./fy,squeeze(S.v(:,:,end)));
    colorbar;
    title('surface v');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    ax(7) = subplot(247);
    contourf(xrmat(:,:,end)./fx,yrmat(:,:,end)./fy,squeeze(S.temp(:,:,end)));
    colorbar;
    title('SST');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    ax(8) = subplot(248);
    contourf(S.x_rho(:,1)./fx,squeeze(S.y_rho(1,:))./fy,squeeze(S.zeta(:,:,1))');
    colorbar;
    title('SSH (zeta)');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;
    
    linkaxes([ax(1:4)],'xy');
    linkaxes([ax(5:8)],'xy');
end

%% Write to Grid & IC file

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

fprintf('\n\n Files %s | %s ', GRID_NAME,INI_NAME);
if flags.OBC, fprintf('| %s ',BRY_NAME); end
fprintf('written.\n\n');

%% Grid and time step information

dx = min(dx(:)); dy = min(dy(:));

fprintf('\n\n X = %.2f | Y = %.2f | Z = %.2f | dx = %.2f m | dy = %.2f m', ...
                max(xrmat(:)), max(yrmat(:)), max(abs(zrmat(:))),dx, dy);

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
fprintf('\n\n Assuming dt = %.2f, ndtfast = %d, \n\n C_bt = %.3f | C_bc = %.3f  | C_bc7 = %.3f\n\n Min. Ri = %.2f\n\n', dt,ndtfast,Cbt,Cbc,Cbc7,min(Ri(:)));
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
%     dT = abs(Tgrad-avg1(S.Tx,i_as))./max(abs(S.Tx(:))) * 100;
%     yind = S.Mm/2;
%     figure
%     ax(1) = subplot(131);
%     contourf(squeeze(xrmat(:,yind,:)), squeeze(zrmat(:,yind,:)), squeeze(S.Tx(:,yind,:)));
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
