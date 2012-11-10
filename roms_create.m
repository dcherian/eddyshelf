% Creates initial condition and grid files for ROMS runs
% Modified from Arango's d_initial.m

% TODO
%  - Add beta
%  - Add uz to Ri computation

%% Parameters

% names cannot start with number
FOLDER    = 'runs\';
GRID_NAME = 'test_grd_01';
INI_NAME  = 'test_ini_01';

% Grid Parameters
S.spherical = 0; % 0 - Cartesian, 1 - Spherical

S.Lm = 82;
S.Mm = 80;
S.N  = 40;
S.NT = 1; % Number of active tracers

S.Vtransform = 2;
S.Vstretching = 4;
S.theta_s = 1.0;     %  S-coordinate surface control parameter.
S.theta_b = 2.0;     %  S-coordinate bottom control parameter.
S.Tcline  = 10.0;    %  S-coordinate surface/bottom stretching width (m)

% Domain Extent (in m)
X = 400000;
Y = 400000;
Z = 2000;

% coriolis parameters
lat_ref = 45;
f0    = 2 * (2*pi/86400) *sind(lat_ref);
beta  = 2e-11;

% Physical Parameters
N2    = 1e-5;
T0    = 14;
S0    = 35;
R0    = 1027;
TCOEF = 1.7e-4;
g     = 9.81;

% fix file names
GRID_NAME = [FOLDER GRID_NAME '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER GRID_NAME '.nc'];
INI_NAME = [FOLDER INI_NAME '.nc'];% '-' num2str(ceil(X/1000)) 'x' num2str(ceil(Y/1000)) '-' num2str(S.Lm) 'x' num2str(S.Mm) 'x' num2str(S.N) '.nc'];[FOLDER INI_NAME '.nc'];

%% Create Junk IC & Grid Files

% Set other variables
S.hc = S.Tcline;
S.ncname = INI_NAME;
OA_INTERPOLATE = 0;
INIname = INI_NAME;
GRDname = GRID_NAME;

% Create *new* grid file
c_grid(S.Lm+2,S.Mm+2,GRID_NAME,1);

% Create IC file
[status]=c_initial(S);

%  Set attributes for "ocean_time".

avalue='seconds since 0001-01-01 00:00:00';
[status]=nc_attadd(INIname,'units',avalue,'ocean_time');
  
avalue='360.0 days in every year';
[status]=nc_attadd(INIname,'calendar',avalue,'ocean_time');


%---------------------------------------------------------------------------
%  Set grid variables.
%---------------------------------------------------------------------------

[vname,nvars]=nc_vname(GRDname);

%  Horizontal grid variables. Read in for input GRID NetCDF file.

if (S.spherical),
  S.lon_rho = nc_read(GRDname, 'lon_rho');
  S.lat_rho = nc_read(GRDname, 'lat_rho');
  
  S.lon_u   = nc_read(GRDname, 'lon_u');
  S.lat_u   = nc_read(GRDname, 'lat_u');
  
  S.lon_v   = nc_read(GRDname, 'lon_v');
  S.lat_v   = nc_read(GRDname, 'lat_v');
else,  
  S.x_rho   = nc_read(GRDname, 'x_rho');
  S.y_rho   = nc_read(GRDname, 'y_rho');
  
  S.x_u     = nc_read(GRDname, 'x_u');
  S.y_u     = nc_read(GRDname, 'y_u');
  
  S.x_v     = nc_read(GRDname, 'x_v');
  S.y_v     = nc_read(GRDname, 'y_v');  
end,  

%  Read in Land/Sea mask, if appropriate.

for n=1:nvars,
  name=deblank(vname(n,:));
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

[S.s_rho, S.Cs_r]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N, ...
			     0, 1);

[S.s_w,   S.Cs_w]=stretching(S.Vstretching, ...
                             S.theta_s, S.theta_b, S.hc, S.N, ...
			     1, 1);

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

%% Fix grids and bathymetry

% for future use
xmid = ceil(S.Lm/2);
ymid = ceil(S.Mm/2);
zmid = ceil(S.N/2);

% x,y grids
dx = X/S.Lm;
dy = Y/S.Mm;

S.x_rho = repmat([-dx/2:dx:X+dx/2]',[1 S.Mm+2]);
S.y_rho = repmat([-dy/2:dy:Y+dy/2] ,[S.Lm+2 1]);

S.x_u = repmat([0:dx:X]',[1 S.Mm+2]);
S.y_u = repmat([-dy/2:dy:Y+dy/2],[S.Lm+1 1]);

S.x_v = repmat([-dx/2:dx:X+dx/2]',[1 S.Mm+1]);
S.y_v = repmat([0:dy:Y],[S.Lm+2 1]);

S.x_psi = repmat([0:dx:X]',[1 S.Mm+1]);
S.y_psi = repmat([0:dy:Y],[S.Lm+1 1]);

% Bathmetry
S.h = Z*ones(size(S.h)); % constant depth

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
[z_u]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 3, S.h, S.zeta);
[z_v]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 4, S.h, S.zeta);
[z_w]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 5, S.h, S.zeta);
             
% create grid matrices
% [xrmat,yrmat,zrmat] = ndgrid(S.x_rho(:,1),S.y_rho(:,1),z_r(1,1,:));
% [xumat,yumat,zumat] = ndgrid(S.x_u,S.y_u,z_u);
% [xvmat,yvmat,zvmat] = ndgrid(S.x_v,S.y_v,z_v);

xrmat = repmat(S.x_rho,[1 1 S.N]);
yrmat = repmat(S.y_rho,[1 1 S.N]);
zrmat = z_r;

xumat = repmat(S.x_u,[1 1 S.N]);
yumat = repmat(S.y_u,[1 1 S.N]);
zumat = z_u;

xvmat = repmat(S.x_v,[1 1 S.N]);
yvmat = repmat(S.y_v,[1 1 S.N]);
zvmat = z_v;

zwmat = z_w;

clear z_r z_u z_v z_w

hrmat = repmat(S.h,[1 1 S.N]);

% set initial tracers
S.temp = T0*ones(size(S.temp));
S.salt = S0*ones(size(S.salt));

%% Now set initial conditions - all variables are 0 by default S = S0; T=T0

%%%%%%%%%%%%%%%%% options
perturb_zeta = 0; % add random perturbation to zeta
use_thermal_wind = 0; % use thermal wind to calculate balance.
use_radial = 1; % use radial thermal wind balance
%%%%%%%%%%%%%%%%%

tic

% salt
S.salt = S0*ones(size(S.salt));

% temperature
S.temp = T0*ones(size(S.temp));

% Eddy parameters - all distances in m
eddy.dia = 90*1000; % in m
eddy.depth = 500; % depth below which flow is 'compensated'
eddy.Ncos = 14; % no. of points over which the cosine modulates to zero
eddy.tamp = 0.1; % controls gradient
eddy.a = 2;  % alpha in Katsman et al. (2003)
eddy.cx = X/2; % center of eddy
eddy.cy = Y/2; %        " 

% cylindrical co-ordinates
r0 = eddy.dia/2;
[th,r] = cart2pol((S.x_rho-eddy.cx),(S.y_rho-eddy.cy));
rnorm = r./r0;

% temperature field
exponent =  (eddy.a -1)/eddy.a .* (rnorm.^(eddy.a)); % needed for radial calculations later
etemp = eddy.tamp * exp( -1 * exponent ); % in xy plane
zind = find_approx(zrmat(1,1,:),-1 * eddy.depth,1);
eddy.zprof = [zeros(zind-eddy.Ncos,1); (1-cos(pi * [0:eddy.Ncos]'/eddy.Ncos))/2; ones(S.N-zind-1,1)];

% check eddy profile (normalized)
% subplot(131)
% pcolorcen(S.x_rho/1000,S.y_rho/1000,etemp); 
% axis('square');colorbar; xlabel('x (km)'); ylabel('y (km)');
% subplot(132)
% plot(S.x_rho(:,1)/1000,etemp(:,ymid),'*-');xlabel('x (km)');
% subplot(133)
% plot(eddy.zprof,squeeze(zrmat(1,1,:)),'*-'); ylabel('z (m)');

% assign initial stratification
Tz = N2/g/TCOEF * ones(size(zwmat) - [0 0 2]); % at w points except top / bottom face
S.temp(1,:,:) = T0;
for k=2:size(zrmat,3)
    S.temp(:,:,k) = S.temp(:,:,k-1) + Tz(:,:,k-1).*(zrmat(:,:,k)-zrmat(:,:,k-1));
end

% add eddy temperature perturbation
eddy.temp = bsxfun(@times,etemp,S.temp) .* repmat(permute(eddy.zprof,[3 2 1]),[S.Lm+2 S.Mm+2 1]);
eddy.t0 = max(eddy.temp(:));
S.temp = S.temp + eddy.temp;

if use_radial
    % integrated z profile
    int_zprof = trapz(squeeze(zrmat(1,1,:)), ...
        repmat(permute(eddy.zprof,[3 2 1]),[size(S.zeta,1) size(S.zeta,2) 1]),3);

    S.zeta = -TCOEF * eddy.t0 * int_zprof .* (1-exp(-exponent));
    S.zeta = S.zeta + min(S.zeta(:));
    
    zeta1 = -TCOEF;
    
% CHECK UNITS &  WHY IS PROFILE WEIRD NEAR CENTER. 
% ANS =  r d(theta)/dt NOT d(theta)/dt

    % azimuthal velocity = r d(theta)/dt
    rutz = avg1(bsxfun(@times, eddy.temp, ...
                g*TCOEF* 1./f .* (-exponent./r *eddy.a)),3);
    rut = zeros(size(xrmat));
    for i=2:size(xrmat,3)
        rut(:,:,i) = rut(:,:,i-1) + rutz(:,:,i-1).*(zrmat(:,:,i)-zrmat(:,:,i-1));
    end
            
    uu = -1 * bsxfun(@times,rut, sin(th));
    vv = bsxfun(@times, rut, cos(th));
    
    S.u = avg1(uu,1);
    S.v = avg1(vv,2);
         
end

if use_thermal_wind
    % calculate Tx at v points and Ty and u points
    Txv1 = avg1(avg1(diff(S.temp,1,1)./diff(xrmat,1,1),2),1);
    Txv = [Txv1(1,:,:);Txv1;Txv1(end,:,:)];
    clear Txv1
    
    Tyu1 = avg1(avg1(diff(S.temp,1,2)./diff(yrmat,1,2),2),1);
    Tyu = [Tyu1(:,1,:) Tyu1 Tyu1(:,end,:)];
    clear Tyu1;
    
    % v field
    vz = g*TCOEF * bsxfun(@times,avg1(1./f,2),Txv);
    S.v = zeros(size(xvmat));
    for i=2:size(xvmat,3)
        S.v(:,:,i) = S.v(:,:,i-1) + vz(:,:,i).*(zvmat(:,:,i)-zvmat(:,:,i-1));
    end
    
    % u field
    uz = -g*TCOEF * bsxfun(@times,avg1(1./f,1),Tyu);
    S.u = zeros(size(xumat));
    for i=2:size(xumat,3)
        S.u(:,:,i) = S.u(:,:,i-1) + uz(:,:,i).*(zumat(:,:,i)-zumat(:,:,i-1));
    end
    
    % zeta from dynamic height - does not work
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
end

xind = ymid; yind = ymid; zind = 20;

%%
% Add random perturbation
zeta0 = S.zeta; % save for thermal wind check later
if perturb_zeta == 1
    rng('shuffle');
    perturb = randn(size(S.zeta));
    perturb = perturb./max(abs(perturb(:))) * 10^(-4)*3;
    perturb = perturb - mean(perturb(:));
    S.zeta = S.zeta + perturb; % .* ~(xrmat < x1 | xrmat > x2);
end

% remove mean zeta 
S.zeta = S.zeta - nanmean(S.zeta(:));

% ubar & vbar
for ii=1:S.Lm+2
    for jj=1:S.Mm+2
        if ii~=S.Lm+2
            S.ubar(ii,jj) = trapz(squeeze(zumat(ii,jj,:)),S.u(ii,jj,:),3)./hrmat(ii,jj);
        end
        if jj~=S.Mm+2
            S.vbar(ii,jj) = trapz(squeeze(zvmat(ii,jj,:)),S.v(ii,jj,:),3)./hrmat(ii,jj);
        end
    end
end

toc;

% these are needed later too
dz = squeeze(diff(zvmat(1,1,:)));
VZ  = bsxfun(@rdivide,squeeze(diff(S.v(:,yind,:),1,3)),dz');
UZ  = bsxfun(@rdivide,squeeze(diff(S.u(xind,:,:),1,3)),dz');
UZy = bsxfun(@rdivide,squeeze(diff(S.u(:,yind,:),1,3)),dz'); % for Ri

% setup for pv calculation
grid1.xu = xumat(:,1,1);
grid1.yu = squeeze(yumat(1,:,1)');
grid1.zu = squeeze(zumat(1,1,:));

grid1.xv = xvmat(:,1,1);
grid1.yv = squeeze(yvmat(1,:,1)');
grid1.zv = squeeze(zvmat(1,1,:));

grid1.xr = xrmat(:,1,1);
grid1.yr = squeeze(yrmat(1,:,1)');
grid1.zr = squeeze(zrmat(1,1,:));

rho = R0 - TCOEF * (S.temp-T0);

[pv,xpv,ypv,zpv] = pv_cgrid(grid1,S.u,S.v,rho,f0,R0);

pvmin = min(pv(:));
pvmid = pv(xmid,ymid,zmid);

% figure;
% contourf(xpv,zpv,squeeze(pv(:,ymid,:))',40);
% colorbar;
% title(['PV | PV_{min}/PV_{mid} = ' num2str(pvmin/pvmid)]);
% xlabel('x'); ylabel('z (m)');

%% Sanity Checks

if min(S.temp(:)) < 3
    error('Temperature less than or close to 0.');
end
if max(S.zeta(:)) > 1
    error('Zeta > 1m.');
end

if min(S.salt(:)) == 0
    error('Salt set to zero');
end

if max(abs(S.v(:))) > 1.0 || max(abs(S.u(:))) > 1.0
    fprintf('\n max(u) = %.3f m/s | max(v) = %.3f m/s \n',max((abs(S.v(:)))), max((abs(S.u(:)))));
    input('Really high velocities. Are you sure?');
end

Ri = addnan(abs(fillnan(N2./(avg1(VZ.^2,1) + UZy.^2),Inf)),10);
if min(Ri(:)) <= 0.3
    figure; imagesc(Ri'); title('Ri < 10'); colorbar;
    error('Ri <= 0.3.');
else
    fprintf('\n Min. Ri = %.3f\n\n', nanmin(Ri(:)));
end

%% Check plots
%
make_plot = 1;

if make_plot
    % change axes to km if needed
    fx = 1; fy = 1; lx = '(m)'; ly = '(m)';
    if max(abs(S.x_u(:))) > 3500
        fx = 1000; lx = '(km)';
    end
    if max(abs(S.x_u(:))) > 3500
        fy = 1000; ly = '(km)';
    end
    
    limx = [min(xrmat(:,1,1)) max(xrmat(:,1,1))]./fx;
    limy = [min(yrmat(1,:,1)) max(yrmat(1,:,1))]./fy;
    dx = diff(xrmat(:,yind,zind));
    dy = squeeze(diff(yrmat(xind,:,zind)));
    
    
    % check temperature profile
    figure;
    subplot(221)
    contourf(S.x_rho(:,1)/fx,squeeze(zrmat(1,1,:)),squeeze(S.temp(:,ymid,:))',20);
    colorbar;
    title('temp (y=mid)');
    xlabel(['x' lx]); ylabel('z (m)');
    subplot(222)
    contourf(S.y_rho(1,:)/fy,squeeze(zrmat(1,1,:)),squeeze(S.temp(xmid,:,:))',20);
    colorbar;
    title('temp (x=mid)');
    xlabel(['y' ly]); ylabel('z (m)');
    subplot(223)
    contourf(S.x_rho(:,1)/fx,squeeze(zrmat(1,1,:)),squeeze(eddy.temp(:,ymid,:))',20);
    colorbar;
    title('Eddy temp (y=mid)');
    xlabel(['x' lx]); ylabel('z (m)');
    subplot(224)
    contourf(S.y_rho(1,:)/fy,squeeze(zrmat(1,1,:)),squeeze(eddy.temp(xmid,:,:))',20);
    colorbar;
    title('Eddy temp (x=mid)');
    xlabel(['y' ly]); ylabel('z (m)');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots to check thermal wind
    figure;
    subplot(221)
    plot(avg1(xrmat(:,1,1))./fx,diff(zeta0(:,yind))./dx ./ avg1(f(:,yind),1)* g); hold on
    plot(xvmat(:,1,1)./fx,S.v(:,yind,end),'r');
    xlim(limx); xlabel(['x' lx]);
    legend('\zeta_x','v_{z=0}','Location','Best'); 
    title('Check thermal wind @ surface');

    if use_thermal_wind
        subplot(222)
        plot(avg1(xrmat(:,1,1))./fx,diff(S.temp(:,yind,zind),1,1)./dx ./ avg1(f(:,yind),1) * TCOEF*g);
        hold on
        plot(xvmat(:,1,1)./fx,VZ(:,zind-1),'r');
        plot(xrmat(:,1,1)./fx,Txv(:,yind,zind)* TCOEF*g/f0,'g')
        plot(xvmat(:,1,1)./fx,squeeze(vz(:,yind,zind)),'k');
        legend('T_x^{field} * \alpha g/f','v_z','T_x^{imposed} * \alpha g/f_0','v_z (imposed)','Location','NorthEast');
        xlim(limx); xlabel(['x' lx]);
        title('Check thermal wind @ mid level');
        
        subplot(224)
        plot(avg1(yrmat(1,:,1),2)./fy,-squeeze(diff(S.temp(xind,:,zind,1),1,2))./dy./avg1(f(xind,:),2) * TCOEF*g);
        hold on
        plot(yumat(1,:,1)./fy,UZ(:,zind-1),'r');
        plot(yumat(1,:,1)./fy,-Tyu(xind,:,zind)* TCOEF*g/f0,'g')
        plot(yumat(1,:,1)./fy,uz(xind,:,zind),'k');
        legend('-T_y^{field} * \alpha g/f','u_z','-T_y^{imposed} * \alpha g/f_0','u_z (imposed)','Location','NorthWest');
        xlim(limy); xlabel(['y' ly]);
    end

    subplot(223)
    plot(avg1(yrmat(1,:,1),2)./fy,-diff(zeta0(xind,:))./dy./avg1(f(xind,:),2) .* g); hold on
    plot(yumat(1,:,1)./fy,S.u(xind,:,end),'r');
    xlim(limy); xlabel(['y' ly]);
    legend('\zeta_x','u_{z=0}','Location','Best');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot all fields
    figure;
    subplot(241)
    contourf(S.x_u(:,1)./fx,squeeze(zumat(1,1,:)),squeeze(S.u(:,ymid,:))');
    colorbar;
    title('u');
    xlabel(['x ' lx]); ylabel('z (m)');
    
    subplot(242)
    contourf(S.x_v(:,1)./fx,squeeze(zvmat(1,1,:)),squeeze(S.v(:,ymid,:))');
    colorbar;
    title('v');
    xlabel(['x ' lx]); ylabel('z (m)');
    
    subplot(243)
    contourf(S.x_rho(:,1)./fx,squeeze(zrmat(1,1,:)),squeeze(S.temp(:,ymid,:))',20);
    colorbar;
    title('temp');
    xlabel(['x ' lx]); ylabel('z (m)');
    
    subplot(244)
    contourf(xpv./fx,zpv,squeeze(pv(:,ymid,:))',20);
    colorbar;
    title('PV');
    xlabel(['x ' lx]); ylabel('z (m)');

    subplot(245)
    contourf(S.x_u(:,1)./fx,S.y_u(1,:)./fy,squeeze(S.u(:,:,zmid))');
    colorbar;
    title('u');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    subplot(246)
    contourf(S.x_v(:,1)./fx,S.y_v(1,:)./fy,squeeze(S.v(:,:,zmid))');
    colorbar;
    title('v');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    subplot(247)
    contourf(S.x_rho(:,1)./fx,S.y_rho(1,:)./fy,squeeze(S.temp(:,:,zmid))');
    colorbar;
    title('temp');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;

    subplot(248)
    contourf(S.x_rho(:,1)./fx,squeeze(S.y_rho(1,:))./fy,squeeze(S.zeta(:,:,1))');
    colorbar;
    title('zeta');
    xlabel(['x ' lx]); ylabel(['y ' ly]);
    axis square;
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

fprintf('\n\n Files %s & %s written.\n\n', GRID_NAME,INI_NAME);

%% Grid and time step information

dx = min(dx(:)); dy = min(dy(:));

fprintf('\n\n dx = %.2f m | dy = %.2f m', dx, dy);

DX = sqrt(min(dx)^2 + min(dy)^2);

% From Utility/metrics.F
% Barotropic courant number
dt = 50;
ndtfast = 25;
Cbt = sqrt(g*min(S.h(:))) * dt/ndtfast * sqrt(1/dx^2 + 1/dy^2);
Cbc = sqrt(N2)*min(S.h(:))/pi * dt * sqrt(1/dx^2 + 1/dy^2);
Cbc7 = 7 * dt * sqrt(1/dx^2 + 1/dy^2);

fprintf('\n\n Assuming dt = %.2f, ndtfast = %d, \n\n C_bt = %.3f | C_bc = %.3f  | C_bc7 = %.3f\n\n Min. Ri = %.2f\n\n', dt,ndtfast,Cbt,Cbc,Cbc7,min(Ri(:)));
%fprintf('\n\n (dt)_bt < %.2f s | (dt)_bc < %.2f s\n\n', DX/(sqrt(g*min(S.h(:)))), DX/(sqrt(N2)*min(S.h(:))/pi))

%% old stuff for eddyshelf
%     load 'runs/spinup.mat'
%     S.zeta = zeta;
    
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

%% more old stuff

% v
% v0  = zeros(size(xvmat));            
% % vx  = TCOEF * g/f0 *((Tx2-Tx1)/x1 * (xvmat < x1) - (Tx2-Tx1)/x2 * (xvmat > x2) + 0 * ~(xvmat < x1 | xvmat > x2)).*(zvmat+100);
% vz  = TCOEF * g/f0 * (Tx(:,1:end-1,:,:) + Tx(:,2:end,:,:))/2;
% % % S.v = 0 .* ones(size(xvmat));
% % % xref = 0 * (xvmat < x1) + x1 * (xvmat >+ x1 & xvmat < x2) + x2 * (xvmat >=x2);
% % % S.v = v0 + vx.*(xvmat) + vz.*(zvmat+100);
% 
% S.v = zeros(size(xvmat));
% 
% % for i=2:size(xvmat,1)  
% %     S.v(i,:,:) = S.v(i-1,:,:) + (xvmat(i,:,:) - xvmat(i-1,:,:)) .* vx(i,:,:);
% % end
% for k=2:size(xvmat,3)
%     S.v(:,:,k) = S.v(:,:,k-1) + (zvmat(:,:,k)-zvmat(:,:,k-1)) .* vz(:,:,k);
% end
% for j=1:size(xvmat,2)
% end


