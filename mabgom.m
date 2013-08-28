link = 'http://tds.marine.rutgers.edu/thredds/dodsC/roms/mabgom/v6/avg';

%% get grid

ocean_time = ncread(link,'ocean_time');
xr = ncread(link,'lon_rho');
yr = ncread(link,'lat_rho');
xu = ncread(link,'lon_u');
yu = ncread(link,'lat_u');
xv = ncread(link,'lon_v');
yv = ncread(link,'lat_v');
mr = ncread(link,'mask_rho');
mu = ncread(link,'mask_rho');
mv = ncread(link,'mask_v');
h  = ncread(link,'h');
f  = ncread(link,'f');

S.Vtransform = ncread(link,'Vtransform');
S.Vstretching = ncread(link,'Vstretching');
S.hc = ncread(link,'hc');
S.N  = ncread(link,'N');
S.theta_s = ncread(link,'theta_s');
S.theta_b = ncread(link,'theta_b');

S.zeta = zeros(size(S.h));

% make z grid
[zr]=set_depth(S.Vtransform, S.Vstretching, ...
                S.theta_s, S.theta_b, S.hc, S.N, ...
                1, h, S.zeta,0);
[zu]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 3, h, S.zeta,0);
[zv]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 4, h, S.zeta,0);
[zw]=set_depth(S.Vtransform, S.Vstretching, ...
                 S.theta_s, S.theta_b, S.hc, S.N, ...
                 5, h, S.zeta,0);
             
save mabgom-grid.mat S ocean_time xr yr zr xu yu zu xv yv zv

%%