dir0 = '/scratch/dcherian/runew-06-scylla/';
dirbig = '/scratch/dcherian/runew-06-big/';
dirsp = '/scratch/dcherian/runew-06-bigsp/';

%rgrid0 = roms_get_grid([dir0 '/ocean_his.nc'],[]);
%rgridbig = roms_get_grid([dirbig '/ocean_his.nc'],[]);
%rgridsp = roms_get_grid([dirsp '/ocean_his.nc'],[]);

IX0 = size(rgrid0.x_rho,2);IY0 = size(rgrid0.y_rho,1);
IXbig = size(rgridbig.x_rho,2);IYbig = size(rgridbig.y_rho,1);
IXsp = size(rgridsp.x_rho,2);IYsp = size(rgridsp.y_rho,1);

zeta0 = dc_roms_read_data(dir0,'zeta',[],{'x' IX0-20 ...
                    IX0-20; 'y' IY0-20 IY0-20},[],rgrid0);
zetabig = dc_roms_read_data(dirbig,'zeta',[],{'x' IX0-20 ...
                    IX0-20; 'y' IY0-20 IY0-20},[],rgridbig);
zetasp = dc_roms_read_data(dirsp,'zeta',[],{'x' IX0-20 ...
                    IX0-20; 'y' IY0-20 IY0-20},[],rgridsp);