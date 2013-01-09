file = 'runs/ocean_rst.nc';

[xrmat,yrmat,zrmat,~,~,~] = dc_roms_var_grid(file,'rho');

ncwrite(file,'dye_01',repmat(xrmat,[1 1 1 2 2]));
ncwrite(file,'dye_02',repmat(yrmat,[1 1 1 2 2]));
ncwrite(file,'dye_03',repmat(zrmat,[1 1 1 2 2]));

