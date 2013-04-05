track_file = 'runs/eddyshelf.nc';
out_file = 'runs/runeddy-01-1his/ocean_his.nc';

lat = ncread(track_file,'lat');
lon = ncread(track_file,'lon');
depth = ncread(track_file,'depth');
age = ncread(track_file,'age');

h = ncread(out_file,'h');
%rgrid = roms_get_grid(out_file,out_file,1,1);

%%

figure
surf(rgrid.x_rho(1,:)/1000,rgrid.y_rho(:,1)/1000,-h');
shading flat
hold on
plot3(lon/1000,lat/1000,depth,'k');
xlabel('x'); ylabel('y'); zlabel('z');