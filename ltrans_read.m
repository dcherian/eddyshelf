track_file = 'ltrans-eddy-04.nc';
out_file = 'runs/runeddy-04-1his/ocean_avg.nc';

lat = ncread(track_file,'lat');
lon = ncread(track_file,'lon');
depth = ncread(track_file,'depth');
age = ncread(track_file,'age')/86400;

h = ncread(out_file,'h');
rgrid = roms_get_grid(out_file,out_file,1,1);

%%

figure
surf(rgrid.x_rho(1,:)/1000,rgrid.y_rho(:,1)/1000,-h');
shading flat
hold on
for i=1:size(lon,1)
    plot3(lon(i,:)/1000,lat(i,:)/1000,depth(i,:),'k');
end
plot3(lon(:,1)/1000,lat(:,1)/1000,depth(:,1),'bo');
xlabel('x'); ylabel('y'); zlabel('z');

%% animate

zeta = double(ncread(out_file,'zeta'));
time = double(ncread(out_file,'ocean_time'))/86400;

%%

% ADD FACTOR
hf = figure;
for i=41:size(zeta,3)
    figure(hf);
    clf
    contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,i)'); shading flat
    axis image
    hold on
    plot(lon(:,1)/1000,lat(:,1)/1000,'bx','MarkerSize',16);
    for j = 1:size(lon,1)
         %plot(lon(j,i)/1000,lat(j,i)/1000,'k.','MarkerSize',16);
         plot(lon(j,1:i*6)/1000,lat(j,1:i*6)/1000,'k.','MarkerSize',16);
    end
    title(['time = ' num2str(time(i)) ' days']);
    pause(0.01)
end