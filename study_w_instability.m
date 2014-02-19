dir1 = 'runs/topoeddy/runew-06-scylla/';
dir2 = 'runs/topoeddy/runew-06-stretch/';
rgrid6 = roms_get_grid([dir1 '/ocean_avg.nc'], [dir1 '/ocean_avg.nc']);
rgrid6str = roms_get_grid([dir2 '/ocean_avg.nc'], [dir2 '/ocean_avg.nc']);

%% assume run6 and run6str are loaded
xind = 70; yind = 55;
w = dc_roms_read_data('runs/topoeddy/runew-06-big/', ...
                'w',[],{'x' xind xind; 'y' yind yind});
wstr = dc_roms_read_data('runs/topoeddy/runew-06-stretch/', ...
                'w',[],{'x' xind xind; 'y' yind yind});


z = rgrid6.z_w(:,yind,xind);
zstr = rgrid6str.z_w(:,yind,xind);


figure;
subplot(311)
plot(avg1(z),diff(z)); hold on
plot(avg1(zstr),diff(zstr),'r');
ylabel('\Delta z (m)');
legend('old','new');
liney(50);
xlim([-2000 0]);

subplot(312);
plot(repmat(z,[1 size(w,2)]),w);
limy = ylim; ylabel('w (m/s)');
xlim([-2000 0]);

subplot(313);
plot(repmat(zstr,[1 size(wstr,2)]),wstr);
ylim(limy); ylabel('w (m/s) - new stretching');
xlabel(' z (m)');
xlim([-2000 0]);
