dir1 = 'runs/topoeddy/runew-06-scylla/';
dir2 = 'runs/topoeddy/runew-06-z60/';
rgrid6 = roms_get_grid([dir1 '/ocean_avg.nc'], [dir1 '/ocean_avg.nc']);
rgrid6str = roms_get_grid([dir2 '/ocean_avg.nc'], [dir2 '/ocean_avg.nc']);

%% assume run6 and run6str are loaded
xind = 140; yind = 60;
w = dc_roms_read_data(dir1, ...
                'w',[],{'x' xind xind; 'y' yind yind});
wstr = dc_roms_read_data(dir2, ...
                'w',[],{'x' xind xind; 'y' yind yind});


z = rgrid6.z_w(:,yind,xind);
zstr = rgrid6str.z_w(:,yind,xind);

%%
figure;
subplot(311)
plot(avg1(z),diff(z)); hold on
plot(avg1(zstr),diff(zstr),'r');
ylabel('\Delta z (m)');
legend('old','z60');
liney(50);
xlim([min([z; zstr]) 0]);

subplot(312);
plot(repmat(z,[1 size(w,2)]),w);
limy = ylim; ylabel('w (m/s)');
xlim([min([z; zstr]) 0]);

subplot(313);
plot(repmat(zstr,[1 size(wstr,2)]),wstr);
%ylim(limy); ylabel('w (m/s) - z60');
xlabel(' z (m)');
xlim([min([z; zstr]) 0]);
