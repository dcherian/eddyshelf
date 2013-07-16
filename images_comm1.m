%% compare N-S & E-W isobaths

%runns = runs('runs/topoeddy/runteb-04-hires-9');
%runew = runs('runs/topoeddy/runbathysouth-02/ocean_avg_001.nc');

fontsize = [14 14 16];
figure
subplot(121)
hold on
pcolorcen(runns.rgrid.xr/1000,runns.rgrid.yr/1000,runns.bathy.h);
xlabel('X (km)'); ylabel('Y (km)');
plot(runns.eddy.cx/1000,runns.eddy.cy/1000,'Color','b','LineWidth',2);
colorbar
clim = caxis;
plot(runns.eddy.we/1000,runns.eddy.cy/1000,'k');
axis image
title('Bathymetry (color), center (blue), west edge (black)');
beautify(fontsize);

subplot(122)
hold on
pcolorcen(runew.rgrid.xr/1000,runew.rgrid.yr/1000,runew.bathy.h);
xlabel('X (km)'); ylabel('Y (km)');
plot(runew.eddy.cx/1000,runew.eddy.cy/1000,'Color','b','LineWidth',2);
plot(runew.eddy.cx/1000,runew.eddy.se/1000,'k');
axis image
caxis(clim);
title('Bathymetry (color), center (blue), south edge (black)');
beautify(fontsize);

export_fig E:/Work/notes/research/eddyshelf/images/compare-ns-ew.png

%% convergence

% run eddytrackres.m