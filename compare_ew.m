run1 = runs('runs/topoeddy/runbathysouth-03-new/'); % no mean flow and farther away
run2 = runs('runs/topoeddy/runbathysouth-02/ocean_avg_001.nc'); % mean flow and closer
% 1.5 x 1.5 km

%%
dt = 16; % days
figure;
subplot(4,1,[1 2])
hold on
plot(run1.eddy.cx/1000, run1.eddy.cy/1000,'b', run2.eddy.cx/1000,run2.eddy.cy/1000,'r');
legend('no mean flow','w/ mean flow + closer + 1.5km','Location','NorthWest');
[cc,hh] = contour(run1.rgrid.xr/1000,run1.rgrid.yr/1000,run1.rgrid.h','k');
clabel(cc,hh);
plot(run1.eddy.cx(1)/1000, run1.eddy.cy(1)/1000,'kx', ...
     run2.eddy.cx(1)/1000, run2.eddy.cy(1)/1000,'kx','MarkerSize',12);
axis image; ylabel('Y (km)'); xlabel('X (km)');
ylim([0 150]);

subplot(4,1,3)
plot(run1.eddy.t,run1.eddy.cy/1000,'b',run2.eddy.t + dt,run2.eddy.cy/1000,'r')
xlabel('Time (days)'); ylabel('y-center (km)');

