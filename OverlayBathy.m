% overlay bathymetry image for JHU APL AVHRR SST

load 'data/etopo2_extract.mat'

%%
ix = find_approx(topo.x, -76, 1):find_approx(topo.x,-66,1);
iy = find_approx(topo.y, 34, 1):find_approx(topo.y, 44, 1);

contour(topo.x(ix), topo.y(iy), topo.z(ix,iy)',[0 0], 'k', 'LineWidth', 2);
hold on
[cc,hh] = contour(topo.x(ix), topo.y(iy), topo.z(ix,iy)',-[100 200 2000], ...
    'k' , 'LineWidth', 1);
clabel(cc,hh);
axis image;