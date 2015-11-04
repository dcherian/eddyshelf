loc = 'data/leebrink2010/';
g3 = load([loc 'globec3.mat']);

g3uv = load([loc 'g3_t2_intensive_uv']);
g3.xmat = repmat([1:size(g3.pres,1)]',[1 size(g3.pres,2)]);

%% 

subplot(3,1,1)
pcolorcen(g3.xmat,g3.pres,g3.theta);
set(gca, 'ydir', 'reverse');

%%
figure
pcolorcen(g3.xmat,g3.pres,g3.sal);
set(gca, 'ydir', 'reverse');

figure
pcolorcen(hypot(g3uv.longitude, g3uv.latitude), g3uv.depth, g3uv.v');
%%
scatter(g3uv.longitude,g3uv.latitude,12,g3uv.water_depth(:,1));

%% representative bathymetry

g3.bathy = smooth(g3uv.water_depth(157:243),5);
plot(-g3.bathy)

%% 3d
figure;
