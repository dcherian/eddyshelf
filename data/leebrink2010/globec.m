g3 = load('globec3');

g3uv = load('g3_t2_intensive_uv');
g3.xmat = repmat([1:size(g3.pres,1)]',[1 size(g3.pres,2)]);

%% 

subplot(3,1,1)
pcolorcen(g3.xmat,g3.pres,g3.theta)
set(gca, 'ydir', 'reverse');

subplot(3,1,2);

%%
scatter(g3uv.longitude,g3uv.latitude,12,g3uv.water_depth(:,1));

%% representative bathymetry

g3.bathy = smooth(g3uv.water_depth(157:243),5);
plot(g3.bathy)
