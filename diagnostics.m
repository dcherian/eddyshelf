%% load data

dir = 'runs/runeddy-02-1his/';
fname = [dir 'ocean_his.nc'];

rgrid = roms_get_grid(fname,fname,1,1);
zrmat = permute(rgrid.z_r,[3 2 1]);
time = ncread(fname,'ocean_time')/86400;
h = rgrid.h';
isb = find_approx(h(:,1),200,1);

%% calculate Obuko-Weiss parameter & relative vorticity

u = double(squeeze(ncread(fname,'u',start,count)));
v = double(squeeze(ncread(fname,'v',start,count)));

[xu,yu,zu,~,~,~] = roms_var_grid(fname,'u');
[xv,yv,zv,~,~,~] = roms_var_grid(fname,'v');

% calculate required derivatives and average to interior rho points
ux = bsxfun(@rdivide,diff(u,1,1),diff(xu,1,1));
uy = avg1(avg1(bsxfun(@rdivide,diff(u,1,2),diff(yu,1,2)),1),2);
vx = avg1(avg1(bsxfun(@rdivide,diff(v,1,1),diff(xv,1,1)),1),2);
vy = bsxfun(@rdivide,diff(v,1,2),diff(yv,1,2));

ux = ux(:,2:end-1,:);
vy = vy(2:end-1,:,:);

% obuko weiss parameter
w = 4 * (vx.*uy - ux.*vy) + (ux+vy).^2;
vor = vx-uy;
thresh = 0.2 * std(std(w)); % Isern-Fontanet et al. (2006)

%% get eddy tracks

eddy = track_eddy(dir);

% compare tracks of centroid and local maximum
figure;
[C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
clabel(C,hc); hold on
plot(eddy.cx/1000,eddy.cy/1000,'b',eddy.mx/1000,eddy.my/1000,'r'); 
legend('Bathy','Weighted Center','Maximum');

%% temperature anomalies (linearly stratified background state)

temp = roms_read_data(dir,'temp');

% figure out background stratification
tedge = squeeze(temp(end,end,:,1));
tanom = bsxfun(@minus,temp,temp(:,end,:,1));

%%
for tt = 1:size(tanom,4)
   %pcolorcen(rgrid.x_rho/1000,rgrid.y_rho/1000,tanom(:,:,end,tt)'); axis image; 
   pcolorcen(squeeze(tanom(isb,:,:,tt))');
   title(['t = ' num2str(time(tt)) ' days']);
   shading interp; colorbar; caxis([-0.2 0.2]);
   pause(0.01);
end

%% offshore transport calculation

%[sbreak,isb] = find_shelfbreak(fname);

tmax = 80;
u = squeeze(double(ncread(fname,'u',[isb 1 1 1],[1 Inf Inf tmax])));
% vertically integrated cross-shore vel.
U = squeeze(trapz(squeeze(zrmat(isb,1,:)),u,2));

%% hovmoeller plot - offshore transport

yhov = repmat(rgrid.y_u(:,1),[1 tmax]);
thov = repmat(time(1:tmax)',[length(rgrid.y_u(:,1)) 1]);
figure
pcolorcen(thov,yhov/1000,U); shading flat
hold on
plot(time(1:tmax),eddy.my(1:tmax)/1000,'k');
plot(time(1:tmax),(eddy.my(1:tmax)+eddy.dia(1:tmax)/2)/1000,'b-');
plot(time(1:tmax),(eddy.my(1:tmax)-eddy.dia(1:tmax)/2)/1000,'b-');