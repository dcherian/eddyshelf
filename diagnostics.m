%% load data

dir = 'runs/runeddy-03-1his/';
fname = [dir 'ocean_his.nc'];

rgrid = roms_get_grid(fname,fname,1,1);
zrmat = permute(rgrid.z_r,[3 2 1]);
time = ncread(fname,'ocean_time')/86400;
h = rgrid.h';
isb = find_approx(h(:,1),130,1);

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
%%
% compare tracks of centroid and local maximum
figure;
[C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
clabel(C,hc); hold on
plot(eddy.cx/1000,eddy.cy/1000,'b',eddy.mx/1000,eddy.my/1000,'r'); 
plot(eddy.mx/1000,eddy.my/1000,'r-'); 
legend('Bathy','Weighted Center','Maximum');

%% temperature anomalies (linearly stratified background state)

tmax = 80;
%temp = squeeze(double(ncread(fname,'temp',[1 1 1 1],[Inf Inf Inf tmax])));
temp = squeeze(double(ncread(fname,'temp')));

% remove background stratification
tanom = bsxfun(@minus,temp,temp(:,end,:,1));

%%
for tt = 1:size(tanom,3)
   %pcolorcen(rgrid.x_rho/1000,rgrid.y_rho/1000,tanom(:,:,end,tt)'); axis image; 
   pcolorcen(squeeze(tanom(isb,:,:,tt))');
   hold on; 
   plot(eddy.my(1:tt)/1250,20,'k.','MarkerSize',16);
   title(['t = ' num2str(time(tt)) ' days']);
   shading flat; colorbar; caxis([-0.2 0.2]);
   pause(0.01);
end

%%
figure
for tt = 1:size(zeta,3)
   %pcolorcen(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,tt)); axis image; 
   %pcolorcen(squeeze(tanom(isb,:,:,tt))');
   pcolorcen(xr/1000,yr/1000,zeta(:,:,tt)); axis image; 
%    hold on; 
%    plot(eddy.my(1:tt)/1250,20,'k.','MarkerSize',16);
%    title(['t = ' num2str(time(tt)) ' days']);
   shading flat; colorbar; caxis([-0.05 0.05]);
   pause(0.01);
end

%% offshore transport calculation

%[sbreak,isb] = find_shelfbreak(fname);

tmax = 160;
u = squeeze(double(ncread(fname,'u',[isb 1 1 1],[1 Inf Inf tmax])));
% vertically integrated cross-shore vel.
U = squeeze(trapz(squeeze(zrmat(isb,1,:)),u,2));

v = squeeze(double(ncread(fname,'v',[isb 1 1 1],[1 Inf Inf tmax])));
% vertically integrated cross-shore vel.
V = squeeze(trapz(squeeze(zrmat(isb,1,:)),v,2));

%% hovmoeller plots - depth integrated velocity

figure
ax(1) = subplot(131);
yhov = repmat(rgrid.y_u(:,1),[1 tmax]);
thov = repmat(time(1:tmax)',[length(rgrid.y_u(:,1)) 1]);
pcolorcen(thov,yhov/1000,U); shading flat
hold on
plot(time(1:tmax),eddy.my(1:tmax)/1000,'k');
plot(time(1:tmax),eddy.ne(1:tmax)/1000,'b-');
plot(time(1:tmax),eddy.se(1:tmax)/1000,'b-');
xlabel(' Time (days) ');
ylabel(' Y (km) ');
title('U');


ax(2) = subplot(132);
yhov = repmat(rgrid.y_v(:,1),[1 tmax]);
thov = repmat(time(1:tmax)',[length(rgrid.y_v(:,1)) 1]);
pcolorcen(thov,yhov/1000,V); shading flat
hold on
plot(time(1:tmax),eddy.my(1:tmax)/1000,'k');
plot(time(1:tmax),eddy.ne(1:tmax)/1000,'b-');
plot(time(1:tmax),eddy.se(1:tmax)/1000,'b-');
xlabel(' Time (days) ');
ylabel(' Y (km) ');
title('V');

subplot(133);
[C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
clabel(C,hc); hold on
plot(eddy.mx/1000,eddy.my/1000,'k-'); 
plot(eddy.we/1000,eddy.my/1000,'b-'); 
plot(eddy.ee/1000,eddy.my/1000,'b-'); 
legend('Bathy','Maximum');
for i=1:10:length(time)
    plot(eddy.mx(i)/1000,eddy.my(i)/1000,'r*');
    text(eddy.mx(i)/1000 + 5,eddy.my(i)/1000,num2str(time(i)));
end

pause
spaceplots


axes(ax(1)); center_colorbar;
axes(ax(2)); colorbar

linkaxes(ax,'xy');


