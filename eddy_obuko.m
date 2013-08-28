% detects eddes using an obuko-weiss parameter threshold

%% read data
fname = 'runs/runte-11/ocean_his.nc';

start = [1 1 40 1];
count = [Inf Inf 1 10];

thresh = -2e-12;

u = double(squeeze(ncread(fname,'u',start,count)));
v = double(squeeze(ncread(fname,'v',start,count)));
zeta = double(squeeze(ncread(fname,'zeta',[start(1) start(2) start(4)], ...
                        [count(1) count(2) count(4)])));

[xu,yu,~,~,~,~] = roms_var_grid(fname,'u');
[xv,yv,~,~,~,~] = roms_var_grid(fname,'v');
[xr,yr,~,~,~,~] = roms_var_grid(fname,'zeta');

zeta = zeta(2:end-1,2:end-1,:);
xr = xr(2:end-1,2:end-1);
yr = yr(2:end-1,2:end-1);

%% Obuko Weiss Parameter

% calculate required derivatives and averagse to interior rho points
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

for i=1:size(w,3)
    w(:,:,i) = smooth2a(w(:,:,i),3,3);
end
zmean = mean(mean(zeta,1),2);
mask = fillnan(bsxfun(@lt,w,thresh) & (vor<0) & ...
        (bsxfun(@minus,zeta,zmean) > 0),0);

% find centroid in km
zmasked = repnan(zeta .* mask,0);
xc = squeeze(nanmean(trapz(xr(:,1),bsxfun(@times,zmasked,xr),1) ...
            ./ trapz(xr(:,1),zmasked,1),2)) ./ 1000;
yc = squeeze(nanmean(trapz(yr(1,:),bsxfun(@times,zmasked,yr),2) ...
            ./ trapz(yr(1,:),zmasked,2),1)) ./ 1000;

% get radius (in km)   
dx = xr(2,1)-xr(1,1);
dy = yr(1,2)-yr(1,1);
for i = 1:size(mask,3)
    mask1 = mask(:,:,i);
    npoints = sum(~isnan(mask1(:)));
    area = npoints * dx * dy;
    radius(i) = sqrt(area./pi)/1000;
end
        
%% movie

h = figure;
for tt=1:size(u,3)
    clf
    pcolorcen(xr./1000,yr./1000,zeta(:,:,tt) .* mask(:,:,tt));
    hold on
    caxis([0 0.1]);
    freezeColors; cbfreeze
    [C,h] = contour(xr./1000,yr./1000, w(:,:,tt) * 1e10,[0 1 2 3],'k');
    clabel(C,h);
    plot(xc(tt),yc(tt),'x','MarkerSize',16);
    plot(xc(1:tt),yc(1:tt),'k-');
    title(['zeta (color, meter) and W * 10^{10} in black, tt=' num2str(tt)]);
    beautify; axis image; xlabel('X (km)'); ylabel('Y (km)');
    pause
end