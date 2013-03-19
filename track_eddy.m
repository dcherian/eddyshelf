% detects eddes using an obuko-weiss parameter threshold

%% read data
fname = 'runs/runte-11/ocean_his_0011.nc';

start = [1 1 40 1];
count = [Inf Inf 1 10];

thresh = -2e-12;

% u = double(squeeze(ncread(fname,'u',start,count)));
% v = double(squeeze(ncread(fname,'v',start,count)));
zeta = double(squeeze(ncread(fname,'zeta',[start(1) start(2) start(4)], ...
                        [count(1) count(2) count(4)])));

% [xu,yu,~,~,~,~] = roms_var_grid(fname,'u');
% [xv,yv,~,~,~,~] = roms_var_grid(fname,'v');
[xr,yr,~,~,~,~] = roms_var_grid(fname,'zeta');

zeta = zeta(2:end-1,2:end-1,:);
xr = xr(2:end-1,2:end-1);
yr = yr(2:end-1,2:end-1);

dx = xr(2,1)-xr(1,1);
dy = yr(1,2)-yr(1,1);

for tt=1:size(zeta,3)
    temp = eddy_diag(zeta(:,:,tt),dx,dy);
    
    eddy.cx(tt)  = temp.cx;
    eddy.cy(tt)  = temp.cy;
    eddy.amp(tt) = temp.amp;
    eddy.dia(tt) = temp.dia;
    eddy.mask(:,:,tt) = temp.mask;
end

%% movie

h = figure;
for tt=1:size(zeta,3)
    clf
    pcolorcen(xr./1000,yr./1000,zeta(:,:,tt) .* eddy.mask(:,:,tt));
    hold on
    caxis([0 0.05]);
    %freezeColors; cbfreeze
    %[C,h] = contour(xr./1000,yr./1000, w(:,:,tt) * 1e10,[0 1 2 3],'k');
    %clabel(C,h);
    plot(eddy.cx(tt)./1000,eddy.cy(tt)./1000,'x','MarkerSize',16);
    plot(eddy.cx(1:tt)./1000,eddy.cy(1:tt)./1000,'k-');
    title(['zeta (color, meter), tt=' num2str(tt)]);
    beautify; axis image; xlabel('X (km)'); ylabel('Y (km)');
    pause
end