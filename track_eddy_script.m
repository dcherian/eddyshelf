% on 08 Apr 2013, this script was converted into the function track_eddy.m
% this is backup

% detects eddes using an obuko-weiss parameter threshold

%% read data
dir = 'runs/runeddy-02-1his/';
fname = [dir 'ocean_his.nc'];

start = [1 1 40 1];
count = [Inf Inf 1 Inf];

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

%% SSH threshold detection

% search region for tracking eddies (in addition to detected diameter)
limit_x = 40*1000;
limit_y = 40*1000;
                                      
zeta = roms_read_data(dir,'zeta');
[xr,yr,zr,~,~,~] = roms_var_grid(fname,'zeta');
h = ncread(fname,'h');
t = ncread(fname,'ocean_time');
dt = t(2)-t(1);

zeta = zeta(2:end-1,2:end-1,:);
xr = xr(2:end-1,2:end-1);
yr = yr(2:end-1,2:end-1);

dx = xr(2,1)-xr(1,1);
dy = yr(1,2)-yr(1,1);

sz = size(zeta(:,:,1));

% detect shelfbreak
sbreak = find_shelfbreak(fname);

% remove background flow contribution to zeta
zeta_bg = zeta(:,end,1);

tic;
for tt=1:size(zeta,3)
    if tt == 1,
        mask = ones(sz);
        d_sbreak = Inf;
    else 
        mask = nan*ones(sz);
        lx = eddy.dia(tt-1)/2 + limit_x;
        ly = eddy.dia(tt-1)/2 + limit_y;
        ix1 = find_approx(xr(:,1),eddy.cx(tt-1)-lx,1);
        ix2 = find_approx(xr(:,1),eddy.cx(tt-1)+lx,1);
        iy1 = find_approx(yr(1,:),eddy.cy(tt-1)-ly,1);
        iy2 = find_approx(yr(1,:),eddy.cy(tt-1)+ly,1);
        
        mask(ix1:ix2,iy1:iy2) = 1;
        
        d_sbreak = eddy.cx-sbreak;
    end
    temp = eddy_diag(bsxfun(@minus,zeta(:,:,tt),zeta_bg) .* mask, ...
                        dx,dy,sbreak,w(:,:,tt));
                    
    % [cx,cy] = location of weighted center (first moment)
    eddy.cx(tt)  = temp.cx;
    eddy.cy(tt)  = temp.cy;
    % [mx,my] = location of maximum closest to shelfbreak
    eddy.mx(tt)  = temp.mx;
    eddy.my(tt)  = temp.my;
    % amplitude defined as max - mean amplitude around perimeter
    eddy.amp(tt) = temp.amp;
    % diameter of circle with same area
    eddy.dia(tt) = temp.dia;
    % mask for diagnostic purposes
    eddy.mask(:,:,tt) = temp.mask;
    % number of pixels in eddy
    eddy.n(tt) = temp.n;
end
toc;
disp('Done.');

%% compare tracks of centroid and local maximum

figure;
[C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
clabel(C,hc); hold on
plot(eddy.cx/1000,eddy.cy/1000,'b',eddy.mx/1000,eddy.my/1000,'r'); 
legend('Bathy','Weighted Center','Maximum');

%% movie

hfig = figure;

%var = zeta;
var = double(squeeze(ncread(fname,'temp',[1 1 40 1],[Inf Inf 1 Inf])));
var = var(2:end-1,2:end-1,:); varname = 'temp';
%var = vor; varname = 'surface rel. vor.';

xuz = repmat(xu(:,1), [1 size(zu,3)]);
xrz = repmat(xr(:,1), [1 size(zu,3)]);

%%
figure
for tt=40:size(zeta,3)
    clf
%     subplot(121);
    pcolorcen(xr./1000,yr./1000,zeta(:,:,tt));% .* eddy.mask(:,:,tt));
    shading flat;
    hold on
    caxis([nanmin(zeta(:)) nanmax(zeta(:))]); colorbar
    freezeColors; cbfreeze;   
    [C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
    clabel(C,hc);
    plot(eddy.cx(tt)./1000,eddy.cy(tt)./1000,'x','MarkerSize',16);
    plot(eddy.mx(tt)./1000,eddy.my(tt)./1000,'bo','MarkerSize',16);
    plot(eddy.cx(1:tt)./1000,eddy.cy(1:tt)./1000,'k-');
    contour(xr./1000,yr./1000,eddy.mask(:,:,tt),'k');
    title(['zeta | t =' num2str((tt-1)*dt/86400) ' days']);
    %beautify; 
    xlabel('X (km)'); ylabel('Y (km)');
    axis image
    
%     subplot(122)
%     pcolorcen(xr./1000,yr./1000,var(:,:,tt));
%     shading flat;
%     hold on
%     caxis([nanmin(var(:)) nanmax(var(:))]); colorbar
%     freezeColors; cbfreeze
%     [C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
%     clabel(C,hc);
%     plot(eddy.cx(tt)./1000,eddy.cy(tt)./1000,'x','MarkerSize',16);
%     plot(eddy.cx(1:tt)./1000,eddy.cy(1:tt)./1000,'k-');
%     contour(xr./1000,yr./1000,eddy.mask(:,:,tt),'k');
%     title([varname ' | t =' num2str((tt-1)*dt/86400) ' days']);
%     %beautify; 
%     xlabel('X (km)'); ylabel('Y (km)');
%     axis image
    
%     subplot(133)
%     pcolorcen(xuz./1000,squeeze(zu(:,1,:)),squeeze(u(:,floor(eddy.my(tt)/dy),:,tt)));
%     shading flat; colorbar
%     caxis([-0.15 0.15]);
    
    pause(0.01);
end