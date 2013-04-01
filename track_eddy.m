% detects eddes using an obuko-weiss parameter threshold

%% read data
dir = 'runs/runeddy-01-1his/';
fname = [dir 'ocean_his.nc'];

start = [1 1 40 1];
count = [Inf Inf 1 10];

%thresh = -2e-12;

% search region for tracking eddies (in addition to detected diameter)
limit_x = 20*1000;
limit_y = 20*1000;
                                      
zeta = roms_read_data(dir,'zeta');
[xr,yr,~,~,~,~] = roms_var_grid(fname,'zeta');
h = ncread(fname,'h');
t = ncread(fname,'ocean_time');
dt = t(2)-t(1);

zeta = zeta(2:end-1,2:end-1,:);
xr = xr(2:end-1,2:end-1);
yr = yr(2:end-1,2:end-1);

dx = xr(2,1)-xr(1,1);
dy = yr(1,2)-yr(1,1);

sz = size(zeta(:,:,1));

for tt=1:size(zeta,3)
    if tt == 1,
        mask = ones(sz);
    else 
        mask = zeros(sz);
        lx = eddy.dia(tt-1)/2 + limit_x;
        ly = eddy.dia(tt-1)/2 + limit_y;
        ix1 = find_approx(xr(:,1),eddy.cx(tt-1)-lx,1);
        ix2 = find_approx(xr(:,1),eddy.cx(tt-1)+lx,1);
        iy1 = find_approx(yr(1,:),eddy.cy(tt-1)-ly,1);
        iy2 = find_approx(yr(1,:),eddy.cy(tt-1)+ly,1);
        
        mask(ix1:ix2,iy1:iy2) = 1;
    end
    temp = eddy_diag(zeta(:,:,tt) .* mask,dx,dy);
    
    eddy.cx(tt)  = temp.cx;
    eddy.cy(tt)  = temp.cy;
    eddy.amp(tt) = temp.amp;
    eddy.dia(tt) = temp.dia;
    eddy.mask(:,:,tt) = temp.mask;
    eddy.n(tt) = temp.n;
end

%% movie

hfig = figure;
for tt=1:size(zeta,3)
    clf
    pcolorcen(xr./1000,yr./1000,zeta(:,:,tt));% .* eddy.mask(:,:,tt));
    shading flat;
    hold on
    caxis([0 nanmax(zeta(:))]); colorbar
    freezeColors; cbfreeze
    [C,hc] = contour(xr./1000,yr./1000,h(2:end-1,2:end-1),[200 500 750 1100],'k');
    clabel(C,hc);
    plot(eddy.cx(tt)./1000,eddy.cy(tt)./1000,'x','MarkerSize',16);
    plot(eddy.cx(1:tt)./1000,eddy.cy(1:tt)./1000,'k-');
    contour(xr./1000,yr./1000,eddy.mask(:,:,tt),'k');
    title(['zeta (color, meter), t =' num2str((tt-1)*dt/86400) ' days']);
    beautify; xlabel('X (km)'); ylabel('Y (km)');
    axis image
    pause(0.01)
end