%% load data

dir = 'runs/runeddy-04-1his/';
fname = [dir 'ocean_avg.nc'];

rgrid = roms_get_grid(fname,fname,1,1);
zrmat = permute(rgrid.z_r,[3 2 1]);
time = ncread(fname,'ocean_time')/86400;
h = rgrid.h';
xr = rgrid.x_rho';
yr = rgrid.y_rho';

% locate shelfbreak
isb = find_approx(h(:,1),120,1);
xsb = xr(isb,1);

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
[C,hc] = contour(xr./1000,yr./1000,h,[200 500 750 1100],'k');
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

%% quiver movie

usurf = double(squeeze(ncread(fname,'u',[1 1 40 1],[Inf Inf 1 Inf])));
vsurf = double(squeeze(ncread(fname,'v',[1 1 40 1],[Inf Inf 1 Inf])));

% subtract out inflow - doesn't work very well
%usurf = bsxfun(@minus,usurf(:,:,:),usurf(:,end,:));
%vsurf = bsxfun(@minus,vsurf(:,:,:),vsurf(:,end,:));

% average to interior RHO points
usurf = avg1(usurf(:,2:end-1,:),1);
vsurf = avg1(vsurf(2:end-1,:,:),2);

% decimate for quiver plot
xskip = 4;
yskip = 10;
% buffer points offshore of shelbreak to plot
buffer = 8;

limx = [min(rgrid.x_rho(:)) max(rgrid.x_rho(:))]/1000;
limy = [min(rgrid.y_rho(:)) max(rgrid.y_rho(:))]/1000;

%mask = hypot(usurf,vsurf) < 0.2;
mask = zeros(size(usurf));
mask(1:isb+buffer,:,:) = 1;
for ii=1:size(usurf,3)
    clf
    quiver(rgrid.x_rho(2:yskip:end-1,2:xskip:end-1)'/1000,rgrid.y_rho(2:yskip:end-1,2:xskip:end-1)'/1000, ...
                usurf(1:xskip:end,1:yskip:end,ii).*mask(1:xskip:end,1:yskip:end,ii), ...
                vsurf(1:xskip:end,1:yskip:end,ii).*mask(1:xskip:end,1:yskip:end,ii),2);
    %axis image
    hold on;
    plot_eddy_track(eddy,ii);
    xlim(limx); ylim(limy);
    title(['t = ' num2str(ii) ' days']);
    %axis image
    pause(0.05);
end

%% quiver at shelfbreak movie

index = isb+6;

usb = double(squeeze(ncread(fname,'u',[index 1 1 1],[1 Inf Inf Inf])));
vsb = double(squeeze(ncread(fname,'v',[index 1 1 1],[1 Inf Inf Inf])));

% average to interior RHO points
usb = usb(2:end-1,:,:);
vsb = avg1(vsb,1);

% create grid matrices
yzmat = repmat(yr(1,2:end-1)',[1 size(zrmat,3)])/1000;
zymat = squeeze(zrmat(index,2:end-1,:));
fakew = zeros(size(usb(:,:,1)));

%%
yskip = 20; zskip = 6;

%az = 21.5; el = 34;
%az = 36.5; el = 42;
az = 0; el=90; % 2d plan view

yrange = 1:yskip:size(usb,1)-1;
zrange = 1:yskip:size(usb,2);

figure
for ii=20:size(usb,3)-5
    clf
    quiver3(xr(index,1)/1000 * ones(size(yzmat(yrange,zrange))), ...
                yzmat(yrange,zrange),zymat(yrange,zrange), ...
                usb(yrange,zrange,ii), ...
                vsb(yrange,zrange,ii),fakew(yrange,zrange),'x');
    hold on;
    surf(xr/1000,yr/1000,-rgrid.h'); shading flat
    %plot3(eddy.mx(1:ii)/1000,eddy.my(1:ii)/1000,zeros([1 ii]),'k-');
    plot_eddy_track(eddy,ii);
    xlabel('x (km)'); ylabel('y (km)'); zlabel('z (m)');
    title(['t = ' num2str(time(ii)) ' days']);
    view(az,el);
    if el == 90, axis image; end
    pause(0.03)
end

%% coneplot
xbuf = 20; % get xbuf*2 points about shelfbreak
ybuf = 100; % skip last ybuf points
xstride = 2;
ystride = 4;
tinit = 60;

% first read small amount of data
start = [isb-xbuf 1 1 tinit];
count = [xbuf*2/xstride (size(rgrid.h,1)-ybuf)/ystride Inf Inf];
stride = [xstride ystride 1 1];
xrange = start(1):start(1)+count(1)-1;
yrange = start(2):start(2)+count(2)-1;
u = squeeze(double(ncread(fname,'u',start - [1 0 0 0],count + [1 0 0 0],stride)));
v = squeeze(double(ncread(fname,'v',start - [0 0 0 0],count + [0 1 0 0],stride)));

% average to RHO points
u = avg1(u,1);
v = avg1(v,2);

u = permute(u,[2 1 3 4]);
v = permute(v,[2 1 3 4]);

% make grid matrices
xrmat = repmat(rgrid.x_rho(yrange,xrange)',[1 1 size(rgrid.z_r,1)]);
yrmat = repmat(rgrid.y_rho(yrange,xrange)',[1 1 size(rgrid.z_r,1)]);
zrmat = permute(rgrid.z_r,[2 3 1]);
zrmat = zrmat(yrange,xrange,:);
zrmat = permute(zrmat,[2 1 3]);


%% cone plot locations
cx = xrmat(1:10:end,1:10:end,1:10:end);
cy = yrmat(1:10:end,1:10:end,1:10:end);
cz = zrmat(1:10:end,1:10:end,1:10:end);

figure
for ii = 1:size(u,4)
    coneplot(xrmat,yrmat,zrmat,u(:,:,:,ii),v(:,:,:,ii),zeros(size(u(:,:,:,ii))), ...
             cx,cy,cz);
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
[C,hc] = contour(xr./1000,yr./1000,h,[200 500 750 1100],'k');
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

%% T/S census

R0 = ncread(fname,'R0');
TCOEF = ncread(fname,'Tcoef');
SCOEF = ncread(fname,'Scoef');

info = ncinfo(fname,'temp');

tinit  = double(ncread(fname,'temp',[1 1 1 1],[Inf Inf Inf 1]));
tfinal = double(ncread(fname,'temp',[1 1 1 info.Size(end)],[Inf Inf Inf 1]));
ri = tinit; R0 * (1 - TCOEF * tinit );
rf = tfinal; R0 * (1 - TCOEF * tfinal);

nbins = 40;
bins = linspace(min([ri(:);rf(:)]),max([ri(:);rf(:)]),nbins);
ni   = histc(ri(:),bins);
nf   = histc(rf(:),bins);

figure;
bar(bins,[ni nf],'histc');
legend('initial','final');

%% read ROMS floats output
flt = [dir '/ocean_flt.nc'];

floats = roms_read_floats(flt);
% 
% figure;
% %subplot(121)
% plot3(floats.x,floats.y,floats.z);
% hold on;
% plot3(floats.init(:,1)/1000,floats.init(:,2)/1000,floats.init(:,3),'x','MarkerSize',12);

%contourf(rgrid.x_rho,rgrid.y_rho,zeta(:,:,1)'); colorbar;
%surf(rgrid.x_rho/1000,rgrid.y_rho/1000,-rgrid.h);
% 
% subplot(122)
% hold on
% for ii = 1:size(floats.x,2) 
%     scatter3(floats.x(:,ii),floats.y(:,ii),floats.z(:,ii),16,floats.temp(:,ii))
% end
% colorbar
% caxis([19 20]);
% surf(rgrid.x_rho/1000,rgrid.y_rho/1000,-rgrid.h);

% animate floats

zeta = double(ncread(fname,'zeta'));

ntavg = size(zeta,3);
ntflt = size(floats.x,1) - 1;
floats.fac = ntflt/ntavg;

%% animate ROMS floats


    % need to color floats based on starting position
    colors = distinguishable_colors(size(floats.x,1));

    N = 35;
    n = 5;
    for i=1:n:N
       colarray(i:i+n-1,:,:) = repmat(colors(ceil(i/n),:,:),n,1);
    end

    cmap = flipud(cbrewer('div', 'BrBG', 32));

    figure;
    % plot bathymetry
    ax(2) = subplot(132);
    plot(rgrid.x_rho(1,:)/1000,-rgrid.h(1,:));
    liney(-114,'-114');
    liney(-130,'-130');
    beautify
    ax(3) = subplot(133);
    plot3(floats.x/1000,floats.y/1000,floats.z);
    hold on;
    plot3(floats.init(:,1)/1000,floats.init(:,2)/1000,floats.init(:,3),'x','MarkerSize',12);

    for i=31:ntavg
        ax(1) = subplot(131);
        cla
        contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,i)');
        shading flat; axis image
        colormap(cmap); 
        caxis([min(zeta(:)) max(zeta(:))]);
        hold on
        [C,hc] = contour(rgrid.x_rho./1000,rgrid.y_rho./1000,rgrid.h,[114 500 750 1100],'k');
        clabel(C,hc);
        plot(floats.init(:,1)/1000,floats.init(:,2)/1000,'x','MarkerSize',12);
        for j = 1:size(floats.x,2)
              plot(floats.x(1:floats.fac*i,j)/1000,floats.y(1:floats.fac*i,j)/1000, ...
                     'Color',colarray(j,:,:),'MarkerSize',10);
    %     scatter3(floats.x(1:fac*i,j),floats.y(1:fac*i,j),floats.z(1:fac*i,j), ...
    %                 10,floats.z(1:fac*i,j));
        end
        linkaxes(ax,'x')
        title(['t = ' num2str(i) ' days']);
        pause(0.01)
    end

%% read LTRANS floats

%% read TRACMASS floats
redo = 1;
if ~exist('trac','var') || redo == 1
    trac = tracmass_read('runs/eddytest_run.asc',rgrid);
end
if ~exist('zeta','var')
    zeta = ncread(fname,'zeta');
end

ntavg = size(zeta,3);
ntflt = size(trac.time,1) - 1;
trac.fac = ntflt/ntavg;

i = 1;
figure
hold on
plot(trac.x/1000,trac.y/1000);
plot(trac.x(1,i)/1000,trac.y(1,i)/1000,'x','MarkerSize',12);
title('tracmass');
plot3(trac.init(:,1)/1000,trac.init(:,2)/1000,trac.init(:,3),'x','MarkerSize',12);
linex(xsb/1000)
%animate_floats2d(rgrid,zeta,trac)

%% compare ltrans and ROMS output

ltrans = 'ltrans-eddy-04.nc';

% assume floats contains ROMS float data

[a b] = factor(size(floats.x,2));
figure
for i=1:a
    for j=1:b
        subplot(a,b,sub2ind([a b],i,j))
    end
end


%% 3D movie

az = 19.5;
el = 58;

shift = -800;
figure;
for i=1:ntavg
    clf;
    hzeta = surf(rgrid.x_rho/1000,rgrid.y_rho/1000,shift+zeta(:,:,i)');
    shading flat;
    caxis(shift+[min(zeta(:)) max(zeta(:))]);
    hold on
    hsurf = surf(rgrid.x_rho/1000,rgrid.y_rho/1000,-rgrid.h);
    set(hsurf,'FaceColor',[117 104 104]/255,'FaceAlpha',0.5,'EdgeAlpha',0);
    shading flat;
    for j = 1:size(floats.x,2)
          plot3(floats.x(1:fac*i,j),floats.y(1:fac*i,j),floats.z(1:fac*i,j), ...
                 'k.','MarkerSize',12);
%     scatter3(floats.x(1:fac*i,j),floats.y(1:fac*i,j),floats.z(1:fac*i,j), ...
%                 10,floats.z(1:fac*i,j));
    end
    view(az,el);
    pause(0.01)
end

