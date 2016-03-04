% some quick and dirty estimates for the MAB

woa = load('data/woa05.mat');
load('data/etopo2_extract.mat');

% isobath limits for S_α, β_t distributions
hsb = -250;
hend = -80;

% MAB
latdeep =  39;
latshelf = 41;

londeep = -70;
lonshelf = -70;


% south of NY
%londeep = -71.9;
%lonshelf = -73.9;

%latdeep  = 37.5;
%latshelf = 39;

% george's bank
%latdeep = 40;
%londeep = -62.5;
%latshelf = 41.5;
%lonshelf = -67.5;

% george's bank 2
latdeep = 38.85;
latshelf = 40.7;

londeep = -68.13;
lonshelf = -69.14;

f0 = 2 * (2*pi/86400) * sind(latdeep);

% deep water point
ixdtopo = find_approx(topo.x, londeep);
iydtopo = find_approx(topo.y, latdeep);
ixdwoa = find_approx(woa.X, 360 + londeep);
iydwoa = find_approx(woa.Y, latdeep);

% shelf point
ixstopo = find_approx(topo.x, lonshelf);
iystopo = find_approx(topo.y, latshelf);
ixswoa = find_approx(woa.X, 360 + lonshelf);
iyswoa = find_approx(woa.Y, latshelf);

% figure out indices for bathymetry / N² cross-section
[ixt, iyt] = bresenham(ixdtopo, iydtopo, ixstopo, iystopo);
[ixw, iyw] = bresenham(ixdwoa, iydwoa, ixswoa, iyswoa);

zind = sub2ind(size(topo.z), ixt, iyt);

%% distance for gradient calculation
nsmooth = 3;
dist = sw_dist(topo.y(iyt), topo.x(ixt), 'km')*1000; % Δx
xvec = [0;cumsum(dist)]/1000;
h = topo.z(zind);
dhdx = diff(smooth(h,nsmooth))./dist;
dh2dx2(2:length(dhdx)) = diff(smooth(h,nsmooth), 2)./avg1(dist).^2;
dh2dx2(1) = 0; dh2dx2(end+1) = 0;
isb = find_approx(h, hsb, 1);
iend = find_approx(h, hend, 1);

% bottom slope
alpha = max(abs(dhdx));

%%
% make sure I'm in deep enough water
%topo.z(ixtopo, iytopo)

% 1st index is deep water, last is shelf water
for ii=1:length(ixw)
    N2mat(ii,:) = bfrq(woa.sal(:,iyw(ii),ixw(ii)), ...
                       woa.temp(:,iyw(ii), ixw(ii)), ...
                       woa.Z, latdeep);
end

[Vmode, Hmode, c] = vertmode(N2mat(1,:)', woa.Z, 1, 0);
%subplot(1,3,3);
%liney(abs(zbc1))
znew = avg1(woa.Z);
zbc1 = znew(find_approx(Hmode,0,1));

Ndeep = sqrt(max(N2mat(1,:)));

% now for slope - find N² at location of max dh/dx
[~,imax] = max(dhdx);
% convert indices to woa grid
ixslw = find_approx(woa.X, 360+topo.x(ixt(imax)))-1;
iyslw = find_approx(woa.Y, topo.y(iyt(imax)))-1;
% calculate N2
N2slope = bfrq(woa.sal(:,iyslw, ixslw), woa.temp(:,iyslw, ixslw), ...
               woa.Z, woa.Y(iyslw));
Nslope = sqrt(max(N2slope));

%% shelf
Nmax = sqrt(max(N2mat, [], 2)) % shouldn't be too different
S_sh = abs(dhdx(isb+1:iend) * mean(Nmax) ./ f0);
beta = abs(f0./avg1(h(isb-1:iend)) .* dhdx(isb:iend));

%% find zero crossing of horizontal velocity mode
disp(['Total depth = ' num2str(topo.z(ixdtopo, iydtopo)) ' m | ' ...
      'Zero crossing of BC1 at approx. ' ...
       num2str(zbc1), ' m']);

fprintf('\n Max N2 = %.3e s^{-2} | N/f = %.2f | L_D = %.3f km \n', ...
         max(N2mat(1,:)), Ndeep/f0, c/f0/1000)
fprintf('\n Burger_slope = %.4f | alpha = %.4f \n', Ndeep/f0*alpha, alpha)
fprintf('\n Shelfbreak depth = 130 m | Shelfbreak depth / zero crossing depth = %.4f\n ', ...
         130/zbc1)
fprintf('\n Slope width = 175 km | Slope width / L_D = %.2f \n', 175000/c*f0)
fprintf('\n But above Lsl/Ld scaling doesn''t work for constant slope.\n');
fprintf('for constant max.(dhdx) slope, Lsl = 80km | Lsl/Ld = %.2f \n', 80000/c*f0);

%% plot
figure; maximize; insertAnnotation('MAB.m');
ax(1) = subplot(2,2,1);
contour(topo.x, topo.y, topo.z', [-100 -200 -500 -1000 -2000 -3000 -4000], 'k');
hold on
contour(topo.x, topo.y, topo.z', [0 0], 'k', 'LineWidth', 3);
plot(londeep, latdeep, 'r*', lonshelf, latshelf, 'b.', 'MarkerSize', 16);
plot(topo.x(ixt), topo.y(iyt), 'r-');
beautify; axis image
xlim([-76 -60]);
ylim([30 45]);
title(['(' num2str(londeep) ',' num2str(latdeep) ') to (' ...
       num2str(lonshelf) ',' num2str(latshelf) ') | ' ...
      num2str(hend) 'm to ' num2str(hsb) 'm']);

ax(2) = subplot(2,2,2);
plot(xvec, h,'k-'); hold on
plot(xvec, smooth(h,nsmooth) ,'r*-')
linex(xvec(isb));
%plot(avg1(xvec), dhdx * 1e4, 'b-');
%liney(alpha * -1e4);
legend('raw h', 'smoothed h');
ax(2).YAxis.Exponent = 0;
beautify;

ax(3) = subplot(2,2,3);
hist(S_sh(:),[0:0.05:1.5]);
ax(3).XTick = [0:0.1:1.5];
ax(3).XAxis.TickLabelRotation = 40;
ax(3).XLim = [-0.1 1.5];
title('Distribution of \alpha N/f');
beautify;

ax(4) = subplot(2,2,4);
hist(beta,[0:1:10]*1e-9);
xlim([0 10]*1e-9);
title('Distribution of \beta = f0/h dh/dx');
beautify;