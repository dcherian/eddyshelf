% make idealized circular kink
% integrate that to see what the x and z profiles look like.

x = linspace(-3,3,100);
z = linspace(-1.5,0,200);

[xmat,zmat] = ndgrid(x,z);

v = -xmat .* exp(-xmat.^2) .* (1-erf(-zmat));

x0 = -1;
z0 = -0.15;
kzrad = z0; % kink radius
kxrad = 0.3;
% 0.3 is a good value because max. vel is at 0.7, and it looks like that's where the streamer
% cuts in.

mask = (xmat < 0) & ( ((xmat.^2 + zmat.^2) > 1) ... % eddy shape
       | ( (((xmat-x0)/kxrad).^2 + ((zmat-z0)/kzrad).^2) <= 1 )); % kink

figure;
subplot(3,3,[4 5 7 8]);
pcolorcen(xmat, zmat, v);
clim = caxis;
hold on
contour(xmat, zmat, mask, [1 1], 'k', 'LineWidth', 2);
caxis(clim); center_colorbar;
subplot(3,3,[1 2]);
plot(x, trapz(z,v.*mask, 2));
subplot(3,3,[6 9])
plot(trapz(x,v.*mask, 1), z);
