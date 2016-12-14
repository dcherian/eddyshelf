% elliptical eddy structure

syms x y Lx Ly V0 R0 u0 v0

r = sqrt(x^2 + y^2);

R = R0 * exp(-(x/Lx)^2 - (y/Ly)^2);
u = -diff(R,y);
v = diff(R,x);

ux = diff(u,x);
vy = diff(v,y);

ux + vy

%% sds

xvec = [-3:0.05:3]';
zvec = linspace(-2,0,60);

R = 0.7; % radius to max vel. at surface
Rcyc = 1; % radius to cyclonic "balloon"
hsb = -0.15; % really hsb/Lz
L_fil = R/3; % filament length
a = 2;
videal = abs(bsxfun(@times, -xvec.^(a-1) .* exp(-abs(xvec).^a), 1-erf(-zvec)));

% hypotenuse is z = mx + c
m = hsb/(L_fil);
c = m * R;

mask = repmat(xvec < -Rcyc, [1 length(zvec)]);
trimask = bsxfun(@and, bsxfun(@plus, zvec,  - (m*xvec + c)) < 0, ...
                 zvec > hsb);
shmask = mask | trimask;
vmasked = (videal .* shmask);

figure;
subplot(121);
pcolorcen(xvec, zvec, vmasked');
hold on
contour(xvec, zvec, videal', 40, 'k');
center_colorbar;

subplot(122);
plot(trapz(xvec, vmasked, 1), zvec);
liney(hsb);