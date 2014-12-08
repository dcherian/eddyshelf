run = ew04;

rho = dc_roms_read_data(run.dir, 'rho', [1 1], {}, [], run.rgrid);

ix = find_approx(run.rgrid.x_rho(1,:), run.eddy.mx(1));
iy = find_approx(run.rgrid.y_rho(:,1), run.eddy.my(1));
zvec = squeeze(run.rgrid.z_r(:,iy,ix));
rback = squeeze(rho(end,end,:));
redd = squeeze(rho(ix,iy,:));

N2back = -9.81/1025 * diff(rback)./diff(zvec);
N2edd = -9.81/1025 * diff(redd)./diff(zvec);

figure
hold all
plot(N2back);
plot(N2edd);

Ro = run.eddy.Ro(1)/2;
f = run.params.phys.f0;
alpha = Ro*(1+Ro) * f^2 ./ (N2edd(end) - N2back(end))
Bu = Ro*(1+Ro) ./ (N2edd(end)./N2back(end) - 1)
H = sqrt(alpha) * run.params.eddy.Ldef