%% integral (1 - erf(-z/Lz))

%syms z 
%syms H Lz positive

%I(H,Lz) = int( 1 - erf(-z/Lz), z, -H, 0);

HoLz = [0:0.01:10];

f = HoLz + 1/sqrt(pi) - HoLz .* erf(HoLz) - 1/sqrt(pi) .* exp(-HoLz.^2);
% approx = HoLz - HoLz.^2 / sqrt(pi) - HoLz.^4/6/sqrt(pi)
approx = 1/sqrt(pi) + 2*HoLz - 2/sqrt(pi) * exp(-HoLz.^2);
plot(HoLz, f)
hold on
plot(HoLz, approx)

%%
dir = 'runs/runte-11/';
varname = 'temp';

for i=1:11
    fname = [dir sprintf('ocean_his_00%02d',i) '.nc'];
    if i == 1
        v = ncread(fname,varname);
        rgrid = roms_get_grid(fname,fname,0);
    else
        temp = ncread(fname,varname);
        v(:,:,:,end+1:end+size(temp,4)) = temp;
    end
end
%%
%for ii=1:size(v,4)
%    temp = v(:,:,:,ii);
    metric= domain_integrate(v,rgrid.x_rho',rgrid.y_rho',permute(rgrid.z_r,[3 2 1]));
%end

plot(metric);

%%

beta = 2e-11;
bg.shear = -beta;
yvec = yrmat(1,:,1);
Y = max(yrmat(:));

ubt0 = 0.05;
ubtmax = 0.05;

% u = beta * y^2/2 + C1 * y + C2
% u = 0 @ y = 0 => C2 = 0
% u = ubtmax @ y = Y => 
C1 = ubtmax/Y - bg.shear * Y/2;

ubt = bg.shear * yvec.^2/2 + C1 * yvec;

plot(yvec/1000,ubt)

%% internal wave dispersion relation as modified by discretization

vector = linspace(0,4*pi,30);
[kd,ld] = meshgrid(vector,vector);

RRd = 40;

exact = sqrt(1 + ...
        (RRd)^2*(kd.^2 + ld.^2));
cgrid = sqrt(cos(kd/2).^2 .* cos(ld/2).^2 + ...
         4*(RRd)^2*(sin(kd/2).^2 + sin(ld/2).^2));

subplot(131)
pcolorcen(kd/pi,ld/pi,exact);
title('exact')
colorbar;axis square;
ylabel('l\Delta x/\pi');
xlabel('k\Delta x/\pi');
clim = caxis;
subplot(132)
pcolorcen(kd/pi,ld/pi,cgrid);
title('cgrid');
xlabel('k\Delta x/\pi');
colorbar;axis square; caxis(clim);
subplot(133)
pcolorcen(kd/pi,ld/pi,cgrid./exact);
title('cgrid/exact');
xlabel('k\Delta x/\pi');
colorbar;axis square
maximize; pause(0.2);
spaceplots(0.06*ones([1 4]),0.05*ones([1 2]))

%% float code

figure
for ii=1:size(ltrans.x, 1)
    clf;
    plot(ltrans.x(ii,:)/1000, ltrans.y(ii,:)/1000, '*');
    title(num2str(ii));
    pause(0.05);
end