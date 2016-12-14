syms x y
syms z negative
syms Lx Ly Lz H positive

% profile shapes
R(x,y,Lx,Ly) = exp(-(x/Lx).^2 - (y./Ly).^2);
F(z,Lz) = exp(-(z/Lz).^2);

rho = R(x,y,Lx,Ly) .* F(z,Lz);

% geostrophic azimuthal velocity shear
uz(x,y,Lx,Ly,z,Lz) = diff(rho,y);
ugeo(x,y,Lx,Ly,z,Lz,H) = int(uz, z, -H, z);

AM = int(int(int(y.*ugeo,x,-Inf,Inf),y,-Inf,Inf),z, -H, 0)

syms Us

U = Us * (y/Ly) * exp(-(x/Lx).^2 - (y./Ly).^2) * (1-erf(-z/Lz))
int(int(int(y*U,x,-Inf,Inf),y,-Inf,Inf),z,-H,0)