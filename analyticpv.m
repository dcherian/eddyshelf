syms x y
syms z negative
syms R V0 Lz positive
syms rho0 f

v = V0 * (x) * exp(-(x^2 + y^2)) * (1 - erf(-z));
u = -V0 * (y) * exp(-(x^2 + y^2)) * (1 - erf(-z));
rho = rho0 * exp(-(x^2+y^2 + z^2));

pv(x,y,z,V0,rho0) = (f + diff(v,x) - diff(u,y)) * diff(rho, z) ...
     - diff(v,z)*diff(rho,x) + diff(u,z)*diff(rho,y)