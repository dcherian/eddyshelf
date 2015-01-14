% use expressions to determine energy

function [intTE] = eddy_energy_ideal(runs)

    % physical parameters
    TCOEF = runs.params.phys.T0;
    rho0 = runs.params.phys.R0;
    g  = runs.params.phys.g;
    f0 = runs.params.phys.f0;

    % extract surface temperature anomaly
    Tamp = runs.eddy.T(:,end);

    Lx = runs.eddy.vor.lmaj;
    Ly = runs.eddy.vor.lmin;
    Lz = runs.eddy.Lgauss;

    syms x y z
    syms lx ly lz Ta positive

    % temperature anomaly
    T(Ta,lx,ly,lz,x,y,z) = Ta * exp(- (x/lx)^2 - (y/ly)^2 - (z/lz)^2);
    intT = domain_integrate(T);

    % density anomaly
    rho = -rho0 * TCOEF * T;

    % Potential energy
    PE = int(rho * g * z, z);
    intPE = domain_integrate(PE);

    % kinetic energy
    uz = -TCOEF*g/f0 * diff(T,y);
    vz = TCOEF*g/f0 * diff(T,x);

    u = int(uz, z, -Inf, z);
    v = int(vz, z, -Inf, z);

    KE = u^2 + v^2;
    intKE = domain_integrate(KE);

    tic;
    intTE = double(intKE(Tamp', Lx, Ly, Lz)) + double(intPE(Tamp', Lx, Ly, Lz));
    toc;

    % volume
    vol(lx,ly,lz) = intT / Ta;
    vol0 = vol(Lx, Ly, Lz);

    % get group velocity
    cg = runs.topowaves;

    % plots
    figure;
    plot(intTE./max(intTE));
    hold all
    plot(vol0./max(vol0));
    legend('Total Energy', 'Volume');
    title(runs.name);
    beautify([18 20 22]);

    % check group velocity
    % ∂E/∂t + cg.∇E = 0
    %figure;

    runs.eddy.intTE = intTE;
end

function [intf] = domain_integrate(f)

    syms x y z
    intf = int(int(int(f, x, -Inf, Inf), y, -Inf, Inf),z,-Inf,0);
end