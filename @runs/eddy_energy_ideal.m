% use expressions to determine energy
function [intTE] = eddy_energy_ideal(runs)

    % physical parameters
    TCOEF = runs.params.phys.T0;
    rho0 = runs.params.phys.R0;
    g  = runs.params.phys.g;
    f0 = runs.params.phys.f0;
    beta = runs.params.phys.beta;
    Hsb = runs.bathy.hsb;
    alpha = runs.bathy.sl_slope;

    % extract surface temperature anomaly
    Tamp = runs.eddy.T(:,end)';

    Lx = runs.eddy.vor.lmaj;
    Ly = runs.eddy.vor.lmin;
    Lz = runs.eddy.Lgauss;

    % origin is at eddy center
    % H0 is depth at eddy center
    syms x y z H0
    syms lx ly lz Ta positive
    H = symfun(sym(['H0 + ' num2str(alpha) '*y']), [H0 y])

    % temperature anomaly
    T(Ta,lx,ly,lz,x,y,z,H0) = Ta * exp(- (x/lx)^2 - (y/ly)^2 - (z/lz)^2);
    intT(Ta,lx,ly,lz,H0) = domain_integrate(T, H);

    % density anomaly
    rho = -rho0 * TCOEF * T;

    % Potential energy
    PE = int(rho * g * z, z);
    intPE(Ta,lx,ly,lz) = domain_integrate(PE, H);

    % kinetic energy
    uz = -TCOEF*g/f0 * diff(T,y);
    vz =  TCOEF*g/f0 * diff(T,x);

    u = int(uz, z, -H, z);
    v = int(vz, z, -H, z);

    KE = u^2 + v^2;
    intKE(Ta,lx,ly,lz) = domain_integrate(KE, H);

    err = Inf;
    H0 = 100;
    while err > 1e-3
        H0 = H0 + 1;
        Hh = H(H0, y);

        u0(x,y) = u(Tamp(1), Lx(1), Ly(1), Lz(1), x, y, H0);
        rho1(x,y) = rho(Tamp(1), Lx(1), Ly(1), Lz(1), x, y, z, H0);
        % beta torque
        Fbeta(y) = beta .* domain_integrate(y*u0, Hh);

        % bottom pressure torque
        Fbot(y) = domain_integrate(int(rho1*g./rho0, z, -Hh, 0), Hh);

        Fbeta0 = vpa(Fbeta(Tamp(1), Lx(1), Ly(1), Lz(1), H0))
        Fbot0 = vpa(Fbot(Tamp(1), Lx(1), Ly(1), Lz(1), H0))

        err = Fbeta0 - Fbot0
    end

    tic;
    intTE = double(intKE(Tamp, Lx, Ly, Lz, H0)) + double(intPE(Tamp, ...
                                                      Lx, Ly, Lz, H0));
    toc;

    % volume
    vol(lx,ly,lz,H) = intT / Ta;
    vol0 = vol(Lx, Ly, Lz, H0);

    % get group velocity
    cg = runs.topowaves;

    % plots
    %figure;
    %plot(intTE./max(intTE));
    %hold all
    %plot(vol0./max(vol0));
    %legend('Total Energy', 'Volume');
    %title(runs.name);
    %beautify([18 20 22]);

    % check group velocity
    % ∂E/∂t + cg.∇E = 0
    %figure;

    runs.eddy.energy.vol = vol0;
    runs.eddy.energy.intTE = intTE;
    runs.eddy.energy.intKE = intKE;
    runs.eddy.energy.intPE = intPE;
    runs.eddy.energy.hash = githash([mfilename('fullpath') '.m']);;

    energy = runs.eddy.energy;
    save([runs.dir '/energy.mat'], 'energy');
end

function [intf] = domain_integrate(f, H)

    syms x y z
    intf = int(int(int(f, x, -Inf, Inf), y, -Inf, Inf),z,-H,0);
end