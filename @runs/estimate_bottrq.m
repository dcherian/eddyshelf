% Use symbolic math toolbox to estimate time series of angular
% momentum balance terms

function [angmom,pbot] = estimate_bottrq(runs)

    syms x y
    syms z negative
    syms Lx Ly Lz amp H positive
    syms ra

    % physical params
    beta = runs.params.phys.beta;
    f0 = abs(runs.params.phys.f0);
    r0 = runs.params.phys.rho0;
    g = runs.params.phys.g;
    alpha = runs.bathy.sl_slope;

    % fit
    [~,~,tind] = runs.locate_resistance;

    runs.read_zeta;
    runs.read_rhosurf(1,1);
    for tt = [1, tind]
        ix = find_approx(runs.rgrid.x_rho(1,:), runs.eddy.mx(tt));
        iy = find_approx(runs.eddy.my(tt), runs.rgrid.y_rho(:,1));

        zxvec = double(runs.zeta(ix,2:end-1,tt));
        zyvec = double(runs.zeta(2:end-1,iy,tt)');
        yzvec = runs.rgrid.y_rho(2:end-1,ix)' - ...
                runs.eddy.my(tt);
        xzvec = runs.rgrid.x_rho(iy,2:end-1) - ...
                runs.eddy.mx(tt);
        [z0x(tt), Lxfit(tt)] = gauss_fit(xzvec, zyvec-zyvec(1), 0);
        [z0y(tt), Lyfit(tt)] = gauss_fit(yzvec, zxvec-zxvec(1), 0);
    end

    % eddy params
    ra0 = runs.eddy.drhothresh(1);
    Lx0 = Lxfit(1);
    Ly0 = Lyfit(1);
    Lz0 = runs.eddy.Lgauss(1);
    Lxt = Lxfit(tind);
    Lyt = Lyfit(tind);
    Lzt = runs.eddy.Lgauss(tind);
    Ht = runs.eddy.hcen(tind);
    ampt = z0x(tind);

    % fields
    zeta = amp * exp(-(x/Lx)^2 - (y/Ly)^2);
    rho = ra * exp(-(x/Lx)^2 - (y/Ly)^2 - (z/Lz)^2);
    uz = -g/r0/f0 * diff(rho,y);
    u = int(uz, z, -Inf, z);

    am(ra,Lx,Ly,Lz,H) = beta * integrate(y*u);
    IZETA = int(int(zeta ,x,-Inf,Inf), y,-Inf,Inf);
    pb(ra,Lx,Ly,Lz,H,amp) = -alpha * g* (IZETA +  1/r0 * integrate(rho));

    % initial angular momentum
    angmom = double(am(ra0,Lx0,Ly0,Lz0,10000*Lz0));
    pbot = double(pb(ra0,Lxt,Lyt,Lzt,Ht,ampt));
end

function [out] = integrate(in)

    syms x y z
    syms Lx Ly Lz H positive
    syms ra g r0 f
    syms bta positive

    out = int(int(int(in,x,-Inf,Inf),y,-Inf,Inf),z,-H,0);
end