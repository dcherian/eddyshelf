% numerical version of syms_angmom

% Take scales
% Create gaussian at location (over topography if necessary)
% do the integrations numerically.

function [pbot,angmom] = num_angmom(runs)

    mx = runs.eddy.mx;
    my = runs.eddy.my;

    Lr = runs.eddy.Lfit;
    Lz = runs.eddy.Lgauss;
    f0 = runs.eddy.fcen;

    g = runs.params.phys.g;
    beta = runs.params.phys.beta;
    rho0 = runs.params.phys.rho0;

    rhoamp = - rho0 * runs.params.phys.TCOEF * runs.eddy.T(:,end);
    %if isempty(runs.zeta), runs.read_zeta; end;

    % ubar = dc_roms_read_data(runs.dir, 'ubar', [], {}, [], runs.rgrid);

    % number of radii
    nr = 4;
    % number of timesteps
    nt = length(Lz);
    N = runs.rgrid.N;
    [Lm Mm] = size(runs.rgrid.xr);

    xvec = runs.rgrid.xr(:,1);
    yvec = runs.rgrid.yr(1,:);

    ixm = vecfind(xvec,mx);
    iym = vecfind(yvec,my);
    ixmin = find_approx(xvec, mx(1)-nr*Lr(1)) - ixm(1);
    ixmax = find_approx(xvec, mx(1)+nr*Lr(1)) - ixm(1);
    iymin = find_approx(yvec, my(1)-nr*Lr(1)) - iym(1);
    iymax = find_approx(yvec, my(1)+nr*Lr(1)) - iym(1);

    dix = max(abs([ixmin ixmax]));
    diy = max(abs([iymin iymax]));

    dzmat0 = diff(permute(runs.rgrid.z_w,[3 2 1]), 1, 3);

    betatrq = nan(size(Lz));
    pbot = nan(size(Lz));

    tic;
    for tt=1:nt
        disp(['tt = ' num2str(tt) '/' num2str(nt)]);
        xx = [-1 1]*dix + ixm(tt);
        yy = [-1 1]*diy + iym(tt);

        xmin = max(xx(1),1); xmax = min(xx(2), Lm-2);
        ymin = max(yy(1),1); ymax = min(yy(2), Mm-2);

        % (x,y,1,t)
        xmat = runs.rgrid.xr(xmin:xmax,ymin:ymax) - mx(tt);
        ymat = runs.rgrid.yr(xmin:xmax,ymin:ymax) - my(tt);

        % center in reference frame
        icx = ixm(tt) - xmin + 1;
        icy = iym(tt) - ymin + 1;

        % vectors for integration
        xv = xmat(:,1);
        yv = ymat(1,:);

        % to polar co-ordinates!
        [thmat,rmat] = cart2pol(xmat, ymat);
        zrmat = permute(runs.rgrid.z_r,[3 2 1]);

        zrmat = zrmat(xmin:xmax,ymin:ymax,:);

        % water depth
        hmat = runs.rgrid.h(ymin:ymax,xmin:xmax)';
        slbot = diff(hmat,1,2)./diff(runs.rgrid.yr(xmin:xmax,ymin:ymax)- ...
                                     my(tt) ,1,2);
        slbot(:,end+1) = slbot(:,end);

        % dz
        dzmat = dzmat0(xmin:xmax, ymin:ymax, :);

        % for plots
        yz = repmat(yv',[1 72]);
        zy = squeeze(zrmat(icx,:,:));

        % actual calculations

        % density
        rho = rhoamp(tt) .* bsxfun(@times, exp(- (rmat/Lr(tt)).^2), exp(- ...
                                                          (zrmat/Lz(tt)).^2));

        % drho / dr
        drhodr = bsxfun(@times, rho, -2* rmat./Lr(tt)/Lr(tt));
        % thermal wind shear
        vgeo = cumsum(dzmat .* (-9.81./rho0./f0(tt) * drhodr), 3);
        % gradient wind
        rfb2 = rmat.*f0(tt) ./ 2;
        sdisc = sqrt(1 + bsxfun(@times,vgeo,2./rfb2));% sqrt(discriminant)
        if isreal(sdisc) % gradient wind doesn't always work with anticyclones
            rut = bsxfun(@times,(-1 + sdisc), rfb2);
        else
            error(['gradient wind calculated complex v! - ' ...
                   'Ro > 0.25']);
        end
        u = -1 * bsxfun(@times, rut, sin(thmat));

        % angular momentum
        betatrq(tt) = integrate(xv, yv, beta .* ...
                               bsxfun(@times, sum(u.*dzmat, 3), ...
                                      yv));

        % sshmask = runs.eddy.mask(xmin:xmax, ymin:ymax, tt);

        % U = ubar(xmin:xmax, ymin:ymax, tt);
        % Ugeo = 1./hmat .* sum(dzmat .* u, 3);
        % ubot{tt} =  U-Ugeo;
        % ubamp(tt) = max(max(abs(ubot{tt} .* sshmask)));

        % pbamp(tt) = max(max( cumtrapz(yv, ubot{tt}, 2),[], 1), [], 2) ...
        %     .* rho0 .* f0(tt) / rho0; % divide by rho0 to agree
        %                               % with Flierl (1987)

        % pbottom{tt} = pbamp(tt) * exp(-(rmat/Lr(tt)).^2);
        % test_ub = 1./f0(tt)./rho0 .* bsxfun(@rdivide, diff(pbottom{tt},1,2), diff(yv));

        % % free-surface
        % zetageo = sum(dzmat .* (-rho./rho0), 3);
        % zamp(tt) = max(zetageo(:));
        % zeta = runs.zeta(xmin:xmax, ymin:ymax, tt);

        % pbotgeo = integrate(xv, yv, g * zetageo .* slbot);
        % pbotzetafull = integrate(xv, yv, g*zeta .*slbot);
        % pbotzetamask = integrate(xv, yv, g*zeta.*sshmask .* slbot);
        % pbotamp(tt) = integrate(xv, yv, pbottom{tt}.*slbot);
        % pbot(tt) = pbotzetafull-pbotgeo;
        % pbotmask(tt) = pbotzetamask-pbotgeo;
    end
    toc;

    angmom.betatrq = betatrq;
    angmom.hash = githash([mfilename('fullpath') '.m']);

    save([runs.dir '/angmom.mat'], 'angmom');
    keyboard;
    % save to file
end

function [out] = integrate(xvec, yvec, in)
    out = trapz(xvec, trapz(yvec, in, 2), 1);
end