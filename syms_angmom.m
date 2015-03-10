function [slope, pbot, angmom] = syms_angmom(runs)

    ticstart = tic;
    syms r z
    syms Lr Lz theta positive
    syms f H bta positive
    syms ra negative
    syms vb

    bta0 = runs.params.phys.beta;
    g = runs.params.phys.g;
    r0 = runs.params.phys.rho0;

    %btrq = load([runs.dir '/btrqinterim.mat']);

    if ~isfield(runs.eddy, 'z0')
        if isempty(runs.usurf), runs.read_velsurf; end
        if isempty(runs.zeta), runs.read_zeta; end
        % fit gaussian to surface velocity
        tic;
        for tt=1:length(runs.eddy.cy)

            ix = find_approx(runs.rgrid.x_rho(1,:), runs.eddy.mx(tt));
            uvec = double(runs.usurf(ix,:,tt)');
            zvec = double(runs.zeta(ix,2:end-1,tt)');
            yuvec = runs.rgrid.y_u(:,ix) - runs.eddy.my(tt);
            yzvec = runs.rgrid.y_rho(2:end-1,ix) - ...
                    runs.eddy.my(tt);

            %ubvec = double(runs.ubot(ix,:,tt)');
            %ubsave(:,tt) = ubvec;
            %[Vb0(tt), runs.eddy.Lfitub(tt)] =
            %exp_fit(yuvec,ubvec,0);

            %pbvec = double(btrq.pbot(ix,:,tt)) .* runs.eddy.mask(ix,1: ...
            %end-3,tt);
                                                            %pbvec = fillnan(pbvec,0);
                                                            %pbvec = smooth(repnan(pbvec - nanmean(pbvec(:)),0),5);

            %pbsave(:,tt) = pbvec;

            ypvec = runs.rgrid.y_rho(2:end-1,ix) - ...
                    runs.eddy.my(tt);

            %[pbot0(tt), Lp(tt)] = exp_fit(ypvec(1:end-3), pbvec,1);

            %hold off; plot(yvec, uvec);
            %ylim([-1 1]* 0.2); xlim([-3 3]*1e5);
            %linex(0); liney(0);
            %pause(0.1);
            [V0(tt), runs.eddy.Lfit(tt)] = exp_fit(yuvec, uvec,0);
            [runs.eddy.z0(tt), runs.eddy.Lfitz(tt), z1(tt)] = gauss_fit(yzvec, ...
                                                              zvec-zvec(1), 0);
        end
        toc;

        eddy = runs.eddy;
        save([runs.dir '/eddytrack.mat'], 'eddy');
    end
    zmax = squeeze(max(max(runs.zeta,[],1),[],2) - min(min(runs.zeta, ...
                                                  [],1),[],2));
    zmax = squeeze(max(max(runs.zeta,[],1),[],2))';
    Lfitz = abs(runs.eddy.Lfitz);
    z0 = runs.eddy.z0;

    %xvec = xvec - runs.eddy.mx(1);
    %yvec = yvec - runs.eddy.my(1);
    %U = avg1(runs.bathy.h, 1).*ubar;
    %AM = trapz(yvec, trapz(avg1(xvec,1), bta0 * bsxfun(@times, U, yvec'), 1), 2);

    ra0 = -r0 * runs.params.phys.TCOEF * runs.eddy.T(:,end)';
    Lr0 = runs.eddy.Lfit; %runs.eddy.vor.dia/2;
    Lz0 = runs.eddy.Lgauss;
    H0 = runs.eddy.hcen';
    f0 = runs.eddy.fcen';
    %vb0 = runs.eddy.Vb;

    alpha = runs.bathy.sl_slope;

    %    bta0 = 1; ra0 = -1; Lr0 = 1; Lz0 = 1; H0 = Inf; vb0 = 0; f0 = ...
    %       1; g = 1; r0 = 1;

    % profile shapes
    R(r,Lr) = exp(-(r/Lr)^2);
    F(z,Lz) = exp(-(z/Lz)^2);

    rho(ra, r,Lr,z,Lz) = ra .* R(r,Lr) .* F(z,Lz);

    % geostrophic azimuthal velocity shear
    uz = -g/r0/f * diff(rho,r);
    ugeo(ra,r,Lr,z,Lz,f,H) = int(uz, z, -H, z);

    % agrees well with ubar
    %Ugeo = int(ugeo,z,-H,0);

    % 0.8L is where gaussian is maximum
    % this agrees with value of 'vgeo' used in roms_create
    usurf = double(ugeo(ra0,Lr0*0.8,Lr0,0,Lz0,f0,H0));

    vbr(vb,r,Lr) = vb * diff(R,r);

    amgeo(Lr,Lz,f,ra,H) = integrate_r(r.*sin(theta) .* ugeo.*-sin(theta));
    ambot(vb,Lr) = int(int(r*sin(theta)*vbr*-sin(theta), r, ...
                                  -Inf, Inf), theta, 0, 2*pi);

    am = double(amgeo(Lr0,Lz0,f0,ra0,H0));

    % % for some reason, it doesn't do the cancellation in pbc+pbt correctly.
    % %pbc(ra,Lr,Lz,H) = integrate_r(g/r0 * rho);
    % %pbt = int(g*zeta, r, -Inf, Inf);
    % %pbotr(ra,r,Lr,Lz,H) = g*(zeta + 1/r0 * int(rho, z,-H, 0));

    % % alternatively,
    % % read zeta field, fit exponential - get zeta amplitude;
    % % (eddy.amp turns out to be good enough)
    % % calculate zeta due to geostrophy - get amplitude = zampgeo.
    % % subtract - get zeta anomaly
    % % that's unbalanced barotropic pressure
    % zetageo(r,Lr,Lz,f,ra,H) =  ( - 1/r0 * int(rho,z,-H,0));
    % zampgeo = double(zetageo(0, Lr0,Lz0,f0,ra0,H0));

    % % ssh mask is definitely appropriate here
    % %zetamask = runs.zeta(2:end-1,2:end-1,:) .* runs.eddy.mask;
    % % define r and θ matrices
    % % [thmat,rmat] = cart2pol(bsxfun(@minus, runs.eddy.xr, permute(runs.eddy.mx,[3 ...
    % %                     1 2])), ...
    % %           bsxfun(@minus, runs.eddy.yr, permute(runs.eddy.my, [3 ...
    % %                     1 2])));
    % % zetageo0 = double(zetageo(rmat,Lr0,Lz0,f0,ra0,H0));

    % % geostrophically balanced bottom pressure
    % %pbottom(f,vb,Lr) = int(f * vb * Lr * R, r, -Inf, Inf);
    % syms pamp amp rend
    % pbottom(amp,Lr, rend) = int(g*amp*R,r,-Inf, rend);
    % pbottominf(amp,Lr) = int(g*amp*R,r,-Inf, Inf);
    % vbot(r,amp,Lr,f) = g/f*amp*diff(R,r);

    % %pbottom(pamp,Lr,rend) = int(pamp*diff(R,r), -Inf, rend);

    % % find how much of eddy is over slope
    % rmat = bsxfun(@minus, runs.rgrid.y_rho(:,1), runs.eddy.my);
    % rend0 = rmat(runs.bathy.isl,:);

    % % zeta anomaly
    % zanom = z0 - zampgeo;
    % zanom = zanom - zanom(3);

    % % bottom velocity magnitude
    % %vb1 = double(vbot(0.8*Lfitz,zanom,Lfitz,f0));

    % %pbot = double(pbc(ra0, Lr0, Lz0, H0) + ...
    % %              pbt(Lr0, Lz0, f0, vb0, ra0, H0));

    % trans(Lr,Lz,f,bta,ra,H) = integrate_r((f + bta*r*sin(theta))* ...
    %                                       zetageo)./H;

    % cfh = runs.eddy.cvy .* double(trans(Lr0,Lz0,f0,bta0,ra0,H0));
    % pbot =  double(pbottom(zanom, Lfitz, rend0));
    % %pbot = double(pbottom(pbot0, Lp, rend0));
    % pbotinf = double(pbottominf(zanom, Lfitz));

    % slope = fillnan(angmom./pbot, Inf);

    plots = 0;
    if plots
        [~,~,tind] = runs.locate_resistance;
        figure;
        subplot(211)
        plot(bta0*am);hold all
        plot(pbot*alpha);
        plot(cfh); ylim([-1 1]* 6)
        %plot(bta0*angmom./pbot); liney(alpha);
        linex([tind runs.traj.tind]);
        title(runs.name);
        legend('\beta * angmom', '\alpha * pbot', 'cfh');
        subplot(212)
        plot(runs.eddy.cvy);
        title('cross-isobath translation velocity');

        export_fig('-painters', ['images/angmom-' runs.name '.png']);
    end

    angmom = runs.angmom;
    angmom.sym_betatrq = bta0 .* am;
    angmom.hash = githash([mfilename('fullpath') '.m']);

    runs.angmom = angmom;

    save([runs.dir '/angmom.mat'], 'angmom');
end

function [out] = integrate_r(in)

    syms r z
    syms theta Lr Lx Ly Lz positive
    syms ra f
    syms H

    out = int(int(int(r*in,r,0,Inf),theta,0,2*pi),z,-H,0);
end

% I want to fit y = y_0 tanh(x/X)
function [y0, X] = exp_fit(x, y, plot_flag, test)

    if ~exist('test', 'var'), test = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if test
        test_fit();
        return;
    end

    initGuess(1) = max(y(:));
    initGuess(2) = max(x(:))/2;

    opts = optimset('MaxFunEvals',1e7);
    [fit2,~,exitflag] = fminsearch(@(fit) fiterror(fit,x,y), ...
                                   initGuess,opts);

    y0 = fit2(1);
    X = fit2(2);

    if plot_flag
        figure;
        plot(x,y,'k*'); hold all
        plot(x, y0*(2*x/X) .* exp( -(x/X).^2));
        liney([-1 1]* y0);
        linex([-1 1]*X);
    end
end

function [E] = fiterror(fit,x,y)
% x = (T0,H,a)
    y0 = fit(1); X = fit(2);

    E = sum((y - y0 .* (2*x/X) .* exp(-(x/X).^2)).^2);
end

function [] = test_fit()
    x = [-4:0.05:4];
    X = 1;
    y0 = 1;
    y = y0 * (2*x/X) .* exp(-(x/X).^2);

    [yy,xx] = exp_fit(x,y,1);

    disp(['y0 = ' num2str(yy) ' | Original = ' num2str(y0)]);
    disp(['X = ' num2str(xx) ' | Original = ' num2str(X)]);
end


% I want to fit y = y_0 tanh(x/X)
function [y0, X, y1] = gauss_fit(x, y, plot_flag, test)

    if ~exist('test', 'var'), test = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if test
        test_gauss_fit();
        return;
    end

    initGuess(1) = max(y(:));
    initGuess(2) = max(x(:))/2;
    initGuess(3) = max(y(:));

    opts = optimset('MaxFunEvals',1e7);
    [fit2,~,exitflag] = fminsearch(@(fit) gaussfiterror(fit,x,y), ...
                                   initGuess,opts);

    y0 = fit2(1);
    X = abs(fit2(2));
    y1 = fit2(3);

    if plot_flag
        figure;
        plot(x,y,'k*'); hold all
        plot(x, y0 * exp( -(x/X).^2) + y1*(x/X));
        liney([-1 1]* y0);
        linex([-1 1]*X);
    end
end

function [E] = gaussfiterror(fit,x,y)
% x = (T0,H,a)
    y0 = fit(1); X = fit(2); y1 = fit(3);

    E = sum((y - y0 .* exp(-(x/X).^2) - y1 * (x/X)).^2);
end

function [] = test_gauss_fit()
    x = [-4:0.05:4];
    X = 1;
    y0 = 1;
    y1 = 1;

    y = y0 * exp(-(x/X).^2) + y1 * (x/X);

    [yy,xx,yy1] = gauss_fit(x,y,1);

    disp(['y0 = ' num2str(yy) ' | Original = ' num2str(y0)]);
    disp(['y1 = ' num2str(yy1) ' | Original = ' num2str(y1)]);
    disp(['X = ' num2str(xx) ' | Original = ' num2str(X)]);
end
