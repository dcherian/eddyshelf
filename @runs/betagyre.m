function [] = betagyre(runs)

    runs.read_zeta;

    tic;
    for tt=1:size(runs.zeta,3)
        ix = find_approx(runs.rgrid.x_rho(1,:), runs.eddy.mx(tt));
        iy = find_approx(runs.rgrid.y_rho(:,1), runs.eddy.my(tt));

        amp = runs.zeta(ix,iy,tt);

        % figure out length scales
        zxvec = double(runs.zeta(ix,2:end-1,tt));
        zyvec = double(runs.zeta(2:end-1,iy,tt)');
        yzvec = runs.rgrid.y_rho(2:end-1,ix)' - ...
                runs.eddy.my(tt);
        xzvec = runs.rgrid.x_rho(iy,2:end-1) - ...
                runs.eddy.mx(tt);
        [z0x(tt), Lxfit(tt)] = gauss_fit(xzvec, zyvec-zyvec(1), 0);
        [z0y(tt), Lyfit(tt)] = gauss_fit(yzvec, zxvec-zxvec(1), 0);

        r0 = hypot(Lxfit(tt), Lyfit(tt));
        [~,rmat] = cart2pol(runs.rgrid.x_rho' - runs.eddy.mx(tt), ...
                            runs.rgrid.y_rho' - runs.eddy.my(tt));

        % symmetric part = Gaussian
        zetasym = mean([z0x(tt) z0y(tt)]) .* ...
                  exp(-(rmat./r0).^runs.params.eddy.a);
        runs.eddy.betagyre(:,:,tt) = runs.zeta(:,:,tt) - zetasym;
    end
    toc;
end