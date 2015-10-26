function [] = plot_simplepv(runs)
% this function contours the qgpv approximation of the
% background pv

    if runs.bathy.axis == 'y'
        dhdx = diff(runs.bathy.h,1,2)./diff(runs.rgrid.yr,1,2);
        ax = 2;
    else
        ax = 1;
        dhdx = diff(runs.bathy.h,1,1)./diff(runs.rgrid.xr,1,1);
    end

    beta_t = runs.params.phys.f0 * dhdx/max(runs.rgrid.zr(:));

    q = runs.params.phys.f0 + ...
        (runs.params.phys.beta + beta_t) .* avg1(runs.rgrid.yr,ax);

    clf;
    subplot(211);
    contourf(q');
    subplot(212);
    hold on
    plot(q(2,:));
    plot(-runs.bathy.h(2,:)/max(runs.bathy.h(:)),'k');
    legend('qgpv','bathy');

end
