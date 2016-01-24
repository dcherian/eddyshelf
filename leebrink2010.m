salname = '../leebrink2010/oa-field-salinity.mat';
velname = '../leebrink2010/oa-field-velocity.mat';

S = load(salname);
V = load(velname);

figure;
insertAnnotation('leebrink2010.m');
hax(1) = subplot(1,3,[1 2]);
contourf(V.xgv, V.zgv, V.vgridded, 80, 'EdgeColor', 'none');
hold on
[cc,hh] = contour(S.xgsal, S.zgsal, S.salgridded, [33:0.5:35], 'k');
clabel(cc,hh);
caxis([-1 1]*0.5); hcb = center_colorbar;
hcb.Label.String = '';
set(gcf, 'renderer', 'painters');
ylabel('Z (m)');
xlabel('X (km)');
xlim([0 40]);
title('Lee & Brink(2010) | Negative velocity = offshore');
beautify;

% do the integration
% 1) interpolate velocity onto salinity grid
Vinterp = interp2(V.xgv, V.zgv', V.vgridded, S.xgsal, S.zgsal');
Sinterp = interp2(S.xgsal, S.zgsal', S.salgridded, V.xgv, V.zgv');

hax(2) = subplot(1,3,3); hold on;
lightDarkLines(5);
for S0=33:0.5:35

    Tprofile = trapz(S.xgsal*1e3, ...
                     repnan(Vinterp,0) .* (repnan(S.salgridded,1200) < S0) ...
                     .* (repnan(Vinterp,0) < 0), 2);
    axes(hax(2))
    plot(Tprofile, S.zgsal);
    Tvi = trapz(S.zgsal, Tprofile, 1);

    Tsi = trapz(V.xgv*1e3, trapz(V.zgv, ...
                                 repnan(V.vgridded,0) .* (repnan(Sinterp,1200) < S0) ...
                                 .* (repnan(V.vgridded,0) < 0), 1), 2);

    fprintf('S0 = %2.1f | Tvi = %.2f mSv | Tsi = %.2f mSv \n', ...
            S0, Tvi/1000, Tsi/1000);
end
hax(2).XDir = 'reverse';
hax(2).XAxisLocation = 'top';
legend(cellstr(num2str([33:0.5:35]')), 'Location', 'SouthEast')
beautify
linkaxes(hax, 'y');
%export_fig -r150 -a2 images/leebrink2010.png