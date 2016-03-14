%%
if ~exist('ewshdp', 'var')
    ewshdp = runArray({'ew-64461-9-shallow', 'ew-64461-9-deep', ...
                       'ew-64361-shallow', 'ew-64361-deep'});
end

figure; maximize;
ax(1) = subplot(121);
ewshdp.filter = [1 2];
ewshdp.plot_ts('eddy.hcen', ax(1))
pbaspect([1.618 1 1]);
title('Rh = 20')
legend('off');

ax(2) = subplot(122);
ewshdp.filter = [3 4];
ewshdp.plot_ts('eddy.hcen', ax(2));
pbaspect([1.618 1 1]);
title('Rh = 60');
legend('off');

linkaxes(ax, 'y');
packfig(1,2, 'columns');
ax(2).XTick(1) = [];

[ax(3), htitle] = suplabel('Water depth at eddy center', 't');
ax(3).Position(4) = 0.72;
beautify;

export_fig -r150 -transparent images/shallow-deep-hcen.png

%% ew-34 vertprofile
% density profile changes with time as it adjusts to no-flux boundary condition. For some
% runs (ew-34 etc). this is on the order of Δρ seen due to up/downwelling in cyclone. So, ignore!!!

figure; maximize;
ax(1) = subplot(121);
ew34.plot_vertProfile('rho', 250, 50, [1 20 40 60 80 100], ax(1));
ax(1).Title.String = ['On the slope | ' ax(1).Title.String];
ax(2) = subplot(122);
ew34.plot_vertProfile('rho', 50, 210, [1 100 200 300 350], ax(2));
ax(2).Title.String = ['Just outside sponge | ' ax(2).Title.String];
export_fig -a1 images/ew-34-vert-profile-rho.png

%% ew-2360 vertprofile
% density profile changes with time as it adjusts to no-flux boundary condition. For some
% runs (ew-34 etc). this is on the order of Δρ seen due to up/downwelling in cyclone. So, ignore!!!

figure; maximize;
ax(1) = subplot(121);
ew2360.plot_vertProfile('rho', 250, 170, [1 20 40 60 80 100], ax(1));
ax(1).Title.String = ['On the slope | ' ax(1).Title.String];
ax(2) = subplot(122);
ew2360.plot_vertProfile('rho', 35, 345, [1 100 200 300 350], ax(2));
ax(2).Title.String = ['Just outside sponge | ' ax(2).Title.String];
export_fig -a1 images/ew-2360-vert-profile-rho.png

%% ew-8383 along-shelf velocity
sh.plot_avgProfile('ubar', 'x', 'res', 1, []);