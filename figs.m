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