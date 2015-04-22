%% Figures for paper 1:
% EW isobaths
ew = runArray({ ...
    'runew-4341/', ...
    'runew-36-sym/', ...
    'runew-6362-2', ...
    'runew-64361/', ...
                });
ew.name = {'L_{sl}/L > 1', 'L_{sl}/L ~ 1', 'L_{sl}/L < 1', 'L_{sl}/L < 1'};

% shelfbreak depth
sb = runArray({ ...
    'runew-36', 'runew-2360', 'runew-2361_wider', 'runew-2362_wider', ...
              });

% NS isobaths
folders = { ...
    'runns-64361', 'runns-6341', 'runns-6362-2',...
    'runns-6441', ...
          };
ns = runArray(folders);

%% EW, NS Center-tracks
ew.filter = [];
ns.filter = [];

figure; maximize();
subplot(121);
ew.plot_penetration(gca);
subplot(122);
ns.plot_penetration(gca);
export_fig('images/paper1/centrack.pdf');
%export_fig('images/paper1/EWcentrack.pdf')
%export_fig('images/paper1/NScentrack.pdf');

%% x-y plots of shelfbreak runs
% add colorbar
% fix bathy
% figure out best time indices
% change ticks to mark shelfbreak and slopebreak exactly.
tinds = [200 100 150 100]; % figure out what these are
figure;
for ii=1:sb.len
    run = sb.array(ii);
    run.read_zeta;

    xr = run.rgrid.x_rho'/1000;
    yr = run.rgrid.y_rho'/1000;

    subplot(sb.len/2,2,ii);
    pcolorcen(xr, yr, run.zeta(:,:,tinds(ii)));
    if ii == 1, clim = caxis; end
    hold all
    run.plot_bathy('contour');
    caxis(clim);

    title(['\lambda = ' num2str(run.bathy.hsb/run.eddy.Lgauss(1), ...
                                2) ...
          ' (t = ' num2str(run.eddy.t(tinds(ii)), '%.0f') ' days)']);
    ylabel('Y (km)');
    xlabel('X (km)');
end

export_fig('images/paper1/sb-maps.pdf');