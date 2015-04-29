%% Figures for paper 1:
% EW isobaths
ewall = runArray({ ...
    'runew-4341/', ...
    'runew-36-sym/', ...
    'runew-6362-2', ...
    'runew-64361/', ...
                });
ewall.name = {'L_{edd}/L_{sl} > 1', 'L_{edd}/L_{sl} ~ 1', ...
              'L_{edd}/L_{sl} < 1', 'L_{edd}/L_{sl} < 1'};

ew = runArray({ ...
    'runew-64361', ...
    'runew-6341', ...
    'runew-6362-2', ...
    'runew-6441' });

% shelfbreak depth
sb = runArray({ ...
    'runew-36', 'runew-2360', 'runew-2361_wider', 'runew-2362_wider', ...
              });

% wide slope runs
folders = { ...
%'runew-04', 'runew-6040', 'runew-6041', 'runew-6042-new', ...
    ... %'runew-15', 'runew-6150', 'runew-6151', ... % 6152 sucks
    ... %'runew-34', 'runew-35', 'runew-36', ...
    'runew-6341', 'runew-6342', ...%'runew-6343', ...
    ...%'runew-6353', ...
    'runew-6362', 'runew-6362-1', 'runew-6362-2', ...
    'runew-6371', 'runew-6372', ...% 'runew-6373' ...
    'runew-6452', 'runew-6352', 'runew-64361', 'runew-64361-2', ...
    'runew-64461-3', 'runew-64462', ...
    'runew-b4361', ...
    'runew-64351-cyc', 'runew-6341-fneg', ...
    ... %'runew-05', 'runew-6052', ...
    ... %'runew-06', 'runew-6062', ...
          };
sl = runArray(folders);

% image effect velocity - sb
image = runArray({ 'runew-2360-fb', 'runew-2360', 'runew-2360_wider'});

% NS isobaths
folders = { ...
    'runns-64361', 'runns-6341', 'runns-6362-2',...
    'runns-6441', ...
          };
ns = runArray(folders);

%% initial condition
run = ew.array(2);

%% x-y cross-section
run = runs('../topoeddy/runew-34/');
fontSize = [16 16 18];
tt = [1 350];
figure; ax1 = subplot(221);
run.animate_vorsurf(ax1,tt(1),1);
title('Surface vorticity / f');
xlabel([]);
clim = caxis; beautify(fontSize);
ax2 = subplot(223);
run.animate_vorsurf(ax2,tt(2),1);
title([]); xlabel([]);
caxis(clim); beautify(fontSize);
seqcolor = cbrewer('seq','Reds',12);
ax3 = subplot(222);
run.animate_field('eddye', ax3, tt(1), 1);
title('Dye'); caxis([0 1]); beautify(fontSize);
colormap(ax3,seqcolor);
ax4 = subplot(224);
run.animate_field('eddye', ax4, tt(2), 1);
title([]); caxis([0 1]);
colormap(ax4,seqcolor); beautify(fontSize);
export_fig('images/paper1/xymap.pdf');

%% EW - center tracks - all,
ewall.plot_penetration; maximize();
export_fig('images/paper1/centrack.pdf');

%% EW, NS Center-tracks - wide slope
ew.filter = [];
ns.filter = [];

figure; maximize();
ax1 = subplot(121);
ew.plot_penetration(ax1);
subplot(122);
ns.plot_penetration(gca); drawnow;
ax1 = gca; ax1.XTick = unique([ax1.XTick 1])
export_fig('images/paper1/sl-centrack.pdf');

%% parameterization
sl.print_diag('bottom torque');
title([]); pause(1);
export_fig('images/paper1/penetration-erf-param.pdf');

%% energy decay
sl.plot_dEdt; maximize(); pause(1);
subplot(121); title([]);
pbaspect([1.618 1 1]); xlim([0 400]);
legend('off');
subplot(122); xlim([0 400]);
pbaspect([1.618 1 1]);
export_fig('images/paper1/energy-decay.pdf');

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

%% image effect?
figure; maximize(); hold on;
for ii=1:image.len
    run = image.array(ii);

    tvec = run.time/run.eddy.turnover;
    cvx = run.eddy.cvx;
    cvx(cvx < -0.06) = NaN;
    plot(tvec, cvx);
    ylabel({'Centroid', 'along-isobath', 'velocity', '(km/day)'});
    xlabel('Time / Turnover Time');
end
hl = liney(0); ylim([-0.08 0.05]);
uistack(hl, 'bottom');
legend('0 km shelf', '40 km shelf', '150 km shelf');
beautify;
export_fig('images/paper1/image-effect.pdf');