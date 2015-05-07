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

% NS isobaths
folders = { ...
    'runns-64361', 'runns-6341', 'runns-6362-2',...
    'runns-6441', ...
          };
ns = runArray(folders);

% shelfbreak depth - xy maps
sb = runArray({ ...
    'runew-36', 'runew-2360_wider', 'runew-2361_wider', ...
    'runew-2363_wider', 'runew-2362_wider', ...
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
image = runArray({ 'runew-2360-fb', 'runew-2360-20km', ...
                   'runew-2360', 'runew-2360_wider'});

%% bottom friction
folders = { ...
    ... %'runew-34', 'runew-5341', 'runew-5343'
    'runew-6341', 'runew-56341', 'runew-56341-2', ...
          };
bfrics = runArray(folders);
for ii=1:bfrics.len
    run = bfrics.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    bfrics.name{ii} = ['r = ' ...
                       num2str(bfrics.array(ii).params.misc.rdrg) ' m/s'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial condition
run = ew.array(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x-y cross-section
if ~exist('run', 'var') || ~strcmpi(run.name, 'ew-34')
    run = runs('../topoeddy/runew-34/');
end
fontSize = [16 16 18];
tt = [1 250];
seqcolor = flipud(cbrewer('div','RdYlBu',32));
figure; maximize(); pause(0.5);
ax1 = subplot(221);
run.animate_vorsurf(ax1,tt(1),1);
plot(run.eddy.mx/1000, run.eddy.my/1000, 'k');
plot(run.eddy.mx(tt(1))/1000, run.eddy.my(tt(1))/1000, 'kx');
title('Surface vorticity / f');
xlabel([]);
clim = caxis; beautify(fontSize);
correct_ticks('y',[],[3 5]);

ax2 = subplot(223);
run.animate_vorsurf(ax2,tt(2),1);
plot(run.eddy.mx/1000, run.eddy.my/1000, 'k');
plot(run.eddy.mx(tt(2))/1000, run.eddy.my(tt(2))/1000, 'kx');
title([]); xlabel([]);
caxis(clim); beautify(fontSize);
correct_ticks('y',[],[3 5]);

ax3 = subplot(222);
run.animate_field('eddye', ax3, tt(1), 1);
plot(run.eddy.mx/1000, run.eddy.my/1000, 'k');
plot(run.eddy.mx(tt(1))/1000, run.eddy.my(tt(1))/1000, 'kx');
title('Dye'); caxis([-1 1]); beautify(fontSize);
colormap(ax3,seqcolor);
correct_ticks('y',[],[3 5]);

ax4 = subplot(224);
run.animate_field('eddye', ax4, tt(2), 1);
plot(run.eddy.mx/1000, run.eddy.my/1000, 'k');
plot(run.eddy.mx(tt(2))/1000, run.eddy.my(tt(2))/1000, 'kx');
title([]); caxis([-1 1]);
colormap(ax4,seqcolor); beautify(fontSize);
correct_ticks('y',[],[3 5]);

export_fig('-r450','images/paper1/xymap.png');

%% EW - center tracks - all,
fs  = 18;
ewall.plot_penetration; maximize();
ax = gca;
title([]); pbaspect([1.618 1 1]);
text(-10, 0.2, 'shelfbreak', 'FontSize', fs);
text(-0.3, 1.5, 'L_{edd} > L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(1,:));
text(0.2, 3.8, 'L_{edd} ~ L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(2,:));
text(-9, 2.1, 'L_{edd} < L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(3,:));
text(-6, 3.8, 'L_{edd} < L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(4,:));
%legend('off');
export_fig('images/paper1/centrack.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom friction

bfrics.plot_penetration; maximize();
pbaspect([1.618 1 1]);
export_fig('images/paper1/bfrics-centrack.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameterization
sl.print_diag('bottom torque');
title([]); pause(1);
export_fig('images/paper1/penetration-erf-param.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% energy decay
sl.plot_dEdt; maximize(); pause(1);
subplot(121); title([]);
pbaspect([1.618 1 1]); xlim([0 400]);
legend('off');
subplot(122); xlim([0 400]);
pbaspect([1.618 1 1]);
export_fig('images/paper1/energy-decay.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x-y plots of shelfbreak runs
% add colorbar
% figure out best time indices
tinds = [90 90 100 100]; % figure out what these are
figure; maximize();
sb.filter = [2 4 3 5];
fontSize = [16 16 18];
len = length(sb.filter);
for ii=1:len
    run = sb.array(sb.filter(ii));
    run.read_zeta;

    ix = run.spng.sx1:run.spng.sx2;
    iy = run.spng.sy1:run.spng.sy2;

    xr = run.rgrid.x_rho(iy,ix)'/1000;
    yr = run.rgrid.y_rho(iy,ix)'/1000;

    subplot(floor(len/2),2,ii);
    run.animate_field('eddye', gca, tinds(ii), 1);
    %pcolorcen(xr, yr, run.zeta(ix,iy,tinds(ii)));
    hold on;
    plot(run.eddy.mx/1000, run.eddy.my/1000, 'k');
    plot(run.eddy.mx(tinds(ii))/1000, run.eddy.my(tinds(ii))/1000, ...
         'kx');
    if ii == 1, clim = caxis; end
    hold all
    caxis(clim);

    title(['\lambda = ' num2str(run.bathy.hsb/run.eddy.Lgauss(1), ...
                                2) ...
          ' (H_{sb} = ' num2str(run.bathy.hsb, '%.0f') ' m)']);
    ylabel('Y (km)');
    xlabel('X (km)');
    if ii == 3
        correct_ticks('y',[],[4 6]);
    else
        correct_ticks('y',[],[5]);
    end
    beautify(fontSize);
end

%tic; export_fig('images/paper1/sb-maps.pdf'); toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image effect?
figure; maximize(); hold on;
for ii=1:image.len
    run = image.array(ii);
    tvec = run.time/run.eddy.turnover;
    nsmooth = 3; run.eddy.turnover./diff(run.time(1:2));
    cvx = run.eddy.cvx;
    cvx(cvx < -0.06) = NaN;
    plot(tvec, smooth(cvx, 6*nsmooth));
    ylabel({'Centroid along-isobath', 'velocity (km/day)'});
    xlabel('Time / Turnover Time');
end
legend('0 km shelf', '20 km shelf', '40 km shelf', '150 km shelf');
hl = liney(0); ylim([-0.08 0.05]);
uistack(hl, 'bottom'); axis tight;
beautify; pbaspect([1.618 1 1]);
export_fig('images/paper1/image-effect.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%