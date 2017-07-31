
%% flux time series
if ~exist('ew8341', 'var')
    if exist('sh', 'var')
        ew8341 = sh.array(5);
    else
        ew8341 = runs('../topoeddy/runew-8341/');
    end
end

handles = ew8341.plot_fluxts(1, 1, 1);
handles.htitle.String = ['Offshore flux of shelf water at the shelfbreak'];
handles.htitle.FontSize = 28;
maximize; %drawnow;
          % mark timesteps
axes(handles.hax(1));
ylim([0 12]*1e3);
handles.hax(1).YTick(handles.hax(1).YTick == 5e3) = [];
hold on
handles.htstep = ...
    plot(handles.ts(1).XData(timesteps), handles.ts(1).YData(timesteps), ...
         '.', 'MarkerSize', 28);
linkprop([handles.ts(1) handles.htstep], 'Color');

handles.leghandles = [handles.leghandles handles.htstep];
legstr = handles.hleg.String;
legstr{end+1} = 'Time instants shown in Figure 2';
[handles.hleg, handles.icons] = legend(handles.leghandles, legstr, ...
                                       'Location', 'NorthWest', 'Box', 'off');
hpt = findobj(handles.icons, 'Type', 'patch');
hpt.FaceAlpha = 0.2;

export_fig -r150 -a2 -png images/paper3/flux-diags

%% ubar scale
if ~exist('shfric2', 'var')
    folders = { ...
        'runew-8341', ...
        'runew-583413', 'runew-583411', ...
        'runew-583414', 'runew-583415', ...
              };
    shfric2 = runArray(folders);
end

shfric2.sort('run.params.misc.rdrg');
shfric2.name = cellstr(num2str(shfric2.print_params('run.params.misc.rdrg')', '%1.0e'));
shfric2.name{1} = '0';
handles = shfric2.plot_ts('-run.ubarscale.scale/1000', [], ...
                          'run.ubarscale.time/86400');
ylabel('Distance from shelfbreak (km)');
xlabel('Time (days)');
delete(handles.htind);
delete(handles.hmaxloc);
set(gcf, 'Position', [190 128 1126 810]);
hleg = legend;
hleg.Location = 'SouthEast';
ylim([-40 0]);
htxt = text(0.85, 0.45, 'r_f (m/s)', 'Units', 'normalized');
title('Cross-isobath extent of along-shelf supply jet');
xlim([120 320]);
linex([190 230]);
correct_ticks('x', [], {'200'});

resizeImageForPub('portrait');
pbaspect([2 1 1]);
beautify([12 13 14], 'Times');

export_fig -painters images/paper3/shfric-ubarscale.pdf

%% bottom density anomaly
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-8342-2')
    ew = runs('../topoeddy/runew-8342-2/');
end

tindex = find_approx(ew.time, 311*86400,1);
handles = ew.overlay_section('v', 'rho', tindex, {'y' 1 60}, 's', 1);
handles.h_plot2.LevelList = linspace(21.78, 21.84, 26);
center_colorbar;
[csdbot,xx,yy,~] = dc_roms_read_data(ew, ew.csdname, ...
                                     tindex, {'y' 1 60; 'z' 1 1});
hold on;
[~,h2] = contour(xx/1000,yy/1000,csdbot < ew.bathy.xsb, [1 1], ...
                 'Color', ew.shelfSlopeColor, 'LineWidth', 3);
xlim([230 400]);
title({'Cross-isobath velocity (color, m/s), \rho contours (black) and';  ...
       'shelf-water front (blue contour) at the bottom'});
xlabel('X (km)');
ylabel('Y (km)');
beautify;

export_fig -opengl -r150 -a4 images/paper3/vbot-rhobot.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cyclone
cyc = runs('../topoeddy/runew-34-cyc/');

opt.addvelquiver = 0;
opt.quiverloc = 'bot';

%cyc.animate_field('csdye', [], '300', 1, opt);

clf;
clear handles
hax = packfig(2,2);
depth = -[30 100];
labels = 'ab';
for ii=1:2
    handles{ii} = cyc.animate_zslice('eddye', depth(ii), '300', hax(ii), opt);
    colorbar('off');
    ylim([150 300]);
    xlim([100 510]);
    handles{ii}.htime.String = [labels(ii) ') z = ' num2str(depth(ii)) ' m'];
    handles{ii}.htime.Position(2) = 0.1;
    title('');
    hl(ii) = linex(320);
end

hax(2).YTickLabels = {};
hax(2).YLabel.String = '';
hax(2).XLabel.String = '';
correct_ticks('x', [], {'300'; '500'}, hax(1));
correct_ticks('x', [], {'300'}, hax(2));
correct_ticks('y', [], {'200'}, hax(1));

handles{3} = mod_movie(cyc.dir, 'dye_03', 301, {}, 's', 1, '', hax(3));
handles{3}.h_plot.EdgeColor = 'none';
cyc.plot_bathy('contour');
ylim([150 300]);
xlim([100 510]);
colorbar('off');
ylabel('Y (km)');
xlabel('X (km)');
htxt = text(0.05, 0.1, 'c) at the bottom', 'Units', 'normalized');
title('');
linex(320);
correct_ticks('x', [], {'300'}, hax(3));
hax(3).YTickLabelMode = 'auto';
correct_ticks('y', [], {'200'; '250'}, hax(3));

handles{4} = cyc.PlotSingleYZSection('eddye', '300', '320e3', hax(4));
xlim([150 300]);
colorbar('off');
pbaspect([1.4 1 1]);
hax(4).Position(1) = 0.59;
title('');
handles{4}.hlzlabel.delete;
handles{4}.htlabel.Position(1) = 0.6;
handles{4}.htlabel.Position(2) = 0.4;
correct_ticks('y', [], {'-312'}, hax(4));
correct_ticks('x', [], {'231'}, hax(4));
htxt(2) = text(0.055, 0.1, 'd)', 'Units', 'normalized');

htitle = suplabel('Eddy dye sections at t = 300 days for a cyclone', 't');
htitle.Position(4) = 0.8;
htitle.FontSize = 28;

export_fig -a2 -r150 images/thesis/cyclone-sections.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PV/RV
ew = runs('../topoeddy/runew-34/');
ew.read_pvsurf;

clear handles
opt.csdcontours = 0;
opt.limy = [0 150];
opt.drawtrack = 0;
opt.drawcenter = 0;
opt.rhocontourplot = 0;

figure;
hax = packfig(3,1);

handles(1) = ew.animate_field('csdye', hax(1), 300, 1, opt);
title(''); xlabel('');
htxt(1) = text(0.05, 0.1, 'a) Cross-shelf dye', 'Units', 'Normalized', 'Color', 'white');

handles(2) = ew.animate_field('pv', hax(2), 350, 1, opt);
xlabel(''); title('');
htxt(2) = text(0.05, 0.1, 'b) Surface PV (x 10^{-11})', 'Units', 'Normalized', 'Color', 'white');
caxis([0.4e-11 9.4e-11]/1e-11);
handles(2).hfield.CData = handles(2).hfield.CData/1e-11;

handles(3) = ew.animate_field('rv', hax(3), 300, 1, opt);
htxt(3) = text(0.05, 0.1, 'c) Vorticity/f_0', 'Units', 'Normalized');
title('');
caxis([-1 1]*0.4);

for ii=1:2
    hax(ii).XTickLabel = {};
    correct_ticks('y', [], {'0'}, hax(ii));
end

linkaxes(hax, 'xy');

export_fig -r150 -a2 images/thesis/ew-34-csdye-pv-rv.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xymap
name = 'ew-64361';
if ~exist('run', 'var') || ~strcmpi(run.name, name)
    run = runs(['../topoeddy/run' name '/']);
end
fontSize = [20 22 24];
ms = 12; % marker size
trackcolor = [1 1 1]*0.65;
sbslcolor = trackcolor;
tt = [1 250];
seqcolor = flipud(cbrewer('div','RdYlBu',32));
figure; maximize(); pause(0.5);

ax3 = subplot(121);
run.animate_field('eddye', ax3, tt(1), 1);
limx = xlim;
plot(run.eddy.mx/1000, run.eddy.my/1000, 'Color', trackcolor);
plot(run.eddy.mx(tt(1))/1000, run.eddy.my(tt(1))/1000, 'x', ...
     'MarkerSize', ms, 'Color', trackcolor);
text(0.15*diff(limx), run.bathy.xsl/1000, 'slopebreak', ...
     'VerticalAlignment', 'Bottom', 'FontSize', fontSize(1)-4, ...
     'Color', sbslcolor);
text(0.15*diff(limx), run.bathy.xsb/1000, 'shelfbreak', ...
     'VerticalAlignment', 'Top', 'FontSize', fontSize(1)-4, ...
     'Color', sbslcolor);
title('Dyes and SSH'); caxis([-1 1]); beautify(fontSize);
colormap(ax3,seqcolor);
correct_ticks('y',[],[3 6]);

ax4 = subplot(122);
run.animate_field('eddye', ax4, tt(2), 1);
plot(run.eddy.mx/1000, run.eddy.my/1000, 'Color', trackcolor);
plot(run.eddy.mx(tt(2))/1000, run.eddy.my(tt(2))/1000, 'x', ...
     'MarkerSize', ms, 'Color', trackcolor);
caxis([-1 1]);
colormap(ax4,seqcolor); beautify(fontSize);
title('Dyes and SSH'); ylabel([]);
ax4.YTickLabel = [];
correct_ticks('y',[],[3 6]);

export_fig('-r450','images/xymap-poster.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

tic; export_fig('-r450','images/paper1/sb-maps.png'); toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image effect?
co = image.sorted_colors;
figure; maximize();
ax(2) = subplot(223); hold on;
ax(3) = subplot(224); hold on;
image.flip_colors = 1;

% track
ax(1) = subplot(2,2,[1 2]);
image.plot_penetration(gca);
title([]);
% add timestamps
text(-1.85, 2.75, '50d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-1.85*[1 1], [2.75 2.3], 'k-');
text(-3.45, 2.75, '100d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-3.45*[1 1], [2.75 1.9], 'k-');
text(-3.45, 2.75, '100d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-3.45*[1 1], [2.75 1.9], 'k-');
text(-5.45, 2.75, '200d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-5.45*[1 1], [2.75 1.5], 'k-');

% volume, speed diagnostics
co = image.sorted_colors;
for ii=1:image.len
    run = image.array(ii);
    tvec = run.time/run.eddy.turnover;
    nsmooth = 3; run.eddy.turnover./diff(run.time(1:2));

    subplot(223);
    cvx = run.eddy.mvx;
    %cvx(cvx < -0.06) = NaN;
    plot(tvec, smooth(cvx, 8));
    ylabel({'Eddy center along-isobath', 'velocity (km/day)'});

    subplot(224);
    plot(tvec, run.eddy.vol(:,1)./run.eddy.vol(1,1));
    ylabel('Volume / Volume(t=0)');
    xlabel('Time / Turnover Time');
end
subplot(223);
hl = liney(0); ylim([-0.08 0.05]);
uistack(hl, 'bottom'); axis tight;
ax(2).YTick = sort(unique([ax(2).YTick -2.3 -1.5 -0.8]));
correct_ticks('y', '%.1f', {'-2','-1'});
beautify; pbaspect([1.618 1 1]);

subplot(224);
ylim([0 1]);
%liney([0.75 0.4 0.3 0.25]);
ax(3).YTick = sort(unique([ax(3).YTick 0.3  0.75]));
ax(3).XTick = sort(unique([ax(3).XTick 70]));
correct_ticks('y', '%.2f', '0.8');
beautify; pbaspect([1.618 1 1]);
image.reset_colors(co);

linkaxes(ax(2:3), 'x');
export_fig('images/paper1/image-effect.pdf');

%% image effect velocity - sb
image = runArray({ 'runew-2360-fb', 'runew-2360-20km', ...
                   'runew-2360', 'runew-2360_wider'});
for ii=1:image.len
    sh(ii) = image.array(ii).params.bathy.L_shelf;
    if image.array(ii).params.flags.flat_bottom
        sh(ii) = 0;
    end
    image.name{ii} = [num2str(sh(ii)/1000, '%.0f') ' km shelf'];
end
image.sort(sh);
