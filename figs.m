%% sbssh plots
for ii=1:sh.len
    %sh.array(ii).ShelfbreakSSHScale;
    sbssh = sh.array(ii).sbssh; %
    figure; maximize;
    plot(sbssh.fitobj, sbssh.xvec, sbssh.ssh - mean(sbssh.ssh), 'predobs');
    hax = gca;
    hax.XAxis.Exponent = 3;
    title([sh.array(ii).name ' | ' num2str(sbssh.X*2/1000) ' km']);
    legend('Location' ,'SouthEast');

    %sh.array(ii).avgSupplyJet;
    export_fig(['images/sbssh-' sh.array(ii).name '.png']);
end

%%
% 8383, 8352, 8392

figure; hold on; maximize;
for ii=[5 6 9]
    plot(sh.array(ii).sbssh.xvec, sh.array(ii).sbssh.ssh./sh.array(ii).sbssh.fitobj.y0, ...
         'DisplayName', sh.array(ii).name);
end
legend;

%% ew-34
if ~exist('ew34', 'var')
    ew34 = runs('../topoeddy/runew-34/')
    ew34.read_csdsurf;
end
ew34.makeVideo = 1;
ew34.video_init('csd-plain');
handles = ew34.plot_surf('csdsurf', 'contourf', 330);
handles.EdgeColor = 'none';
hax = gca;
hax.XAxis.Color = 'none';
hax.YAxis.Color = 'none';
colormap(cbrewer('seq', 'Greys', 32));
caxis([30 160]*1e3);
ew34.video_update;
for tt=330:350
    ew34.update_surf('csdsurf', handles, tt);
    ew34.video_update;
end
ew34.video_write;

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

%% sloping shelf runs
sh.filter = [6 7];
sh.plot_avgProfile('u', 'x', 350, 0);


folders = { ...
    'runew-04', 'runew-05', 'runew-06', ...
    'runew-34', 'runew-35', 'runew-36', 'runew-37', ...
    'runew-6370', ... %'runew-6341', 'runew-6371', ...
    ... %'runew-2041', 'runew-2043', ...
    'runew-2340', 'runew-2345', ...
    ... %'runew-82340', ...
    ... %'runew-2341_wider','runew-2361_wider', 'runew-2363_wider',  ...
    'runew-2360_wider', 'runew-2365-75km', ...
     ... %'runew-2362_wider', 'runew-2360-20km', ...
    ... %'runew-4341', 'runew-3340-80km' , 'runew-3341-2-big', ...
    'runew-23340', 'runew-24340', ...
    'runew-4040', 'runew-4050', ...
    'runew-4343', ... %'runew-4344', ...
    'runew-5341', 'runew-5343', ...
    'runew-5040', 'runew-5041', 'runew-5043', ...
    ... %'runew-b4361', ...
    'runew-8041', 'runew-8042',' runew-8150', 'runew-8151', ...
    ... %'runew-82342', 'runew-82343', 'runew-82344' ...
    'runew-8341', 'runew-8352', ...
    'runew-583411', 'runew-583413', ...
    'runew-8383', ...
          };

if ~exist('csf', 'var'), csf = runArray(folders); end

hax = csf.plot_fluxparam('avg flux', 1, 'no_sloping_shelf');

hax(2).Title.String = 'Integrated to shelfbreak depth';
hax(2).Title.Interpreter = 'latex';
hax(2).Title.FontSize = 26;
hax(4).YLabel.FontSize = 24;
hax(8).XLabel.FontSize = 24;
dx = 100;
xmax = 450;
hax(2).XTick = [0:dx:xmax];
for ii=[2 4:9]
    hax(ii).XLim = [0 xmax];
    hax(ii).XTick = [100:dx:xmax];
end
hax(7).XTick = [0 hax(7).XTick];
hax(1).XColor = hax(1).YColor;
hax(1).XTickLabelMode = 'auto';
hax(1).XAxisLocation = 'top';
hax(1).YLim = [0 120];
hax(1).YTick = [0:20:100];
hax(2).YTick = [0:20:100];
hax(end).XTick(end) = [];

axes(hax(3));
hleg = legend;
legstr = hleg.String;
hleg.delete;
correct_ticks('y', '%.2f');
hax(3).YTick(end-1) = [];
hax(3).Position(1) = 0.68;
htxt(1) = text(1.25,0.30,'Slope (m)', 'FontSize', 18);
htxt(2) = text(1.25,0.03,'y-intercept (c)', 'FontSize', 18, ...
               'Color', [1 1 1]*0.55);
ylim([-0.1 0.35]);

% export_fig -r150 -a2 -png -pdf -opengl images/paper2/avgflux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ew-2360 z-slice mosaic

tinds = [63 70 77 95];
handles = ew2360.mosaic_zslice('dye_03', 200, tinds);
for ii=1:length(handles)
    handles(ii).csdsurf.LevelList = [170] * 1e3;
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([265 465]);
ylim([130 250]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 100m | Ro = 0.25';
handles(1).htitle.FontSize = 22;
handles(1).htitle.FontWeight = 'normal';
handles(1).supax.Position(4) = 0.88;
axes(handles(1).hax); correct_ticks('y', '', {'150' '200'});
axes(handles(3).hax); correct_ticks('y', '', {'150' '200'});

xx = [343 346 349]';
yy = [231 237 231]';


colors = brighten(cbrewer('qual', 'Paired', length(xx)), -0.5);
axes(handles(4).hax);
for ii=1:length(xx)
    hpnt(ii) = plot(xx(ii), yy(ii), 'x', 'Color', colors(ii,:), ...
                    'MarkerSize', 16);
end
hpnt(end+1) = plot(ew2360.eddy.mx(tinds(end))/1000, ew2360.eddy.my(tinds(end))/1000, ...
                   'x', 'Color', [1 1 1]*0.45, 'MarkerSize', 16);

export_fig -r150 -a2 images/paper2/ew-2360-mosaic-zslice.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom rho PV

figure;
hplt = ew2360.plotBottomRhoPV(70);
hplt.h_bathy{1}.LineStyle = '--';
hplt.h_bathy{1}.LineWidth = 1;
hplt.h_bathy{2}.LineStyle = '--';
title('Bottom PV with \rho contours');
xlabel('X (km)');
ylabel('Y (km)');
hplt.htext = text(330, 149, 'shelfbreak', 'Color', hplt.h_bathy{2}.Color)
ylim([140 180]);
set(gca, 'YTickMode', 'auto');

export_fig -r150 -a2 images/paper2/ew-2360-bottomrhopv.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cyclone vs. anti-cyclone flux vertical profiles

clear handles;
if ~exist('fluxvert', 'var')
    fluxvert = runArray({'runew-2360_wider', 'runew-34-cyc'});
end

day = {'119'; []; '75'};
names = {'Anticyclone'; 'Cyclone'}; %; 'Anticyclone - flat bottom'};
opt.rhocontours = 1;
figs(2) = 1;
fluxvert.plot_fluxes(4, [],[],figs);
title('');
legend(names, 'Location', 'SouthEast')
set(gcf, 'Position', [1 1 800 821]);
beautify;
hax(4) = gca;
subplot(2,2,[2 4],hax(4));
axes(hax(4));
pbaspect([1 1 1]);
hax(4).Position(1) = 0.65;
hax(4).Position(3) = 0.28;

for ii=1:2
    hax(ii) = subplot(2,2,1 + 2*(ii-1));
    handles(ii) = ...
        fluxvert.array(ii).PlotSingleXZSection('v', 4, day{ii}, opt, hax(ii));
    title(names{ii});
    if ii ~= fluxvert.len
        hax(ii).XTickLabel = {};
        xlabel('');
    end
    handles(ii).htime.Position(2) = 0.11;
end
hax(2).Title.Color = hax(4).Children(1).Color;
hax(1).Title.Color = hax(4).Children(2).Color;

linkaxes(hax(1:fluxvert.len), 'xy')
ylim([-400 0]);
xlim([-100 100]);

axes(hax(1));
correct_ticks('y', [], '-100');
axes(hax(2));
correct_ticks('y', [], {'-100'; '-300'});

handles(1).hcb.Position(1) = 0.5;
handles(2).hcb.Position(1) = 0.5;
handles(2).hcb.Label.String = {'Cross-shelf'; 'velocity (m/s)'};
handles(2).hcb.Label.Position = [1.0 0.135 0];
handles(2).hcb.Label.Rotation = 0;

for ii=1:2
    handles(ii).htext{1}.Units = 'Normalized';
    handles(ii).htext{2}.Units = 'Normalized';
    handles(ii).htext{1}.Position(1) = 0.5;
    handles(ii).htext{2}.Position(1) = 0.5;
end
export_fig -r150 -a2 images/paper2/fluxvertprofile.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max flux

hax = csf.plot_fluxparam('max flux');
hax(2).Title.String = 'Integrated to shelfbreak depth';
for ii=[2 4:9]
    hax(ii).XLim = [0 300];
    hax(ii).XTick = [0:100:300];
end

axes(hax(3));
hleg = legend;
hleg.Position(1) = 0.82;
ylim([-1 1]*55);

export_fig -r150 -a2 images/paper2/maxflux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% western coast

if ~exist('ns', 'var')
    ns = runs('../topoeddy/runns-35/');
end

timesteps = [1 75 150 200 280 400];

opt.csdcontourplot = 0;
opt.csdcontours = ns.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.addzeta = 1;
handles = ns.mosaic_field('csdye', timesteps, opt);
for ii=1:6
    handles.hfield{ii}.hzeta.LevelList = linspace(0, 0.025, 6);
    handles.hfield{ii}.hzetaneg.LevelList = linspace(-0.01,0,6);
    handles.hfield{ii}.htlabel.Position(1) = 0.55;
    handles.hfield{ii}.htlabel.Position(2) = 0.06;
end
handles.hcb.delete;

ylim([min(ylim) 300]);
handles.supax.Position(4) = 0.88;
handles.htitle.String = 'Surface cross-shelf dye (km)';

axes(handles.hax(1));
[hleg,icons] = legend([handles.hfield{1}.hcen, ...
                    handles.hfield{1}.htrack, ...
                    handles.hfield{1}.hrho, ...
                    handles.hfield{1}.hzeta], ...
                      {'Eddy center', 'Center track', 'Eddy core', 'SSH'}, ...
                      'Location', 'NorthWest'); %, 'FontSize', 14);
hleg.Box = 'off';
hleg.Position(1) = 0.26;
hleg.Position(2) = 0.56;
icons(end).Children.Children(1).LineWidth = 1;
icons(end).Children.Children(2).LineWidth = 1;
icons(end).Children.Children(3).LineWidth = 1;

correct_ticks('x', [], {'50'; '100'}, handles.hax(4:6));

export_fig -painters -a2 images/paper2/ns-35-csdsurf.png


%% penetration on shelf
ew = sh.array(3);
opt.dt = 3;
opt.addzeta = 1;
opt.addvelquiver = 1;
[start,stop] = ew.flux_tindices(ew.csflux.off.slope(:,1,1), 0.3);
opt.limy = [0 ew.bathy.xsb/1000 + 40];
ew.animate_field('csdye', [], 80, 1, opt);

%% shelf supply diagnostics
env = ew.csflux.off.slopewater.envelope(:,1);
tvec = ew.csflux.time(~isnan(env));
env = env(~isnan(env));

tanh_fit(tvec/86400, ew.bathy.xsb/1000 - env/1000, 1);
hax = gca;
hax.Children(1).Visible = 'off';
hax.Children(2).Visible = 'off';
hax.Children(6).LineStyle = '-';
for ii=[3 4 5 6]
    hax.Children(ii).YData = hax.Children(ii).YData + 6.619;
end
ylabel('Max. distance of parcels from shelfbreak');
xlabel('Time (days)');
beautify([26 24 28]);
export_fig images/env.png

%%
if ~exist('xy', 'var')
    xy = runArray({'xy-64461-9-shallow', 'xy-64361'});
end

figure;
hold on;
handles = xy.array(1).plot_bathy('contour');
handles{1}.LevelList = [50 100 200 300 400 500 600];
handles{2}.delete;
handles{3}.delete;
for ii=1:2
    hplt(ii) = plot(xy.array(ii).eddy.mx/1000, xy.array(ii).eddy.my/1000);
end
legend(hplt,xy.name);
beautify
export_fig images/xy-eddytrack.png

opt.maskcontour = 0;
xy.array(1).PlotSingleYZSection('rho', '500');
xy.array(3).PlotSingleYZSection('csdye', '200');


%% f/h contours

clim = [0 1e-6];
figure;
run = xy.array(1);
contour(run.rgrid.xr(:,1)/1000, run.rgrid.yr(1,:)/1000, ...
          (run.rgrid.f'./run.bathy.h)', 'b');
hold on;
contour(run.rgrid.xr(:,1)/1000, run.rgrid.yr(1,:)/1000, run.bathy.h', 'k');
plot(run.eddy.mx/1000, run.eddy.my/1000);
legend('f/h', 'h', 'eddy track');
export_fig -a4 -r150 images/xy-shallow-f-h.png

% figure;
% run = xy.array(3);
% pcolorcen(run.rgrid.xr(:,1), run.rgrid.yr(1,:), ...
%           (run.rgrid.f'./run.bathy.h)');
% colorbar;
% caxis(clim);

xy.filter = [1 3];
xy.plot_ts('run.eddy.fcen/run.eddy.hcen''');


%%
figure;
hax(1) = subplot(221);
handles = xy.array(1).animate_field('csdye', hax(1), '200', 1);
handles.hbathy{3}.delete;
handles.hbathy{1}.LevelList = linspace(min(xy.array(1).bathy.h(:))+5, max(xy.array(1).bathy.h(:)), ...
                                       10);

hax(2) = subplot(222);
% cla
% run = xy.array(1);
% cross = hypot(avg1(run.ubar(:,2:end-1,:),1), avg1(run.vbar(2:end-1,:,:),2));
% contourf(run.rgrid.xr(2:end-1,2:end-1)/1000, run.rgrid.yr(2:end-1,2:end-1)/1000, ...
%          cross(:,:,102));
handles = xy.array(1).animate_field('csdye', hax(2), '300', 1);
handles.hfield.LevelList = [handles.hfield.LevelList 1e-3];
handles.hbathy{3}.delete;
handles.hbathy{1}.LevelList = linspace(min(xy.array(1).bathy.h(:))+5, max(xy.array(1).bathy.h(:)), ...
                                       10);

hax(3) = subplot(2,2,[3 4]);
plot(xy.array(1).eddy.t, xy.array(1).eddy.mvx*1000/86400);
liney(-1*0.004,'ubar on slope');
ylim([-0.03 0]);

%% shelfbc - sh from figs_paper3
sh.sort(sh.print_params('params.misc.rdrg'));
sh.filter = [5 15:18];
sh.sorted_colors;
bcind = 2;
figure;
hax = packfig(2,1);
for ff=1:length(sh.filter)
    ii = sh.filter(ff);
    run = sh.array(ii);
    axes(hax(1)); hold on;
    hplt1(ff) = plot(run.csflux.time/86400, run.csflux.off.slope(:,1,1), ...
                    'DisplayName', [run.name ' | r = ' num2str(run.params.misc.rdrg,'%.1e')]);

    axes(hax(2)); hold on;
    hplt2(ff) = plot(run.shelfbc.time/86400, run.shelfbc.shelf(:,bcind));
end

axes(hax(1)); legend(hplt1); beautify;
axes(hax(2)); beautify;

hplt1(1).Color = 'k';
hplt2(1).Color = 'k';
hax(1).YLabel.String = 'Shelf water offshore flux (m^3/s)';
hax(2).YLabel.String = ['BC_{' num2str(run.shelfbc.thresh(bcind), '%0.1f') '}'];
hax(2).XLabel.String = 'Time (days)';
hax(1).YLim(1) = 0;
hax(2).YLim = [0 1];
linkaxes(hax, 'x');
insertAnnotation('figs.m');
export_fig images/shelfbc-shfric.png
startup

%%
run = shfric2.array(5);

run.makeVideo = 0
opt.addvelquiver = 0;
opt.addzeta = 1;
opt.rhocontour = 0;
opt.commands = 'caxis([-1 1]*1e-3)';
run.animate_field('ubot', [],  '170', 1, opt);

run.animate_field('zeta', [],  '220', 1, opt);

%% sh bfric
sh.filter = [5 15:18];
shfric2.plot_ts(['(run.supply.zeta.zetameanUnmasked-run.supply.zeta.zetameanUnmasked(1))' ...
                 './(run.supply.zeta.zetameanUnmasked(end)-' ...
                 'run.supply.zeta.zetameanUnmasked(1))']);

%%
% for ii=1:5
%     shfric2.array(ii).read_velbar;
% end

figure;
clf; hold on;
shfric2.sorted_colors;
for ii=1:5
    run = shfric2.array(ii);
    %run.read_velbar;
    vec = nanmean(run.ubar(420,1:50,170:180),3);
    %vec = vec - vec(1);
    plot(vec, [1:50]);
end
liney([38 22]);
title('normalized ubar @ x = 425 km')
insertAnnotation('figs.m');
beautify;
legend(shfric2.name, 'Location', 'southwest');
export_fig images/shfric-ubar-downstream.png

%%
if ~exist('shfric', 'var')
    folders = { ...
        'runew-34', 'runew-5341', 'runew-8341', ...
        'runew-583413', 'runew-583411', ...
        ... %'runew-583414', 'runew-583415', ...
              };
    names = { ...
        '  0      0'; ...
        '  0    5e-3'; ...
        '  0.05     0'; ...
        '  0.05   5e-4'; ...
        '  0.05   5e-3'; ...
            };
    shfric = runArray(folders, names);
end

clear handles;
handles = shfric.PlotFluxSummary(1);

axes(handles.hax(1))
handles.hax(1).Position(2) = 0.70;
ylim([0 15]);
delete(handles.hleg);
handles.hleg(1) = legend(handles.hflux(1:2), shfric.name{1:2});
handles.hleg(2) = legend(handles.henv(3:5), shfric.name{3:5});
handles.hleg(1).Box = 'off';
handles.hleg(2).Box = 'off';
handles.hleg(1).Position(1) = 0.69;
handles.hleg(1).Position(2) = 0.79;
pos1 = handles.hleg(1).Position;
handles.hleg(2).Position(1) = handles.hleg(1).Position(1) + pos1(3) + 0.02;
handles.hleg(2).Position(2) = pos1(2) + pos1(4) - handles.hleg(2).Position(4);

tblx = [0.685 0.905] + 0.006;
htable(1) = text(1.21, pos1(2)+0.08, 'S_{sh}   r', 'Units', 'normalized', 'FontSize', 18);
htable(2) = text(1.45, pos1(2)+0.08, 'S_{sh}     r', 'Units', 'normalized', 'FontSize', 18);
htable(3) = annotation('line', tblx, [1 1]*0.890, 'LineWidth', 1);
htable(4) = annotation('line', tblx, [1 1]*0.848, 'LineWidth', 1);
htable(5) = annotation('line', tblx, [1 1]*0.760, 'LineWidth', 1);
linkprop(htable, 'Color');
htable(1).Color = handles.hax(1).XAxis.Color;
linkprop(handles.hleg, 'TextColor');
handles.hleg(1).TextColor = htable(1).Color;

% avg ssh
linkaxes(handles.hax, 'off');
linkaxes(handles.hax(1:2), 'x');
axes(handles.hax(3));
cla(handles.hax(3), 'reset');
handles2 = shfric.plot_avgProfile('zeta', 0, 'y', 'sb', 1, handles.hax(3));
handles2.htxt.delete;
xlim([-200 200]);
ax = handles.hax(3);
dy = ax.YTick(2);
for ii= [1 2]
    handles2.hpltsh(ii).YData = handles2.hpltsh(ii).YData + dy;
    handles2.hplted(ii).YData = handles2.hplted(ii).YData + dy;
end

handles2.hl.YData = ylim;
legend('off');
title(''); ylabel('');
handles2.htxt(1) = text(0.05,0.85, 'Mean SSH at shelfbreak (m)', ...
                     'Units', 'Normalized', 'FontSize', handles.htxt(1).FontSize);
handles2.htxt(2) = text(-150, 7e-4, 'Flat shelf', 'Units', 'data', ...
                        'HorizontalAlignment', 'center', 'FontSize', 18);
handles2.htxt(3) = text(-150, 2e-4, 'Sloping shelf', 'Units', 'data', ...
                        'HorizontalAlignment', 'center', 'FontSize', 18);
handles2.htxt(4) = text(-70, 1e-3, 'Shelf water', 'Units', 'data', ...
                        'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold');
handles2.htxt(5) = text(150, 1e-3, 'Eddy water', 'Units', 'data', ...
                        'HorizontalAlignment', 'center',  'FontSize', 18);
linkprop(handles2.htxt(4:5), 'Color');
handles2.htxt(4).Color = [1 1 1]*0.7;

% colors
colors = flip(cbrewer('seq', 'Reds', 4));
colors(2,:) = colors(3,:);
colors(3:7,:) = flip(cbrewer('seq', 'Blues', 5));

for ii=1:length(handles.hflux)
    handles.hflux(ii).Color = colors(ii,:);
    handles.henv(ii).Color = colors(ii,:);
    handles.hprofile(ii).Color = colors(ii,:);
    handles2.hpltsh(ii).Color = colors(ii,:);
    handles2.hplted(ii).Color = colors(ii,:);
    handles2.hpltsh(ii).LineStyle = '-';
    handles2.hplted(ii).LineStyle = '-';
end
handles2.htxt(2).Color = colors(2,:);
handles2.htxt(3).Color = colors(5,:);
uistack(handles2.hplted(3), 'top');

%L_β annotation
axes(handles.hax(2))
run = shfric.array(3);
[start,stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
[V0,L0,Lz0] = run.EddyScalesForFlux(start,stop);
betash = run.params.phys.f0 / run.bathy.hsb * run.params.bathy.sl_shelf;
Lbeta = sqrt(V0/betash);
hlbeta = liney(-Lbeta/1000); uistack(hlbeta, 'bottom');
correct_ticks('y', [], '-10');
hanno = annotation('doublearrow', [1 1]*.55, [0 0.085]+0.530);
hanno.LineWidth = 1; hanno.Color = hlbeta.Color;
hanno.HeadSize = 8; hanno.HeadStyle = 'vback3';
htext = text(430, -Lbeta/2000, 'L_\beta', 'Color', hlbeta.Color);
linkprop([htext hanno hlbeta], 'Color');
htext.Color = handles2.htxt(3).Color;

export_fig -r150 -a2 -png images/paper3/sb-flux-summary

%%
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-8352-2')
    ew = runs('../topoeddy/runew-8352-2/');
end
surfvar = 'csdye';

handles = ew.PlotCrossSBVelocity('csdye', '180');

handles.hax(1).XLim = [150 400];
handles.hax(1).YLim = [0 150];
handles.hax(3).CLim = [-1 1]*400;
handles.hax(4).CLim = [-1 1]*400;

correct_ticks('x', [], {'200'}, handles.hax(3:4));

export_fig -r150 images/cross-shelfbreak-flows-ew-8352-2.png

%% cyclone xzsections
cyc = runs('../topoeddy/runew-34-cyc/');
handles = cyc.plot_xzsection(3, '229');
handles.hrunname.delete;
linkaxes(handles.hax([1 3 4]), 'xy');
handles.hax(1).YLim = [-300 0];
handles.hax(1).XLim = [-120 80];
handles.hw.delete;
for ii=1:3
    for jj=1:4
        handles.hline(jj).hl{ii}.delete;
        handles.hline(jj).htxt{ii}.delete;
    end
end
axes(handles.hax(3));
colorbar('off');

handles.hax(2).YLim = handles.hax(1).YLim;
handles.hax(2).YTick = handles.hax(1).YTick;

correct_ticks('y', [], {'-103'}, handles.hax([1 3]));
handles.hax(2).YTickLabelMode = 'auto';
correct_ticks('y', [], {'-103'}, handles.hax(2));
handles.hax(4).YTick(4) = [];

handles.hcb(4).Ticks = sort(handles.hcb(4).Ticks * -1);
handles.hcb(4).TickLabels{2} = 'Slope Water';
handles.hcb(4).TickLabels{3} = 'Eddy Water';

export_fig -r150 -a4 -opengl images/thesis/cyclone-xzsections.png

%% ew-4341 evolution
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-4341')
    ew = runs('../topoeddy/runew-4341/');
end

opt = [];
opt.rhocontourplot = 0;
handles = ew.animate_field('rho', [], '205', 1, opt);
handles.hbathy{2}.Color = 'w';
handles.hbathy{3}.Color = 'w';
ylim([0 400]);
xlim([50 600]);
correct_ticks('y', [], {'100'});
title('Surface \rho | \lambda \approx 0.1');

export_fig -r150 -a2 images/thesis/ew-4341-waves.png

%% hovmoeller 04, 8041
ew = runArray({'ew-04', 'ew-8041'});

figure;
hax = packfig(1,2);
ew.array(1).hovmoeller('zeta', 'y', 'sb', [], hax(1));
ew.array(2).hovmoeller('zeta', 'y', 'sb', [], hax(2));
linkaxes(hax, 'xy');
ylim([0 350]);
axes(hax(1)); caxis([-1 1]*2.5e-3); colorbar('off');
axes(hax(2)); caxis([-1 1]*2.5e-3); legend('off');
hax(2).YTickLabel = {};
hax(2).YLabel.String = '';

export_fig images/ew-04-8041-sbssh-hov.png


%% avg *unmasked* sbssh at shelfbreak
ew.plot_avgProfile('zeta', 0, 'y', 'sb', 0);

%% sbssh regression
sh.print_diag('sbssh');

%% eddy + slope inflow water
ew = runArray({'ew-04', 'ew-8041', 'ew-34', 'ew-8341'});

ew.plot_ts('run.csflux.on.slope(:,1,1)');

%%
ew5 = runArray({'ew-8350-2'; 'ew-8352-2'});

%%

bvel = runArray({'ew-64461-5', 'ew-64461-8', 'ns-6362-2', 'xy-64461-8-shallow', 'ew-64361'});
for ii=1:bvel.len
    bvel.array(ii).calc_eddy_velbot;
end

bvel.plot_ts('eddy.Vb./eddy.V')

%%

N = sqrt(ew.params.phys.N2);
beta = ew.params.phys.beta;
f0 = ew.params.phys.f0;
H = ew.traj.H;
K = 1./(ew.params.eddy.dia/2);
m = 1./(N*H/(pi/2)/f0);

beta/(K^2 + m^2)*100 % cm/s ≈ km/day

%%

mm=11;
sh.array(mm).ShelfBaroclinicity;
figure;
plot(sh.array(mm).shelfbc.shelf(:,2)); hold on;
plot(sh.array(mm).shelfbc.farfield(:,2));
plot(sh.array(mm).shelfbc.sbreak.shelf(:,2));
title(sh.array(mm).name)

%% surface dye plots for gordon poster
timesteps = [1 230];

if ~exist('ew8341', 'var')
    if exist('sh', 'var')
        ew8341 = sh.array(5);
    else
        ew8341 = runs('../topoeddy/runew-8341/');
    end
end
opt = [];
opt.rhocontourplot = 0;
opt.csdcontourplot = 1;
opt.csdcontours = ew8341.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.addzeta = 1;
opt.topdown = 0;
opt.bathycolor = 'w';
handles = ew8341.mosaic_field('csdye', timesteps, opt);
handles.hcb.delete;
for ii=1:3
    handles.hfield{ii}.hzeta.LevelList = linspace(0,0.013, 4);
    handles.hfield{ii}.hzetaneg.LevelList = [-1.6e-3 -1.25e-3 -1.1e-3];
end
ylim([0 250]);
correct_ticks('y', [], {'50', '100'}, handles.hax);

handles.hsuptitle.String = ...
    ['Surface cross-shelf dye (km) | Initial Ro = 0.1 | Initial eddy ' ...
     'scales = (25 km, 400m)'];

axes(handles.hax(1));
hleg = legend([handles.hfield{1}.hcen, ...
               handles.hfield{1}.htrack, ...
               ... %handles.hfield{1}.hrho, ...
               handles.hfield{1}.hzeta], ...
              {'Eddy center', 'Track of eddy center', ...%'Eddy core',
               'SSH'}, 'Location', 'NorthWest'); %, 'FontSize', 14);
hleg.Box = 'off';

resizeImageForPub('portrait');
hleg.Position(2) = 0.50;
hleg.Position(1) = 0.13;
handles.supax.Position(4) = 0.55;

fs = [12 13 14];
for ii=1:length(handles.hax)
    axes(handles.hax(ii));
    handles.hax(ii).Color =  'none';
    handles.hfield{ii}.htlabel.FontSize = fs(1);
    handles.hfield{ii}.htlabel.FontName = 'Times';
    handles.hsuptitle.FontSize = fs(end);
    handles.hsuptitle.FontName = 'Times';
    beautify(fs, 'Times')
end
handles.hax(2).YTickLabels = {};
hleg.FontName = 'Times'
hleg.FontSize = fs(1);

export_fig -transparent -r300 -a2 -opengl images/grs-ew-8341-surface-csdye.png
