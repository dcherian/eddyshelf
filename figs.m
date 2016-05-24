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
sh.filter = [5 14:17];
figure;
hax = packfig(2,1);
for ff=1:length(sh.filter)
    ii = sh.filter(ff);
    run = sh.array(ii);
    axes(hax(1)); hold on;
    hplt1(ff) = plot(run.csflux.time/86400, run.csflux.off.slope(:,1,1), ...
                    'DisplayName', [run.name ' | r = ' num2str(run.params.misc.rdrg,'%.1e')]);

    axes(hax(2)); hold on;
    hplt2(ff) = plot(run.shelfbc.time/86400, run.shelfbc.shelf(:,2));
end

axes(hax(1)); legend(hplt1); beautify;
axes(hax(2)); beautify;

hplt1(1).Color = 'k';
hplt2(1).Color = 'k';
hax(1).YLabel.String = 'Shelf water offshore flux (m^3/s)';
hax(2).YLabel.String = 'BC_{0.2}';
hax(2).XLabel.String = 'Time (days)';
hax(1).YLim(1) = 0;
hax(2).YLim = [0 1];
linkaxes(hax, 'x');

insertAnnotation('figs.m');
export_fig images/bc-shfric.png