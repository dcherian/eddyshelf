%% load runs
cyc = runs('../topoeddy/runew-34-cyc/');
ew = runArray({'runew-34', 'runew-36'});
ew34 = ew.array(1);
ew36 = ew.array(2);
ew2360 = runs('../topoeddy/runew-2360_wider/');

%% surface dye field - anticyclone and cyclone
figure;
ax(1) = subplot(211);
ew34.animate_field('dye_03', ax(1), 320, 1);
title('Anticyclone moving southwards');
ax(2) = subplot(212);
cyc.animate_field('dye_03', ax(2), 290, 1);
title('Cyclone moving northwards');
export_fig images/ew-34-cyc-acyc-eddye.png

%% eddye y-z section
tind = 320; limx = [0 150];
ew36.plot_yzsection(tind)
ylim([-400 0]);
xlim(limx);
hax = gca;
hax.YTick = sort([hax.YTick -floor(ew36.bathy.hsb)]);
suplabel('', 't'); title('Eddy dye');
text(0.85*limx(2), -1 * ew36.eddy.Lgauss(tind), ...
         {'vertical','scale'}, 'VerticalAlignment', 'Bottom', ...
     'HorizontalAlignment','Center');
% set(gcf, 'Position', [515 262 650 537]);
export_fig images/ew-36-eddye-yz-320.png

%% mosaics
%% ew-34 zslice mosaic
handles = ew34.mosaic_zslice('dye_03', 200, [220 230 240 248]);
for ii=1:length(handles)
    handles(ii).csdsurf.LevelList = [60] * 1e3;
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off'
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([200 360]);
ylim([30 130]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.1';
handles(1).htitle.FontWeight = 'normal';
handles(1).htitle.FontSize = 22;
handles(1).supax.Position(4) = 0.88;

export_fig images/comm-6/ew-34-mosaic-zslice-1.png

%% ew-34 zslice mosaic 2
handles = ew34.mosaic_zslice('dye_03', 200, [300 310 320 330]);
for ii=1:length(handles)
    handles(ii).csdsurf.LevelList = [60] * 1e3;
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([120 320]);
ylim([30 130]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.1';
handles(1).htitle.FontSize = 22;
handles(1).htitle.FontWeight = 'normal';
handles(1).supax.Position(4) = 0.85;

export_fig images/comm-6/ew-34-mosaic-zslice-2.png

%% ew- 36 zslice mosaic 1
handles = ew.array(2).mosaic_zslice('dye_03', 200, [130 140 150 160]);
for ii=1:length(handles)
    handles(ii).csdsurf.LevelList = [60] * 1e3;
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([250 400]);
ylim([30 130]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.25';
handles(1).htitle.FontWeight = 'normal';
handles(1).htitle.FontSize = 22;
handles(1).supax.Position(4) = 0.89;

export_fig images/comm-6/ew-36-mosaic-zslice-1.png

%% ew- 36 zslice mosaic 2
handles = ew.array(2).mosaic_zslice('dye_03', 200, [160 175 200 220]);
for ii=1:length(handles)
    handles(ii).csdsurf.LevelList = [60] * 1e3;
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([200 400]);
ylim([30 130]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.25';
handles(1).htitle.FontWeight = 'normal';
handles(1).htitle.FontSize = 22;
handles(1).supax.Position(4) = 0.85;

export_fig images/comm-6/ew-36-mosaic-zslice-2.png

%% ew-2360 z-slice mosaic
handles = ew2360.mosaic_zslice('dye_03', 200, [63 70 77 95]);
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
handles(1).supax.Position(4) = 0.87;
axes(handles(1).hax); correct_ticks('y', '', {'150' '200'});
axes(handles(3).hax); correct_ticks('y', '', {'150' '200'});

export_fig images/ew-2360-mosaic-zslice.png

%% v, density cross-section
ew34.plot_fluxes(4, 225);
subplot(221); correct_ticks('y', [], '-50');
subplot(223); correct_ticks('y', [], '-50');
subplot(224); correct_ticks('y', [], '-50');
xlim([-2 0.3]);
ylim([-320 0]);
export_fig images/comm-6/ew-34-flux-sections-223.png

% cyclone
%
cyc.plot_fluxes(4, 225);
subplot(221); correct_ticks('y', '', '-100');
subplot(222); correct_ticks('y', '', '-100');
subplot(223); correct_ticks('y', '', '-100');
subplot(224); correct_ticks('y', '', '-100')
xlim([-2.3 0.85]);
ylim([-330 0]);
export_fig images/comm-6/ew-34-cyc-flux-sections-225.png