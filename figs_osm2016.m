%% options
fontsize = [20 24 26];
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

%% model schematic
opt.addzeta =1;
opt.csdcontourplot = 0;
hbathy = ew34.plot_bathy('pcolor');
hax = gca;
handles = ew34.animate_field('csdye', hax, 1, 1, opt);
handles.hfield.delete;
handles.hbathy{1}.delete;
handles.hbathy{2}.delete;
handles.hbathy{3}.delete;
handles.subax.Visible = 'off';
colormap(cbrewer('seq','Blues', 20));
caxis([-100 1800]);
hcb = colorbar;
hcb.Limits = [50 1200];
handles.htlabel.delete;
hcb.Label.String = 'Water Depth (m)';
hcb.Label.FontSize = 25;
hcb.FontSize = 18;
handles.htrack.Color = [227 54 71]/255;
export_fig -r150 -a2 -painters images/osm2016/model-setup.png

%% ew-34-surface csdye
opt.addzeta = 0
opt.rhocontourplot = 0;
opt.csdcontours = ew34.csflux.x([1 3 5]);
handles = ew34.animate_field('csdye', [], '314', 1, opt);
title('Cross-shelf dye');
colorbar('delete');
xlim([100 450]);
ylim([0 150]);
beautify([22 26 28])
xlabel('Along-shelf (km)');
ylabel('Cross-shelf (km)');
handles.htrack.Color = [227 54 71]/255;
export_fig -painters -r150 -a2 images/osm2016/ew-34-csdsurf.png

%% xzsection
opt.onlyvel = 1;
handles = ew34.plot_xzsection(3, '223', opt);
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(1));
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(2));
linkaxes(handles.hax, 'y');
ylim([-320 0]);
handles.hrunname.delete;
handles.htime.delete;
handles.hax(1).XLim = [-1 1]*80;
handles.hax(1).XLabel.String = 'Along-shelf (km)';
handles.hline(1).hl{1}.delete;
handles.hline(1).htxt{1}.delete;
handles.hline(2).hl{1}.delete;
handles.hline(2).htxt{1}.delete;
handles.hline(1).htxt{2}.Position(1) = 20;
handles.hline(2).htxt{2}.Position(1) = 50;
handles.hline(1).htxt{2}.String = 'shelfbreak depth';
handles.hline(2).htxt{2}.String = 'shelfbreak depth';
handles.hcb(1).TickLabels{1} = 'Onshore';
handles.hcb(1).TickLabels{end} = 'Offshore';
export_fig -r150 -opengl -a2 images/osm2016/ew-34-xzsection.png

%% xzsection
handles = ew34.PlotSingleXZSection('csdye', 3, '223');
hax = gca;
xlim([-1 1]*80);
hax.DataAspectRatio = [80 160 1];
beautify(fontsize);
correct_ticks('y',[], {'-50'});
colorbar('off');
xlabel('Along-shelf (km)');
ylabel('Depth (m)');
handles.htext{2}.Position(1) = 10;
handles.htext{2}.String = 'shelfbreak depth';
handles.htext{2}.Color = 'k';
handles.hline{2}.Color = 'k';
handles.htime.delete;
handles.hline{1}.delete;
handles.htext{1}.delete;
[handles.hline{3}, handles.htext{3}] = linex(0, 'eddy center', 'k');
handles.htext{3}.FontSize = fontsize(1);
handles.htext{2}.FontSize = fontsize(1);
export_fig -r150 images/osm2016/ew-34-xzsection-dye.png

%% video
opt.dt = 2;
opt.nocolorbar = 1;
opt.addzeta = 1;
opt.stopzeta = 100;
opt.csdcontours = ew34.csflux.x([1 3 5]);
opt.vecplot = 1;
opt.csfluxplot = 2;
opt.csfluxIsobath = [1];
opt.AnimateZoom = 1;
opt.ZoomStart = 80;
opt.ZoomEnd = 240;
opt.ZoomXLimStart = [];
opt.ZoomXLimEnd = [100 450];
opt.ZoomYLimStart = [];
opt.ZoomYLimEnd = [0 180];
opt.fontsize = [22 24 26];
opt.csfluxFinalize = 1;

ew34.makeVideo = 1;
ew34.animate_field('csdye', [], 1, 170, opt);

%% ew2360 xz
isobath = 3;
handles = ew2360.plot_xzsection(isobath, 225);
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(1));
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(2));
delete(handles.hrunname);
handles.hax(3).Visible = 'off';
handles.hax(4).Visible = 'off';
for ii=1:2
    handles.hline(2).htxt{ii}.Units = 'normalized';
    handles.hline(2).htxt{ii}.Position(1) = 0.3;
end
handles.hax(1).XLim = [-170 150];
handles.hax(2).XLabel.String = 'Vertical profile of offshore transport (m^2/s)';
handles.hax(2).YLabel.String = '';
handles.hax(2).YTickLabel = {};
hanno = annotation('arrow', [0.47 0.61], [1 1]*0.75);
htxt = annotation('textbox', [0.44 mean(hanno.Y) 0.2 0.2], ...
                  'String', 'Integrate horizontally', ...
                  'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                  'EdgeColor', 'none');
linkprop([hanno htxt], 'Color');
hanno.Color = [1 133 113]/255;

%% 3d schematic

annocolor = 'k'; %[1 1 1]*0.5;
annofs = 23;
annolw = 0.5;
annoheadstyle = 'none'; 'cback3';
annofontname = 'Futura Bk BT';

opt.MoveToZLevel = -850;
opt.eddthresh = 0.8;
opt.csdcontours = [];
opt.eddreducepatch = 0.3;
opt.csdreducepatch = 0.3;
opt.finalize = 1;
opt.linefilter = 1;
opt.sect = 'y';
opt.nolabels = 1;

% t= 200
ticstart = tic;
opt.x = [260, 420] * 1000;
opt.y = [300, 0] * 1000;
opt.csdcontours = []; ew34.bathy.xsb + 5000;
handles = ew34.animate_3d('200', opt);

hax = gca;
hax.DataAspectRatio = [1 1 7];
view(-130,28)

hanno = annotation(hax.Parent, 'textarrow', [0.7 0.59], [0.5 0.58], ...
                   'String', {'There is a';  'persistent bulge'; ...
                    'in the eddy'; 'below'; 'shelfbreak depth'}, ...
                   'LineWidth', annolw, 'Color', annocolor, ...
                   'HeadStyle', annoheadstyle, 'TextMargin', 0.05, ...
                   'FontName', annofontname, 'FontSize', annofs, ...
                   'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
drawnow;
toc(ticstart);
export_fig -nocrop -r120 -a4 images/osm2016/ew-34-3d-200-nocsd.png

% t= 215
opt.x = [260, 420] * 1000;
opt.y = [300, 0] * 1000;
opt.csdcontours = [];
handles = ew34.animate_3d('215', opt);

handles.hlight.Position(2) = 700;
hax = gca;
hax.DataAspectRatio = [1 1 7];
view(-130,28)

hanno = annotation(hax.Parent, 'textarrow', [0.72 0.62], [0.81 0.65], ...
                   'String', {'There is a';  'persistent bulge'; ...
                    'in the eddy'; 'below'; 'shelfbreak depth'}, ...
                   'LineWidth', annolw, 'Color', annocolor, ...
                   'HeadStyle', annoheadstyle, 'TextMargin', 0.05, ...
                   'FontName', annofontname, 'FontSize', annofs, ...
                   'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

hanno = annotation(hax.Parent, 'textarrow', [0.74 0.59], [0.55 0.63], ...
                   'String', {'Cyclonic';  'unstable wave'}, ...
                   'LineWidth', annolw, 'Color', annocolor, ...
                   'HeadStyle', annoheadstyle, 'TextMargin', 0.05, ...
                   'FontName', annofontname, 'FontSize', annofs, ...
                   'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

export_fig -nocrop -r120 -a4 images/osm2016/ew-34-3d-215-nocsd.png

% t= 220
opt.x = [240, 410] * 1000;
opt.y = [300, 0] * 1000;
opt.csdcontours = ew34.bathy.xsb+5000;
handles = ew34.animate_3d('220', opt);

hax = gca;
hax.DataAspectRatio = [1 1 7];
view(-130,28)

hanno = annotation(hax.Parent, 'textarrow', [0.7 0.59], [0.5 0.58], ...
                   'String', {'There is a';  'persistent bulge'; ...
                    'in the eddy'; 'below'; 'shelfbreak depth'}, ...
                   'LineWidth', annolw, 'Color', annocolor, ...
                   'HeadStyle', annoheadstyle, 'TextMargin', 0.05, ...
                   'FontName', annofontname, 'FontSize', annofs, ...
                   'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

hanno = annotation(hax.Parent, 'textarrow', [0.7 0.58], [0.8 0.66], ...
                   'String', {'The cyclonic wave'; 'moves around'; 'the eddy.'}, ...
                   'LineWidth', annolw, 'Color', annocolor, ...
                   'HeadStyle', annoheadstyle, 'TextMargin', 0.05, ...
                   'FontName', annofontname, 'FontSize', annofs, ...
                   'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

handles.hcsd.Visible = 'off';
handles.hcsdsurf{1}.Visible = 'off';
handles.hsect.Visible = 'off';
handles.hsectoutline.Visible = 'off';
handles.hplane.Visible = 'off';
export_fig -nocrop -r120 -a4 images/osm2016/ew-34-3d-220-nocsd.png

hanno.String = {'Shelf water is advected'; 'over the bulge'; 'and the cyclonic wave'};
handles.hcsd.Visible = 'on';
handles.hcsdsurf{1}.Visible = 'on';
handles.hsect.Visible = 'on';
handles.hsectoutline.Visible = 'on';
handles.hplane.Visible = 'on';
export_fig -nocrop -r120 -a4 images/osm2016/ew-34-3d-220-csd.png

%% churchill (1986) section
% Oct 19-20, 1983
handles = PlotOleanderSection(114, [17 11]);
handles.hax(2).Title.String = 'Temperature (C) | Oleander XBT | 19-20 Oct 1983';
axes(handles.hax(1));
handles.hax(1).YLabel.String = 'Depth (m)';
correct_ticks('x', [], {'650'; '750'; '700'}, handles.hax(1));
axes(handles.hax(3));
handles.hax(3).YTickLabels = {};
handles.hax(3).XAxis.TickLabelGapMultiplier = 0;
handles.hax(3).XAxis.TickLabelGapOffset = 0;    
handles.hleg.delete;
colors = get(groot, 'DefaultAxesColorOrder');
handles.hprofile(2).Color = colors(6,:);
for ii=1:2
    handles.xlines{ii}.Color = handles.hprofile(ii).Color;
    handles.xlines{ii}.LineWidth = 2;
end
export_fig -r150 -a2 images/osm2016/oleander-oct1983.png
