%% options
fontsize = [18 22 26];
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

%% ew-34-surface csdye
opt.csdcontours = ew34.csflux.x([1 3 5]);
ew34.animate_field('csdye', [], '314', 1, opt);
title('Cross-shelf dye');
colorbar('delete');
xlim([100 450]);
ylim([0 150]);
beautify(fontsize)
export_fig -r150 -a2 images/osm2016/ew-34-csdsurf.png

%% xzsection
isobath = 3;
handles = ew34.plot_xzsection(3, '223');
correct_ticks('y', [], {'-50'; '-300'}, handles.hax([1 3 4]));
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(2));
linkaxes(handles.hax, 'y');
ylim([-320 0]);
delete(handles.hrunname);
export_fig('-r150' ,'-a2', '-painters', 'images/osm2016/ew-34-csdye-xz.png');
delete(handles.hax([3 4]));
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
export_fig -r150 -opengl -a2 images/osm2016/ew-34-xzsection.png

%% xzsection
ew34.PlotSingleXZSection('csdye', 3, '223');
colorbar;
xlim([-1 1]*80);

%% video
opt.dt = 2;
opt.nocolorbar = 1;
opt.addzeta = 1;
opt.stopzeta = 100;
opt.csdcontours = ew34.csflux.x([1 3 5]);
opt.vecplot = 1;
opt.csfluxplot = 2;
opt.csfluxIsobath = [1 3];
opt.AnimateZoom = 1;
opt.ZoomStart = 80;
opt.ZoomEnd = 240;
opt.ZoomXLimStart = [];
opt.ZoomXLimEnd = [100 450];
opt.ZoomYLimStart = [];
opt.ZoomYLimEnd = [0 200];

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
annoheadstyle = 'cback3';
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
opt.csdcontours = ew34.bathy.xsb + 5000;
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

% t= 210
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
handles = PlotOleanderSection(114, [11 17]);
handles.hax(2).Title.String = 'Temperature (C) | Oleander XBT | 19-20 Oct 1983';
handles.hax(2).Title.FontSize = 26;
axes(handles.hax(1));
beautify([22 24 30]);
correct_ticks('x', [], {'650'; '750'; '700'}, handles.hax(1));
axes(handles.hax(3));
beautify([22 24 30]);
export_fig -r150 -a2 images/osm2016/oleander-oct1983.png
