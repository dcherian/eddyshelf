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
opt.csdcontours = ew34.csflux.x([1 3 5]);
ew34.makeVideo = 1;
ew34.animate_field('csdye', [], 1, [], opt);

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
opt.MoveToZLevel = -850;
opt.eddthresh = 0.8;
opt.csdcontours = ew34.bathy.xsb+5000;
opt.eddreducepatch = 0.3;
opt.csdreducepatch = 0.3;
opt.finalize = 1;
opt.linefilter = 1;
opt.sect = 'y';
opt.x = [200, 410] * 1000;
opt.y = [300, 0] * 1000;

handles = ew34.animate_3d('220', opt);
handles.hcsd.FaceAlpha = 0.8;
handles.hcsd.Visible = 'off';
handles.hcsd.FaceColor = [107 174 214]/255;
handles.hedd.FaceColor = brighten([215 48 31]/255, 0.1);

hax = gca;
hax.Title.String = '';
hax.XAxis.Visible = 'off';
hax.YAxis.Visible = 'off';
hax.ZAxis.Visible = 'off';
hax.Box = 'off';
hax.DataAspectRatio = [1 1 5];
zlim([-850 1.1]);
view(-130,28)


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
