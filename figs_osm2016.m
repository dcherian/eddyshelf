%% options
fontsize = [20 22 24];

%% ew-34-surface csdye
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

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
handles = ew34.plot_xzsection(isobath, 225);
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(1));
correct_ticks('y', [], {'-50'; '-300'}, handles.hax(2));
delete(handles.hrunname);
delete(handles.hax(3:4))
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
