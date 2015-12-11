% describe fluxes
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end
factor = 2;
isobath = 4;

%% field map
opt.csdcontourplot = 1;
opt.csdcontours = ew34.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.zetaOnFirstPlot = 1;
handles = ew34.mosaic_field('csdye', [1 150 200 250 300 380], opt);
handles.hfield{1}.hzeta.LevelList = handles.hfield{1}.hzeta.LevelList(1:2:end);
handles.hcb.delete;

ylim([0 250]);
correct_ticks('y', [], {'50', '100'}, handles.hax([1 4]));

handles.supax.Position(4) = 0.70;
handles.htitle.String = 'Surface cross-shelf dye (km)';

hleg = legend(handles.hax(1), [handles.hfield{1}.hcen, ...
                    handles.hfield{1}.htrack, ...
                    handles.hfield{1}.hzeta, ...
                    handles.hfield{1}.hrho], ...
              {'Eddy center', 'Track of eddy center', 'SSH', 'Eddy core'}, ...
              'Location', 'NorthWest', 'FontSize', 14);
hleg.Box = 'off';
hleg.Position(2) = hleg.Position(2) - 0.03;
export_fig -a1 images/paper2/ew-34-surface-csdye.png

%% center tracks
set(groot,'DefaultAxesColorOrder', cbrewer('qual', 'Paired', 8));
if ~exist('ew', 'var')
    folders = {'runew-34', 'runew-8341', 'runew-5343', 'runew-5341', ...
               'runew-2340', 'runew-2341_wider'};
    names = {'Base case', 'S_{sh} = 0.05', 'r = 5e-4', 'r = 5e-3', ...
             'H_{sb} = 100m', 'H_{sb} = 300m'};
    ew = runArray(folders, names, 0, 1); % reduced version - only eddy diags
end
ew.plot_penetration;
title('');
pbaspect([1.618 1 1]);
set(gcf, 'Position', [675 175 1080 924]);
columnlegend(3, names, 'Location', 'NorthWest');
export_fig -a1 -png -pdf images/paper2/centracks

%% flux diagnostics
handles = ew34.plot_fluxts(factor, 3, 3);
handles.htitle.String = ['Flux of water above ' num2str(factor) 'H_{sb}'];
maximize; drawnow;
export_fig -a1 images/paper2/flux-diags.png

%% x-z sections
handles = ew34.plot_xzsection(isobath, 225);
correct_ticks('y', [], 3, handles.hax([1 3 4]));
correct_ticks('y', [], 4, handles.hax([2]));
drawnow;
export_fig -a1 images/paper2/ew-34-xzsection.png

%% secondary eddy
if ~exist('ew2360', 'var') | ~strcmpi(ew2360.name, 'ew-2360_wider')
    ew2360 = runs('../topoeddy/runew-2360_wider/');
    xx = [343 346 349]';
    yy = [231 237 231]';
end
opt.addvelquiver = 0;
opt.csdcontours = ew2360.csflux.x(isobath);
handles = ew2360.secondary_vortices(95, [xx yy], opt);
correct_ticks('y', [], '198', handles.hax(1));
handles.hax(1).Title.String = 'Cross shelf dye - X_{sb} (km)';
export_fig -a1 images/paper2/ew-2360-secondary-cyclone.png

%% z-slice
opt.csdcontours = ew34.csflux.x(isobath);
handles = ew34.mosaic_zslice('dye_03', 100, [225 235 245 255], opt);
for ii=1:length(handles)
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([200 400]);
ylim([0 150]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.10';
handles(1).htitle.FontSize = 22;
handles(1).htitle.FontWeight = 'normal';
handles(1).supax.Position(4) = 0.87;
axes(handles(1).hax); correct_ticks('y', '', '0');
axes(handles(3).hax); correct_ticks('y', '', '0');
export_fig -a1 images/paper2/ew-34-mosaic-zslice.png