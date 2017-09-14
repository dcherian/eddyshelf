%% shelf slopes
folders = { ...
    'ew-04', 'ew-8040', 'ew-8041', 'ew-8042', ...
    ... %'ew-8150', 'ew-8151', ...
    'ew-82340', 'ew-82342', 'ew-82343', ... %'ew-82344' ...
    'ew-34', 'ew-8341', ...%'ew-8361', ...
    'ew-8350-2', 'ew-8351-2', 'ew-8352-2', ...
    'ew-8352', 'ew-8342-2', ...
    'ew-8380', 'ew-8381', 'ew-8384', 'ew-8385', ...
    'ew-8383', 'ew-8392'... %, 'ew-8346', ...
    ... %'ew-583411', 'ew-583413', ...
    ... %'ew-583414', 'ew-583415', ...
    ...    % 'ew-36-sym', 'ew-8361'
          };
sh = runArray(folders);
fontSizes = [12.5 13 13];

% Notes
% -----
% 1. 8042 is very different from 8040, 8041!!!
% 2. 8361 is bad because it seems to be undergoing a different
%    type of instability. more like
%    the really deep eddies. things are weird!

for ii=14:15
    if isempty(strfind(sh.array(ii).supply.hash, 'af8c112'))
        disp(sh.array(ii).name)
        sh.array(ii).avgSupplyJet;
    end
    % sh.array(ii).radius = sh.array(ii).csflux.radius;
    % sh.array(ii).ShelfbreakSSHScale;
    %    sh.array(ii).csfluxes;
    %    sh.array(ii).avgStreamerVelSection(1);
    %sh.array(ii).VolumeBudget;
    %    sh.array(ii).ShelfBaroclinicity;
end

%% inflow/outflow vertical profiles

sh.filter = [1:sh.len];
handles = sh.PlotFluxVertProfiles(fontSizes);
hfig = gcf; hfig.Resize = 'off';
hfig.Position(2) = 50;
hfig.Position(4) = 800;
resizeImageForPub('portrait');
pos = handles.hcbar(1).Position;
hanno(1) = annotation('line', ...
                      [1 1]*pos(1), pos(2) + [0 pos(4)], ...
                      'Color', 'k', 'Units', 'normalized');
pos = handles.hcbar(2).Position;
hanno(2) = annotation('line', ...
                      [1 1]* (pos(1) + pos(3)/(3-1)*(1.5-1)), pos(2) + [0 pos(4)], ...
                      'Color', 'k', 'Units', 'normalized');

export_fig -c[Inf,0,Inf,0] images/paper3/sb-vert-profiles.pdf

%% fluxes

cmds = 'no_name_points; small_points;';

figure; hax(1) = gca;
indices{1} = [14:17]; % ew-838x
indices{2} = [9:11]; % ew-835x-2
indices{3} = [7 8]; % ew-834x
indices{4} = [1 2 3]; % ew-804x
jitter = linspace(-0.04, 0.04, 4);
% colors ={[1 1 1]*0.65;... %[51, 102, 170]/255;
%          [51, 102, 170]/255;
%          [238, 51, 51]/255;
%          [0 0 0]};
markers = ['o', '.', '^', 'x']
for ff=[3 2 4 1]
    sh.filter = indices{ff};
    sh.print_diag('flux slope', [1, markers(ff), jitter(ff)], hax(1), cmds);
    disp(' ')
end
title('');
xlim([0 1.1])
ylim([0 1.3])
line45;
axis square;
beautify(fontSizes, 'Times')

resizeImageForPub('onecolumn');

% hax(2) = subplot(122);
% warning off;
% sh.filter = [];
% sh.print_diag('avg flux', [], hax(2), cmds);
% title('')
% beautify(fontSizes, 'Times')
% pbaspect([1 1 1])

export_fig -c[Inf,0,Inf,0] images/paper3/sh-avgflux-slope.pdf

%% eddy water on shelf scale

figure;
hax(1) = subplot(121);
hax(2) = subplot(122);

[diags, plotx, err, norm, color, rmse, P, Perr, h1] = ...
    sh.print_diag('supply', [], hax(1), 'no_name_points, small_points');
[diags, plotx, err, norm, color, rmse, P, Perr, h2] = ...
    sh.print_diag('eddyonshelf', [], hax(2), 'no_name_points, small_points');

hax(1).Title.String = '(a)';
hax(2).Title.String = '(b)';

hax(1).XLim(1) = 0;
hax(2).YLim(1) = 0;

resizeImageForPub('portrait')
axes(hax(1)); beautify(fontSizes, 'times')
axes(hax(2)); beautify(fontSizes, 'times')

hax(1).XLim(end) = 35;
hax(1).YLim(end) = 45;
hax(2).XLim(end) = 25;
h1.htext.FontSize = 11;
h2.htext.FontSize = 11;

export_fig -c[Inf,0,Inf,0] images/paper3/parameterizations.pdf

%% surface dye panels
timesteps = [1 100 150 200 300 350];

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
opt.bathycolor = 'w';
handles = ew8341.mosaic_field('csdye', timesteps, opt);
handles.hcb.delete;
for ii=1:3
    handles.hfield{ii}.hzeta.LevelList = linspace(0,0.013, 4);
    handles.hfield{ii}.hzetaneg.LevelList = [-1.6e-3 -1.25e-3 -1.1e-3];
end
for ii=4:6
    handles.hfield{ii}.hzeta.delete;
    handles.hfield{ii}.hzetaneg.delete;
end
ylim([0 250]);
correct_ticks('y', [], {'50', '100'}, handles.hax([1 4]));
for ii=1:6
    axes(handles.hax(ii));
    beautify([12 13 14], 'Times')
    handles.hfield{ii}.hcen.MarkerSize = 8;
    handles.hfield{ii}.hbathy{1}.ShowText = 'off';
end

hfig = gcf;
hfig.Position(4) = 450;
resizeImageForPub('portrait')
moveSubplotsCloserInY(2, 3, handles.hax, -1/8-1/16)

handles.supax.Position(4) = 0.71;
handles.hsuptitle.String = ...
    ['Surface cross-shelf dye (km) | Initial Ro = 0.1 | Initial eddy ' ...
     'scales = (25 km, 400m)'];
handles.hsuptitle.FontName = 'Times';
handles.hsuptitle.FontSize = 12;

% hleg.delete;
% axes(handles.hax(1));
% [hleg,icons] = legend([handles.hfield{1}.hcen, ...
%                     handles.hfield{1}.htrack, ...
%                     handles.hfield{1}.hzeta], ...
%                       {'Eddy center', 'Track of eddy center', ...
%                     'SSH'}, 'Location', 'NorthWest', 'FontName', 'Times');
% hleg.Box = 'off';
% hleg.Position(2) = 0.65;
% hleg.FontName = 'Times'
% icons(end).Children.Children(1).LineWidth = 1;
% icons(end).Children.Children(2).LineWidth = 1;
% icons(end).Children.Children(3).LineWidth = 1;
% labels = findobj(icons, 'type', 'text');
% for ii=1:length(labels)
%     labels(ii).FontName = 'Times';
%     labels(ii).FontSize = 8;
% end

% axes(handles.hax(5))
% hline = linex(350, [], 'k');

export_fig -r300 -a2 -opengl images/paper3/ew-8341-surface-csdye.png

%% shelf baroclinicity
thresh = 2;

figure;
hax(1) = subplot(121);
set(gca, 'Color', 'none')
sh.print_diag('shelfbc', thresh, hax(1), 'no_name_points');
title('Supply jet on shelf');
ylabel('BC');
axis square;
xlim([0 1]);
ylim([0 1]);
linex(0.3);
for ii=1:length(hax(1).Children)
    hax(1).Children(ii).MarkerSize = 12;
end
beautify(fontSizes, 'Times');

hax(2) = subplot(122);
set(gca, 'Color', 'none')
sh.print_diag('sbreakbc', thresh, hax(2), 'no_name_points');
title('Outflow at shelfbreak');
ylabel('BC');
axis square;
xlim([0 1]);
ylim([0 1]);
linex(0.3);
for ii=1:length(hax(2).Children)
    hax(2).Children(ii).MarkerSize = 12;
end
beautify(fontSizes, 'Times');

set(hax(1), 'TickLength', [1 1]*0.02)
set(hax(2), 'TickLength', [1 1]*0.02)
resizeImageForPub('portrait');

export_fig -c[Inf,0,Inf,0] images/paper3/shelfbc.pdf

%% shelfbc time series
% phi = sh.print_params('bathy.hsb./(V0./bathy.S_sh/sqrt(phys.N2))');
% sh.sort(phi);
% phi = sh.print_params('bathy.hsb./(V0./bathy.S_sh/sqrt(phys.N2))');

% figure; maximize; hold on;
% corder_backup = sh.sorted_colors;
% for ii=1:sh.len
%     run = sh.array(ii);
%     handles(ii) = plot(smooth(run.shelfbc.shelf,10), ...
%                        'DisplayName', ['\phi = ' num2str(phi(ii), '%.2f')]);
%     if run.params.misc.rdrg ~= 0
%         handles(ii).LineStyle = '--';
%         handles(ii).Color = 'k';
%     end
% end
    % sh.reset_colors(corder_backup);

%% volume budget
handles = ew8341.VolumeBudget;
[tstart,tstop] = ew8341.flux_tindices(ew8341.csflux.off.slope(:,1,1));
xlim([ew8341.time([tstart-120])/86400 450]);
handles.hax.XLabel.Position(1) = 450;
t = ew8341.time([tstart; tstop])/86400
[hlines, htimes] = linex(t, {'t_{start}'; 't_{stop}'});
handles.hax.XTickMode = 'auto'

for ii=1:2
    hlines(ii).YData = [-0.3e4 0.95e4];
    hlines(ii).LineStyle = '--'
    htimes(ii).Rotation = 0;
    htimes(ii).Position(2) = 0.8e4;
    htimes(ii).Position(1) = htimes(ii).Position(1) - 5;
    htimes(ii).HorizontalAlignment = 'right';
    handles.htxt(ii).Position(1) = 0.8;
end

delete(handles.residual);
delete(handles.lowsponge);
legend('-dynamiclegend')
handles.hisponge.DisplayName = 'Along-shelf at eastern boundary';
handles.hispongesh.DisplayName = 'Along-shelf at eastern boundary: shelf water';
handles.lowsponge.DisplayName = 'Along-shelf at western boundary';
title('Volume budget for the shelf');
resizeImageForPub('portrait')
pbaspect([1.6 1 1])
beautify([12 13 14], 'Times')

export_fig -transparent images/paper3/ew-8341-volume-budget.pdf

%% flat and sloping, sbssh
if ~exist('ew04', 'var')
    ew04 = runArray({'ew-04', 'ew-8041'});
    ew04 = runArray({'ew-8350-2', 'ew-8352-2'});
end

ew04.name = {'Flat shelf'; 'Sloping shelf'};
figure;
% hax(1) = subplot(121); hold on;
% ew34.plot_avgProfile('zeta', 1, 'y', 'sb', 1, hax(1));
% title('Mean SSH at shelfbreak (normalized)');
% pbaspect([1.617 1 1]);
% hax(2) = subplot(122); hold on;
cla('reset');
hax = gca;
ew04.plot_avgProfile('zeta', 0, 'y', 'sb', 1, hax, 'center');
ylabel('');
title('   Mean shelfbreak SSH (m)');
beautify([12 12.5 13]);
pbaspect([4 1 1]);
resizeImageForPub('portrait');

export_fig -a1 -r600 -c[Inf,0,Inf,0] images/paper3/sb-ssh-ew04-flat-sloping.png

%% compare flat and sloping
if ~exist('ew04', 'var')
    % ew04 = runArray({'ew-04', 'ew-8041'});
    ew04 = runArray({'ew-8350-2', 'ew-8352-2'});
    % ew04 = runArray({'ew-8380', 'ew-8381'});
end

time = 220;
opt.drawtrack = 0;
opt.drawcen = 0;
opt.addvelquiver = 1;
opt.rhocontourplot = 0;
opt.csdcontourplot = 0;
opt.dxi = 8; opt.dyi = 5;

ew04.array(1).read_velsurf(find_approx(ew04.array(1).time, time*86400, 1), 1);
ew04.array(1).read_csdsurf(find_approx(ew04.array(1).time, time*86400, 1), 1);

ew04.array(2).read_velsurf(find_approx(ew04.array(2).time, time*86400, 1), 1);
ew04.array(2).read_csdsurf(find_approx(ew04.array(2).time, time*86400, 1), 1);

clear handles
figure; maximize;
hax = packfig(2,2);
for ii=1:2
    run = ew04.array(ii);

    handles(ii) = run.animate_field('csdye', hax(ii), num2str(time), 1, opt);
    xlim([170 400]);
    ylim([0 120]);

    handles(ii).hbathy{1}.ShowText = 'off';
    handles(ii).hquiv.Color = [1 1 1]*0;
    handles(ii).hquiv.LineWidth = 1;
    handles(ii).hquiv.AutoScaleFactor = 3.5;
    handles(ii).htlabel.Position(2) = 0.1;
    handles(ii).hbathy{1}.Color = [1 1 1]*0;
    handles(ii).hbathy{2}.Color = [1 1 1]*0;
    handles(ii).hcb.delete;
    hax(ii).XTickLabel = {};
    hax(ii).XLabel.delete;
    caxis([-50 250])

    axes(hax(ii+2))
    [hyy(ii,:), hplt(ii,1), hplt(ii,2)] = ...
        plotyy(run.rgrid.x_rho(1,:)/1000, run.vsurf(:, run.bathy.isb, tind), ...
               run.rgrid.x_rho(1,:)/1000, ...
               (run.csdsurf(:, run.bathy.isb, tind) - run.bathy.xsb)/1000);
    hplt(ii,1).Color = 'k';
    hold all;

    dx = 1;
    uavg = avg1(run.usurf(:, run.bathy.isb, tind));
    lowerquiv(ii) = ...
        quiver(run.rgrid.x_rho(1,2:dx:end-1)'/1000, zeros(size(run.rgrid.x_rho(1,2:dx:end-1)')), ...
               zeros(size(run.rgrid.x_rho(1,2:dx:end-1)')), ...
               run.vsurf(2:dx:end-1, run.bathy.isb, tind), 0, 'Color', 'k');

    hyy(ii,1).XLim = [240 400];
    hyy(ii,2).XLim = [240 400];

    hyy(ii,1).Box = 'off';
    hyy(ii,2).Box = 'off';

    axes(hyy(ii,1));
    hyy(ii,1).Color = 'none';
    liney(0);
end

hyy(1,1).YLim = [-1 1]*0.05;
hyy(2,1).YLim = [-1 1]*0.05;
hyy(1,2).YLim = [-40 180];
hyy(2,2).YLim = [-40 180];

hyy(1,1).YTick = [-0.05 -0.03 0 0.03 0.05];
linkprop(hyy(:,1), 'YTick');
hyy(1,2).YAxis.Visible = 'off';
hyy(2,1).YAxis.Visible = 'off';
hyy(2,1).Box = 'off';
hyy(2,1).YTickLabel = {};
hyy(1,2).YTickLabel = {};

hax(2).YLabel.delete;
hax(2).YTickLabel = {};
hax(1).Title.String = 'a) Flat shelf | S_{sh} = 0';
hax(2).Title.String = 'b) Sloping shelf | S_{sh} = 0.05';
correct_ticks('y', [], {'0'}, hax(1));

handles(1).supax.Position(4) = 0.73;
handles(1).supax.Title.String = 'Surface cross-shelf dye (km)';
handles(1).supax.Title.FontSize = 13;
handles(1).supax.Title.FontName = 'Times';

hax(3).YAxis.Color = 'k';
hax(3).YLabel.String = {'Cross-isobath'; 'surface velocity (m/s)'};
hax(3).YLabel.Color = hplt(1,1).Color;
hax(3).XLabel.String = 'X (km)';
hyy(2,2).YLabel.String = 'Cross-shelf dye - Y_{sb}';
hyy(2,2).YLabel.Color = hplt(1,2).Color;
hax(4).XLabel.String = 'X (km)';

resizeImageForPub('portrait');

axes(hyy(1,1))
htxt(1) = text(250, 0.02, 'offshore', 'Color', hplt(ii,1).Color);
htxt(2) = text(360, -0.0225, 'onshore', 'Color', hplt(ii,1).Color);

axes(hyy(2,2));
htxt(3) = text(340, 150, 'eddy water', 'Color', hplt(ii,2).Color);
htxt(4) = text(250, 20, 'shelf water', 'Color', hplt(ii,2).Color);

axes(hyy(1,1)); htxt(5) = text(0.05, 0.9, 'c)', 'Color', 'k', 'units', 'normalized', 'FontSize', 11);
axes(hyy(2,2)); htxt(6) = text(0.05, 0.9, 'd)', 'Color', 'k', 'units', 'normalized', 'FontSize', 11);

for ii=1:4
    axes(hax(ii)); beautify([12.5 13 14], 'Times');
    htxt(ii).FontSize = 12;
end
for ii=1:2
    handles(ii).htlabel.FontSize = 13;
end
axes(hyy(2,2)); beautify([12.5 12.5 13], 'Times');

for ii=1:2
    for jj=1:2
        axes(hyy(ii,jj));
        pbaspect([2 1 1]);

        hyy(ii,jj).Position(2) = 0.35;
    end
end

axes(hax(2))
hl = liney(37.5-12*1.37);
hl.Color = [1 1 1]*0;
hannotxt = text(max(xlim), hl.YData(1), ' 1.37L_\beta', 'FontSize', ...
                11, 'FontName', 'Times');

export_fig -a1 -c[Inf,0,Inf,0] -r600  images/paper3/sbsnapshot-flat-sloping.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flux summary
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
        '  0.05   3e-3'; ...
            };
    shfric = runArray(folders, names);
end
handles = shfric.PlotFluxSummary2(1);

% colors
colors = flip(cbrewer('seq', 'Reds', 4));
colors(2,:) = colors(3,:);
colors(3:7,:) = flip(cbrewer('seq', 'Blues', 5));

for ii=1:length(handles.hflux)
    handles.hflux(ii).Color = colors(ii,:);
    handles.henv(ii).Color = colors(ii,:);
    handles.hprofile(ii).Color = colors(ii,:);
end

% handles.hax(1).Position(2) = 0.515;
% handles.hax(3).Position(1) = 0.72
% handles.hax(3).Position(4) = 0.75;
handles.hax(1).YLim = [0 15];

resizeImageForPub('portrait');

axes(handles.hax(1));
delete(handles.hleg);
handles.hleg(1) = legend(handles.hflux(1:2), shfric.name{1:2});
handles.hleg(2) = legend(handles.henv(3:5), shfric.name{3:5});
handles.hleg(1).Box = 'off';
handles.hleg(2).Box = 'off';
handles.hleg(1).Position(1) = 0.63;
handles.hleg(1).Position(2) = 0.78;
pos1 = handles.hleg(1).Position;
handles.hleg(2).Position(1) = handles.hleg(1).Position(1) + pos1(3) + 0.02;
handles.hleg(2).Position(2) = pos1(2) + pos1(4) - handles.hleg(2).Position(4);
pos2 = handles.hleg(2).Position;

handles.hax(3).Position(1) = 0.74;
handles.hfig.Position(4) = 600;
handles.hax(1).PlotBoxAspectRatio = [2 1 1];
handles.hax(2).PlotBoxAspectRatio = [2 1 1];
handles.hax(2).Position(2) = 0.2;
handles.hax(3).Position(2) = 0.25;

tblx = [0.63 0.98];
dy = 0.05;
htable.delete;
htable(1) = text(1.16, 0.885+dy, 'S_{sh}   r_f', 'Units', 'normalized', ...
                 'FontSize', 12, 'FontName', 'Times');
htable(2) = text(1.53, 0.885+dy, 'S_{sh}     r_f', 'Units', 'normalized', ...
                 'FontName', 'Times', 'FontSize', 12);
htable(3) = annotation('line', tblx, [1 1]*0.84+dy, 'LineWidth', 0.5);
htable(4) = annotation('line', tblx, [1 1]*0.79+dy, 'LineWidth', 0.5);
htable(5) = annotation('line', tblx, [1 1]*0.70+dy, 'LineWidth', 0.5);
linkprop(htable, 'Color');
htable(1).Color = handles.hax(1).XAxis.Color;
linkprop(handles.hleg, 'TextColor');
handles.hleg(1).TextColor = htable(1).Color;

handles.hax(3).XLabel.Position(2) = 0.14;
handles.hax(3).XAxis.Exponent = 0;

labels = findobj(handles.legobj, 'type', 'text');
for ii=1:length(labels)
    labels(ii).FontName = 'Times';
end

set(gcf, 'renderer', 'painters');
export_fig images/paper3/sb-flux-summary.pdf

%%
sh.print_diag('params table', 'images/paper3/sh-params-table.org')

%% friction mosaic
if ~exist('shfric2', 'var')
    folders = { ...
        'runew-8341', ...
        'runew-583413', 'runew-583411', ...
        'runew-583414', 'runew-583415', ...
              };
    shfric2 = runArray(folders);
end

for ii=1:5
    shfric2.array(ii).CalcUbarCrossIsobathScale;
    %     shfric2.array(ii).usurf = [];
    %     shfric2.array(ii).vsurf = [];
    %     shfric2.array(ii).read_velsurf;
end

opt.addvelquiver = 1;
opt.quivercolor = [0 150 136]/255;
opt.scale = 5;
opt.limy = [0 150];
shfric2.filter = [1 2 3 4];
handles = shfric2.mosaic_field('csdye', [350 370 280 280], opt, [-40 150]);

%% friction mosaic
sh.filter = [1 2 5:13];
%sh.plot_avgProfile('vbar');
figs = [0 0 0 0 1];
sh.plot_fluxes(1,1,1,figs);

%% friction plot
if ~exist('shfric2', 'var')
    folders = { ...
        'runew-8341', ...
        'runew-583413', 'runew-583411', ...
        'runew-583414', 'runew-583415', ...
              };
    shfric2 = runArray(folders);
end

shfric2.sort(shfric2.print_params('run.params.misc.rdrg'));
rdrg = shfric2.print_params('run.params.misc.rdrg');
for ii=1:shfric2.len
    shfric2.name{ii} = num2str(rdrg(ii), '%.0e');
    shfric2.name{ii}(4) = [];
end
shfric2.name{1} = '0';
shfric2.array(1).read_zeta;
shfric2.array(1).read_velsurf;

shfric2.array(end).read_zeta;
shfric2.array(end).read_velsurf;

opt.limy = [0 150];
opt.limx = [200 480];
opt.quivercolor = [217 33 32]/255;
opt.drawtrack = 1;
%[136 34 85]/255; [0 150 136]/255; shfric2.array(1).shelfSlopeColor('dark'); %[1 0 0];
opt.addvelquiver = 1;
opt.addzeta = 0;
opt.rhocontourplot = 0;
opt.csdcontourplot = 0;
opt.addzeta = 1;
% opt.commands = 'caxis([-1 1]*1e-3)';
opt.dyi = 8;
opt.dxi = 12;
opt.scale = 2;
opt.uref = 5e-4; opt.vref = opt.uref;
str = 'acbd';
trackcolor = 'k';
zetacolor = [68 70 153]/255; [136 204 238]/255;[244 67 54]/255; [1 1 1]*0.6;
var = 'zeta';
figure;
clf; clear handles;
hax = packfig(4, 2);
for ii=1:2
    if ii == 1
        run = shfric2.array(1);
    else
        run = shfric2.array(end);
    end
    axes(hax(ii)); %hax(ii) = subplot(2,3,ii);
    handles(ii,1) = run.animate_field(var, hax(ii),  '190', 1, opt);
    handles(ii,1).hcb.delete;
    if ii==1
        title(['r_f = 0']);
    else
        title(['r_f = ' num2str(run.params.misc.rdrg, '%1.1e')]);
    end
    % handles(ii,1).htlabel.Position(2) = 0.14
    handles(ii,1).htlabel.String = [str(ii) ') ' handles(ii,1).htlabel.String];
    handles(ii,1).htrack.Color = trackcolor;
    xlabel('');
    % handles(ii,1).hfield.Visible = 'off';

    U = handles(ii,1).hquiv.UData;
    V = handles(ii,1).hquiv.VData;
    y = handles(ii,1).hquiv.YData;
    speed = hypot(U,V);
    %nanmask = speed > 0.175*max(speed(:)) | y > 80;
    nanmask = y > 38;
    handles(ii,1).hquiv.UData(nanmask) = NaN;
    handles(ii,1).hquiv.VData(nanmask) = NaN;

    uistack(handles(ii,1).hquiv, 'top');
    beautify([12 13 14], 'Times');

    axes(hax(ii+2)); % = subplot(2,3,ii+3);
    handles(ii,2) = run.animate_field(var, hax(ii+2),  '230', 1, opt);
    title('');
    handles(ii,2).hcb.delete;
    %handles(ii,2).htlabel.Position(2) = 0.14
    handles(ii,2).htlabel.String = [str(ii+2) ') ' handles(ii,2).htlabel.String];
    handles(ii,2).htrack.Color = trackcolor;

    U = handles(ii,2).hquiv.UData;
    V = handles(ii,2).hquiv.VData;
    speed = hypot(U,V);
    %    nanmask = speed > 0.15*max(speed(:)) | y > 80;
    handles(ii,2).hquiv.UData(nanmask) = NaN;
    handles(ii,2).hquiv.VData(nanmask) = NaN;
    % handles(ii,2).hfield.Visible = 'off';

    zpos = [2 4 8 12]*1e-3;
    zneg = [-7 -5 -3 -0.75 -0.5 -0.25 0]*1e-3;

    handles(ii,1).hzeta.Color = zetacolor;
    handles(ii,1).hzetaneg.Color = zetacolor;
    handles(ii,2).hzeta.Color = zetacolor;
    handles(ii,2).hzetaneg.Color = zetacolor;
    handles(ii,1).hzeta.LevelList = zpos;
    handles(ii,1).hzetaneg.LevelList = zneg;
    handles(ii,2).hzeta.LevelList = zpos;
    handles(ii,2).hzetaneg.LevelList = zneg;

    uistack(handles(ii,2).hquiv, 'top');
    beautify([12 13 14], 'Times');
end
hax(2).YLabel.String = '';
hax(2).YTickLabel = {};
hax(4).YLabel.String = '';
hax(4).YTickLabel = {};
hax(1).XTickLabel = {};
hax(2).XTickLabel = {};

resizeImageForPub('portrait');

hax(5) = subplot(4, 2, [5 6]);
handles3 = shfric2.plot_ts('-run.ubarscale.scale/1000', hax(5), ...
                           'run.ubarscale.time/86400');
ylabel({'Distance';  'from shelfbreak (km)'});
xlabel('');
delete(handles3.htind);
delete(handles3.hmaxloc);
hleg = legend;
hleg.Location = 'SouthEast';
htxt(1) = text(0.04, 0.9, 'e)', 'Units', 'Normalized');
ylim([-30 0]);
title('Cross-isobath extent of along-shelf supply jet');
linex([190 230]);
handles3.hplt(1).Color = 'k';
beautify([12 13 14], 'Times');

hax(6) = subplot(4, 2, [7 8]);
hbc = shfric2.plot_ts('run.shelfbc.shelf(:, 2)', hax(6), ...
                      'run.shelfbc.time/86400');
hbc.hplt(1).Color = 'k';
axes(hax(6)); legend('off');
title('');
ylabel('BC(t)');
xlabel('Time (days)');
linex([190 230]);
htxt(2) = text(0.04, 0.9, 'f)', 'Units', 'Normalized');
beautify([12 13 14], 'Times');

for ii=5:6
    hax(ii).XTickMode = 'auto';
    hax(ii).XTickLabelMode = 'auto';
end

hax(3).Position(2) = 0.55
hax(4).Position(2) = hax(3).Position(2);
hax(5).Position(3) = 0.6
hax(6).Position(3) = hax(5).Position(3);
hleg.Position(1) = 0.75;
hleg.Position(2) = 0.27;

axes(hax(5));
htxt(3) = text(1.1, 0.35, 'r_f', 'Units', 'normalized', ...
               'FontName', 'Times', 'FontSize', 12);

correct_ticks('y', [], {'50'; '100'}, hax([1 3]));
linkaxes(hax([5,6]), 'x');
hax(5).XLim = [120 320];

set(gcf, 'Renderer', 'painters');

export_fig -r300 images/paper3/shfric-ssh-scale-bc.png

%% multipanel bottom velocity plot.
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-8342-2')
    ew = runs('../topoeddy/runew-8342-2/');
end

figure;

clf('reset');
time = '311';
xloc = '325e3';
opt = [];
opt.drawtrack = 0;
opt.rhocontourplot = 0;
opt.addvelquiver = 1;
opt.limx = [180 400];
opt.limy = [0 150];
opt.nocolorbar = 0;

%surfvar = 'csdye'; clim = [-50 250];
surfvar = 'rho'; clim = [21.75 21.77];

% surface dye + surface velocity
opt.dxi = 7; opt.dyi = opt.dxi
hax(1) = subplot(311);
opt.quiverloc = 'surf';
h1 = ew.animate_field(surfvar, hax(1), time, 1, opt);
h1.hcb.Location = 'SouthOutside';
h1.hcb.Label.String = 'Surface \rho -1000 kg/m^3';
h1.htlabel.delete;
caxis(clim);
htxt(1) = text(0.05, 0.1, 'a) Surface velocity vectors', ...
               'Units', 'Normalized', 'Color', 'w', ...
               'FontSize', 12);
beautify([12, 13, 14], 'Times')

% below: surface dye bottom velocity
opt.quiverloc = 'bot';
opt.dxi = 5; opt.dyi = opt.dxi;
hax(2) = subplot(312);
h2 = ew.animate_field(surfvar, hax(2), time, 1, opt);
h2.htlabel.delete;
caxis(clim);
htxt(2) = text(0.05, 0.1, 'b) Bottom velocity vectors', ...
               'Units', 'Normalized', 'Color', 'w', ...
               'FontSize', 12);
beautify([12, 13, 14], 'Times')
colorbar('off')

% yz v
hax(3) = subplot(313); cla('reset')
handles = ew.PlotSingleYZSection('v', time, xloc, hax(3));
hold on
contour(handles.hvar.XData, handles.hvar.YData, handles.hvar.ZData, ...
        [1 1]*-1e-4, 'g', 'LineWidth', 2);
hcb = center_colorbar;
hcb.Location = 'southoutside';
caxis([-1 1]*0.015);
title('');
xlim([0 70]);
hax(3).YTickMode = 'auto';
hax(3).XTickMode = 'auto';
ylim([-450 0]);
handles.hcen.delete;
handles.htlabel.Position(2) = 0.3;
htxt(3) = text(0.05, 0.1, 'c) Cross-shelf velocity (m/s)', ...
               'Units', 'Normalized', 'Color', 'w', ...
               'FontSize', 12);
handles.hlz.delete;
beautify([12, 13, 14], 'Times')

linkaxes(hax(1:2), 'xy')
resizeImageForPub('onecolumn');

hax(1).Position = [0.170    0.7101    0.7661    0.2149];
hax(2).Position = [0.165    0.4096    0.7750    0.2157];
hax(3).Position = [0.165    0.1687    0.7750    0.1909];
h1.hcb.Position(1:3) = [0.35, 0.68, 0.4];

hax(1).XTickLabels = []
hax(1).XLabel.String = '';
hax(1).YTick = hax(2).YTick;
hax(1).YTickLabels = hax(2).YTickLabels;

hax(3).Position(2) = 0.13;

hcb.Position(1) = h1.hcb.Position(1);
hcb.Position(3) = h1.hcb.Position(3);

for ii=1:2
    axes(hax(ii))
    linex(str2double(xloc)/1000, [], 'k');
    hax(ii).Title.String = '';
end

for hh=[h1, h2]
    hh.htlabel.Color = 'w'
end

cmap = brighten(cbrewer('seq', 'Reds', 30), 0.4);
colormap(hax(1), cmap)
colormap(hax(2), cmap);

% hack. Don't why this is needed. (╯°□°）╯︵ ┻━┻
% hax(1).Position = [0.3014 0.64    0.4319    0.2203];
% hax(3).Position(1) = 0.18;

export_fig -r600 images/paper3/eddy-inflow.png

%% cross-shelfbreak hovmoeller
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-8342-2')
    ew = runs('../topoeddy/runew-8342-2/');
end

ew.hovmoeller('vbar', 'y', 'sb', 1);
title('Depth averaged cross-isobath velocity at shelfbreak');
% [start,stop] = ew.flux_tindices(ew.csflux.off.slope(:,1,1));
% liney(ew.eddy.t([start stop]));
hax = gca;
hax.YTickLabel{find_approx(hax.YTick, 240)} = 't_{start}';
hax.YTickLabel{find_approx(hax.YTick, 491)} = 't_{stop}';
correct_ticks('y', [], {'500'});
xlim([50 500]);

resizeImageForPub('portrait');
beautify([12, 13, 14], 'Times')
pbaspect([1.618 1 1]);

% hax.Position(3) = 0.5;
htxt(1) = text(1.17, 0.95, 'Offshore', 'Units', 'normalized', ...
               'Rotation', -90, 'FontSize', 12);
htxt(2) = text(1.17, 0.3, 'Onshore', 'Units', 'normalized', ...
               'Rotation', -90, 'FontSize', 12);
hcb = findobj('type', 'colorbar');
hcb.Label.String = '(m/s)';
hcb.Label.Rotation = 0;
hcb.Label.FontSize = 12;
hcb.Label.Position(1) = 2.5;
hcb.Label.VerticalAlignment = 'middle';

export_fig -r200 images/paper3/ew-8342-2-vbar-hov.png

%% buoyancy arrest
c6 = 24;
c5 = 39;
c2 = 1.6;
c1 = 0.6;

s = 0.05;
f = 5e-5;
N = sqrt(1e-5);
r = 3e-3;
Q = 5000;
alpha = 8e-4;
h0 = 50;
U = Q/h0/12e3

d = r/U * N/f
WT = c2 * (Q/N/alpha^2)^(1/3)
WB = c1/s * sqrt(Q/f/h0)

lt = c6/c2 * (1+s^2)/(d*s) * WT/1000
lb = c5/c1 * (1+s^2)/(d*s) * WB/1000

Lbeta = 12e3
f*Lbeta/N/r/86400

%% arrested topographic wave
% along-shelf scale (km) required for flow to spread across entire shelf
shfric2.print_params('phys.f0/misc.rdrg * bathy.sl_shelf  * (bathy.L_shelf)^2/1e3')

%shfric2.print_params('bathy.hsb/misc.rdrg');

kappa = shfric2.print_params('misc.rdrg/(phys.f0 * bathy.sl_shelf)')
sqrt(kappa*33e3)/1e3

%%
chi = sh.print_params(['(2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
                    'V0/Lz0/(bathy.S_sl*sqrt(phys.N2))']);
phii = sh.print_params('(V0./bathy.S_sl/sqrt(phys.N2))/bathy.hsb');
sh.print_params('bathy.sl_slope * sqrt(phys.N2)/phys.f0')
sh.print_params('(1-erf(bathy.hsb/Lz0))');
phii./ (1-chi)

%% ew-82343 - thesis only
% ew = runs('../topoeddy/runew-82343/');

% opt = [];
% opt.addvelquiver = 1;

% handles = ew.animate_field('csdye', [], '414', 1, opt);
% caxis([-50 350]);
% ylim([0 150]);
% xlim([150 580]);
% title('Cross-shelf dye | \lambda = 0.5');
% colorbar('off');
% handles.htlabel.Position(2) = 0.1

% export_fig -r150 -a2 -opengl images/paper3/ew-82343-on-shelf.png

%% ew-8341 - multipanel cross-sections

if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-8341')
    ew = runs('../topoeddy/runew-8341');
end

tindex = '300'
tindex = find_approx(ew.time, str2double(tindex)*86400, 1);
rhocontours = linspace(21.77, 22.3, 60);
xloc = '350e3'; % ''270e3';
eddcolor = ew.eddyeColormap();
eddcolor = eddcolor(end, :);

opt = []; handles = {};
figure;
hax(1) = subplot(3,2,[1, 2]); cla;
handles{1} = ew.animate_field('csdye', hax(1),  tindex, 1, opt);
handles{1}.htlabel.delete;
xlabel(''); title('');
handles{1}.hcb.Label.String = 'Cross-shelf dye (km)';
hl1 = linex(str2double(xloc)/1e3);
hax(1).YLim = [0 150];
hax(1).XTickMode = 'auto';
hax(1).XTickLabels = {};
hax(1).Position(2) = 0.67;

%=============================================================
% along-shelfbreak section
opt = [];
opt.rhocontours = 1;
opt.eddy0 = 0;
hax(3) = subplot(3,2,[3, 4]); cla
hax(3).Position(2) = 0.44;
handles{3} = ew.PlotSingleXZSection('v', 1, tindex, opt, hax(3));
handles{3}.hmask.Color = ew.shelfSlopeColor;
handles{3}.htext.delete;
handles{3}.htime.delete;
handles{3}.hline.delete;
handles{3}.hcb.Position(1) = handles{1}.hcb.Position(1);
handles{3}.hcb.Label.String = {'Cross-shelf velocity (m/s)';  'at shelfbreak'};
title('')
linkaxes(hax([1 3]), 'x');

% axes(hax(3))
% pos = handles{3}.hcb.Position;
% hvel(1) = text(pos(1)+0.37, pos(2)+0.57, 'offshore', ...
%                'FontName', 'Times', 'FontSize', 12, ...
%                'Units', 'normalized');
% hvel(2) = text(pos(1)+0.37, pos(2)-0.43, 'onshore', ...
%                'FontName', 'Times', 'FontSize', 12, ...
%                'Units', 'normalized');


%=============================================================
% two cross-shelf cross-sections
resizeImageForPub('portrait');
opt = [];
opt.csdfront = 1;
opt.rhocontours = 1;

hax(2) = subplot(3,2,5); cla
% handles = ew.overlay_section('v', 'rho', tindex, {'y' 1 40}, ...
%                              'x', xloc, 'center_colorbar');
handles{2} = ew.PlotSingleYZSection('u', tindex, xloc, hax(2), opt);
xlim([0 ew.bathy.xsb + 10e3]/1e3)
ylim([-80 0])
[eddye, ~,yy, zz] = dc_roms_read_data(ew, ew.eddname, tindex, ...
                                      {'x' xloc xloc});
[~, handles{2}.hedd] = ...
    contour(yy/1000, zz, eddye, [1 1]*0.7, 'Color', eddcolor, 'LineWidth', 2);
handles{2}.hrhocont.LevelList = rhocontours;
title('');
hax(2).XTickLabelMode = 'auto';
hax(2).XTick = [0:20:60];
% pbaspect([1.618 1 1])
handles{2}.hcsdcont.Color = ew.shelfSlopeColor;
handles{2}.hleg = legend([handles{2}.hcsdcont handles{2}.hedd], ...
                         'Shelf-water front', 'Eddy water front', ...
                         'Location', 'SouthWest');
handles{2}.hleg.Box = 'off';
handles{2}.hleg.TextColor = [1 1 1];
handles{2}.htlabel.delete;
colorbar('off');
hax(2).YTickMode = 'auto';

%=============================================================
hax(4) = subplot(3,2,6); cla
handles{4} = ew.PlotSingleYZSection('u', tindex, xloc, hax(4), opt);
[eddye, ~,yy, zz] = dc_roms_read_data(ew, ew.eddname, tindex, ...
                                      {'x' xloc xloc});
[~, handles{4}.hedd] = contour(yy/1000, zz, eddye, [1 1]*0.7, ...
                               'Color', eddcolor, 'LineWidth', 2);
handles{4}.hrhocont.LevelList = rhocontours(1:2:end);
title('');
xlim([20 60]);
ylim([-300 0]);
hax(4).XTickLabelMode = 'auto';
hax(4).XTick = [0:20:60];
handles{4}.hcsdcont.Color = ew.shelfSlopeColor;
handles{4}.hleg = legend([handles{4}.hcsdcont handles{4}.hedd], ...
                         'Shelf-water front', 'Eddy water front', ...
                         'Location', 'SouthWest');
handles{4}.hleg.Box = 'off';
handles{4}.hleg.TextColor = [1 1 1];
handles{4}.htlabel.delete;
hax(4).YLabel.String = '';
hax(4).YTickMode = 'auto';
hcb = colorbar;
hcb.Location='SouthOutside';
hcb.Label.String = 'Along-shelf velocity (m/s) at x = 300 km';

hcb.Position(2) = 0.07;
hcb.Position(1) = 0.34;

lab = 'acbd';
color = 'wwkw';
if exist('htxt', 'var'), htxt.delete; end
for ii=1:4
    axes(hax(ii));
    set(hax(ii), 'color', 'none');
    htxt(ii) = text(0.05, 0.32, [lab(ii) ')'], ...
                    'Color', color(ii), 'FontSize', 13, ...
                    'Units', 'normalized')
    beautify([12 13 14], 'Times');
end

hax(2).Position(2) = 0.14;
hax(4).Position(2) = hax(2).Position(2);
hax(2).Position(3) = 0.35;
hax(4).Position(3) = hax(4).Position(3);

handles{2}.hleg.Position(1) = 0.15;
handles{2}.hleg.Position(2) = 0.145;
handles{4}.hleg.Position(2) = handles{2}.hleg.Position(2);
handles{4}.hleg.Position(1) = 0.58;

hax(1).XLim = [100 400];
hax(1).Position(1) = 0.097;
hax(3).Position(1) = hax(1).Position(1);
hax(3).PlotBoxAspectRatio = hax(1).PlotBoxAspectRatio

colormap(hax(2), cbrewer('div', 'BrBG', 11));
colormap(hax(4), cbrewer('div', 'BrBG', 11));

set(gcf, 'renderer', 'painters')
export_fig -a2 -r300 images/paper3/multipanel-cross-sections.png
