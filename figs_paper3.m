%% shelf slopes
folders = { ...
    'ew-8041', 'ew-8042', ...
    ... %'ew-8150', 'ew-8151', ...
    'ew-82342', 'ew-82343', ... %'ew-82344' ...
    'ew-8341', 'ew-8361', ...
    'ew-8352', 'ew-8352-2', 'ew-8342-2', ...
    'ew-8383', 'ew-8384', 'ew-8385', ...
    'ew-8392'... %, 'ew-8346', ...
    'ew-583411', 'ew-583413', ...
    'ew-583414', 'ew-583415', ...
          };
sh = runArray(folders);

for ii=1:sh.len
    sh.name{ii} = sh.name{ii}(4:end);
end

%% inflow/outflow vertical profiles
handles = sh.PlotFluxVertProfiles;

export_fig -r150 images/paper3/sb-vert-profiles.png

%% fluxes
sh.filter = [1 2 5:13];
sh.print_diag('avg flux');
title('');
export_fig images/paper3/sh-avgflux.png

%% eddy water on shelf scale
figure; maximize;
hax(1) = subplot(121);
hax(2) = subplot(122);

sh.print_diag('supply', [], hax(1), 'no_name_points');
sh.print_diag('eddyonshelf', [], hax(2), 'no_name_points');

hax(1).Title.String = '(a)';
hax(2).Title.String = '(b)';

hax(2).YLim(1) = 0;

export_fig -r150 images/paper3/parameterizations.png

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

handles.supax.Position(4) = 0.7;
handles.hsuptitle.String = ...
    'Surface cross-shelf dye (km) | Initial Ro = 0.1 | Initial eddy scales = (25 km, 400m)';

axes(handles.hax(1));
[hleg,icons] = legend([handles.hfield{1}.hcen, ...
                    handles.hfield{1}.htrack, ...
                    ... %handles.hfield{1}.hrho, ...
                    handles.hfield{1}.hzeta], ...
                      {'Eddy center', 'Track of eddy center', ...%'Eddy core',
                    'SSH'}, 'Location', 'NorthWest'); %, 'FontSize', 14);
hleg.Box = 'off';
hleg.Position(2) = hleg.Position(2) - 0.03;
icons(end).Children.Children(1).LineWidth = 1;
icons(end).Children.Children(2).LineWidth = 1;
icons(end).Children.Children(3).LineWidth = 1;

% annofs = 14;
% hanno(1) = annotation('textarrow', [0.54 0.57], [1 1]*0.733, ...
%                       'String', 'wake cyclone','LineWidth', 1, 'FontSize', annofs, ...
%                       'HeadStyle', 'none');
% hanno(2) = annotation('textarrow', [0.33 0.31], [0.262 0.28], ...
%                       'String', 'leakage','LineWidth', 1, 'FontSize', annofs, 'Color', 'w', ...
%                       'HeadStyle', 'none');
% hanno(3) = annotation('textarrow', [1 1]*0.28, [0.36 0.33], ...
%                       'String', {'secondary'; 'cyclone'},'LineWidth', 1, 'FontSize', annofs, ...
%                       'HeadStyle', 'none', 'HorizontalAlignment', 'center');
% hanno(4) = annotation('textarrow', [0.56 0.6], [0.262 0.28], ...
%                       'String', 'anticyclonic eddies', 'LineWidth', 1, ...
%                       'FontSize', annofs, 'Color', 'w', 'HeadStyle', 'none');

export_fig -r150 -a2 -opengl images/paper3/ew-8341-surface-csdye.png

%% shelf baroclinicity
sh.print_diag('shelfbc');
title('');
axis square;
xlim([0 1]);

export_fig -r150 -a2 images/paper3/shelfbc.png

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

%% volume budget
handles = sh.array(4).VolumeBudget;
xlim([0 450]);
handles.hax.XLabel.Position(1) = 400;
handles.hisponge.DisplayName = 'Along-shelf: supply';

export_fig -r150 -a2 images/paper3/ew-8341-volume-budget.png

%% flat and sloping, sbssh
if ~exist('ew34', 'var')
    ew34 = runArray({'ew-34', 'ew-8341'});
end

ew34.name = {'flat shelf'; 'sloping shelf'};
figure; maximize;
hax(1) = subplot(121); hold on;
ew34.plot_avgProfile('zeta', 1, 'y', 'sb', 1, hax(1));
title('Mean SSH at shelfbreak (normalized)');
pbaspect([1.617 1 1]);
hax(2) = subplot(122); hold on;
ew34.plot_avgProfile('zeta', 0, 'y', 'sb', 1, hax(2));
title('Mean SSH at shelfbreak');
pbaspect([1.617 1 1]);

export_fig -r120 images/sb-ssh-ew34-flat-sloping.png

%% compare flat and sloping
if ~exist('ew04', 'var')
    ew04 = runArray({'ew-04', 'ew-8041'});
end

opt.drawtrack = 0;
opt.drawcen = 0;
opt.addvelquiver = 1;
opt.rhocontourplot = 0;
opt.csdcontourplot = 0;
opt.dxi = 8; opt.dyi = 5;

ew04.array(1).read_velsurf(find_approx(ew04.array(1).time, 169*86400, 1));
ew04.array(1).read_csdsurf(find_approx(ew04.array(1).time, 169*86400, 1));

clear handles
figure; maximize;
hax = packfig(2,2);
for ii=1:2
    run = ew04.array(ii);
    tind = find_approx(run.time, 169*86400, 1);

    handles(ii) = run.animate_field('csdye', hax(ii), '169', 1, opt);
    xlim([170 400]);
    ylim([0 120]);

    handles(ii).hbathy{1}.ShowText = 'off';
    handles(ii).hquiv.Color = [1 1 1]*0;
    handles(ii).hquiv.LineWidth = 1.5;
    handles(ii).hquiv.AutoScaleFactor = 4.5;
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
    hyy(ii,1).XLim = [170 400];
    hyy(ii,2).XLim = [170 400];

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
handles(1).hax.Title.FontSize = 26;
handles(2).hax.Title.FontSize = 26;

handles(1).supax.Position(4) = 0.73;
handles(1).supax.Title.String = 'Surface cross-shelf dye (km)';
handles(1).supax.Title.FontSize = 28;

hax(3).YLabel.String = {'Cross-isobath'; 'surface velocity (m/s)'};
hax(3).YLabel.Color = hplt(1,1).Color;
hax(3).XLabel.String = 'X (km)';
hyy(2,2).YLabel.String = 'Cross-shelf dye - Y_{sb}';
hyy(2,2).YLabel.Color = hplt(1,2).Color;
hyy(2,2).XLabel.String = 'X (km)';

axes(hyy(1,1))
htxt(1) = text(200, 0.04, 'offshore', 'Color', hplt(ii,1).Color);
htxt(2) = text(200, -0.04, 'onshore', 'Color', hplt(ii,1).Color);

axes(hyy(2,2));
htxt(3) = text(320, 150, 'eddy water', 'Color', hplt(ii,2).Color);
htxt(4) = text(320, -20, 'shelf water', 'Color', hplt(ii,2).Color);

axes(hax(2))
hl = liney(37.5-12*1.37);
hl.Color = [1 1 1]*0;
hanno = annotation('doublearrow', [0.6 0.6], 0.615 + [0 0.04]);
hanno.Head1Style = 'cback3';
hanno.Head2Style = 'cback3';
hanno.Head1Length = 7;
hanno.Head2Length = 7;
hanno.Head1Width = 7;
hanno.Head2Width = 7;
htxt = text(185, 30, '1.37L_\beta', 'Color', 'k', 'FontSize', 20);

export_fig -opengl -r150 -a2 images/paper3/sbsnapshot.png

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
hanno = annotation('doublearrow', [1 1]*.55, [0 0.075]+0.545);
hanno.LineWidth = 1; hanno.Color = hlbeta.Color;
hanno.HeadSize = 8; hanno.HeadStyle = 'vback3';
htext = text(430, -Lbeta/2000, 'L_\beta', 'Color', hlbeta.Color);
linkprop([htext hanno hlbeta], 'Color');
htext.Color = handles2.htxt(3).Color;

export_fig -r150 -a2 -png images/paper3/sb-flux-summary

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

% for ii=1:4
%     shfric2.array(ii).usurf = [];
%     shfric2.array(ii).vsurf = [];
%     shfric2.array(ii).read_velsurf;
% end

opt.addvelquiver = 1;
opt.quivercolor = [0 150 136]/255;
opt.scale = 5;
opt.limy = [0 150];
shfric2.filter = [1 2 3 4];
handles = shfric2.mosaic_field('csdye', [350 370 280 280], opt, [-40 150]);

%% along-shelf inflow structure
sh.filter = [1 2 5:13];
%sh.plot_avgProfile('vbar');
figs = [0 0 0 0 1];
sh.plot_fluxes(1,1,1,figs);

%% friction plot
shfric2.sort(shfric2.print_params('run.params.misc.rdrg'));
shfric2.array(1).read_zeta;
shfric2.array(1).read_velsurf;

shfric2.array(end).read_zeta;
shfric2.array(end).read_velsurf;

opt.limy = [0 150];
opt.limx = [200 480];
opt.quivercolor = [0 150 136]/255; shfric2.array(1).shelfSlopeColor('dark'); %[1 0 0];
opt.addvelquiver = 1;
opt.addzeta = 0;
opt.rhocontourplot = 0;
opt.csdcontourplot = 0;
opt.addzeta = 1;
opt.commands = 'caxis([-1 1]*1e-3)';
opt.dyi = 8;
opt.dxi = 8;
opt.scale = 0;
opt.uref = 5e-4; opt.vref = opt.uref;
str = 'acbd';
trackcolor = 'k';
zetacolor = [244 67 54]/255; [1 1 1]*0.6;
var = 'zeta';
%figure; maximize;
clf; clear handles;
hax = packfig(2,3);
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
        title(['SSH (m) | r_f = 0']);
    else
        title(['SSH (m) | r_f = ' num2str(run.params.misc.rdrg, '%1.1e')]);
    end
    % handles(ii,1).htlabel.Position(2) = 0.14
    handles(ii,1).htlabel.String = [str(ii) ') ' handles(ii,1).htlabel.String];
    handles(ii,1).htrack.Color = trackcolor;
    xlabel('');
    handles(ii,1).hfield.Visible = 'off';

    U = handles(ii,1).hquiv.UData;
    V = handles(ii,1).hquiv.VData;
    y = handles(ii,1).hquiv.YData;
    speed = hypot(U,V);
    %nanmask = speed > 0.175*max(speed(:)) | y > 80;
    nanmask = y > 38;
    handles(ii,1).hquiv.UData(nanmask) = NaN;
    handles(ii,1).hquiv.VData(nanmask) = NaN;

    axes(hax(ii+3)); % = subplot(2,3,ii+3);
    handles(ii,2) = run.animate_field(var, hax(ii+3),  '230', 1, opt);
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
    handles(ii,2).hfield.Visible = 'off';

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
end
hax(2).YLabel.String = '';
hax(2).YTickLabel = {};
hax(5).YLabel.String = '';
hax(5).YTickLabel = {};

% hcb = colorbar('northoutside');
% xlabel('');
% hcb.Position(1) = 0.3;
% hcb.Position(2) = 0.46;
% hcb.Ruler.Exponent = 0;

hax(3) = subplot(233);
plot(run.csflux.time/86400, run.csflux.off.slope(:,1,1)/1000);
xlim([0 400]);
hold on;
plot(shfric2.array(1).csflux.time/86400, shfric2.array(1).csflux.off.slope(:,1,1)/1000, 'k');
title({'Offshore flux of shelf water (mSv)'})
ylim([0 10]);
hleg = legend('r_f = 3e-3', 'r_f = 0');
hleg.Location = 'NorthWest';
linex([190 242]);
htxt(1) = text(30, 2, 'e)');
hax(3).Position(1) = 0.68;
hax(3).Position(3) = 0.3;
pbaspect([1.618 1 1]);
beautify;

hax(6) = subplot(236);
plot(run.shelfbc.time/86400, run.shelfbc.shelf(:,2));
hold on;
plot(shfric2.array(1).shelfbc.time/86400, shfric2.array(1).shelfbc.shelf(:,2), 'k');
xlim([0 400]);
title('BC_{0.2}');
xlabel('Time (days)');
linex([190 242]);
htxt(2) = text(30, 0.2, 'f)');
hax(6).Position(1) = hax(3).Position(1);
hax(6).Position(3) = hax(3).Position(3);
pbaspect([1.618 1 1]);
beautify;

correct_ticks('y', [], {'50'; '100'}, hax([1 2]));
correct_ticks('x', [], {'200'}, hax([3 6]));

export_fig images/paper3/shfric-ssh-flux-bc.png