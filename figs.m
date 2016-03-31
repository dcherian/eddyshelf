%%
if ~exist('ewshdp', 'var')
    ewshdp = runArray({'ew-64461-9-shallow', 'ew-64461-9-deep', ...
                       'ew-64361-shallow', 'ew-64361-deep'});
end

fontsize = [22 24 28];
figure; maximize;
ax(1) = subplot(121);
ewshdp.filter = [1 2];
ewshdp.plot_ts('eddy.hcen', ax(1))
pbaspect([1.618 1 1]);
title('$$\frac{U}{\beta L^2} \sim 20 $$');
ax(1).Title.Interpreter = 'latex';
legend('off');
beautify(fontsize);

ax(2) = subplot(122);
ewshdp.filter = [3 4];
ewshdp.plot_ts('eddy.hcen', ax(2));
pbaspect([1.618 1 1]);
title('$$\frac{U}{\beta L^2} \sim 60 $$');
ax(2).Title.Interpreter = 'latex';
legend('off');
beautify(fontsize);

linkaxes(ax, 'y');
packfig(1,2, 'columns');
ax(2).XTick(1) = [];

[ax(3), htitle] = suplabel('Water depth at eddy center', 't');
ax(3).Position(4) = 0.82;
ax(3).FontSize = fontsize(end);
beautify(fontsize);

export_fig -r150 images/shallow-deep-hcen.png

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

handles = shfric.PlotFluxSummary(1);

axes(handles.hax(1))
handles.hax(1).Position(2) = 0.68;
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
htable(1) = text(1.21, pos1(2)+0.08, 'S_{sh}   r', 'Units', 'normalized');
htable(2) = text(1.45, pos1(2)+0.08, 'S_{sh}     r', 'Units', 'normalized');
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
handles2 = shfric.plot_avgProfile('zeta', 'y', 'sb', 1, handles.hax(3));
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
                        'HorizontalAlignment', 'center');
handles2.htxt(3) = text(-150, 2e-4, 'Sloping shelf', 'Units', 'data', ...
                        'HorizontalAlignment', 'center');
handles2.htxt(4) = text(-50, 1e-3, 'Shelf water', 'Units', 'data', ...
                        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
handles2.htxt(5) = text(150, 1e-3, 'Eddy water', 'Units', 'data', ...
                        'HorizontalAlignment', 'center');
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

%export_fig -r150 -a2 -painters -pdf -png images/paper2/sb-flux-summary


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
