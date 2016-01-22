%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% field map

if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end
factor = 2;
isobath = 4;
timesteps = [1 150 200 250 300 380];

opt.rhocontourplot = 1;
opt.csdcontourplot = 1;
opt.csdcontours = ew34.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.zetaOnFirstPlot = 1;
handles = ew34.mosaic_field('csdye', timesteps, opt);
handles.hfield{1}.hzeta.LevelList = handles.hfield{1}.hzeta.LevelList(1:2:end);
handles.hcb.delete;

ylim([0 250]);
correct_ticks('y', [], {'50', '100'}, handles.hax([1 4]));

handles.supax.Position(4) = 0.715;
handles.htitle.String = 'Surface cross-shelf dye (km) | Ro = 0.1 | L = 25 km | L^z = 400m';

axes(handles.hax(1));
[hleg,icons] = legend([handles.hfield{1}.hcen, ...
                    handles.hfield{1}.htrack, ...
                    handles.hfield{1}.hrho, ...
                    handles.hfield{1}.hzeta], ...
                      {'Eddy center', 'Track of eddy center', 'Eddy core', 'SSH'}, ...
                      'Location', 'NorthWest'); %, 'FontSize', 14);
hleg.Box = 'off';
hleg.Position(2) = hleg.Position(2) - 0.03;
icons(end).Children.Children(1).LineWidth = 1;
icons(end).Children.Children(2).LineWidth = 1;
icons(end).Children.Children(3).LineWidth = 1;

export_fig -r120 -a2 -painters images/paper2/ew-34-surface-csdye.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
set(gcf, 'Position', [675 175 1080 924]);
columnlegend(3, names, 'Location', 'NorthWest');
export_fig -r120 -painters -a2 -png -pdf images/paper2/centracks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flux diagnostics

handles = ew34.plot_fluxts(factor, 3, 3);
handles.htitle.String = ['Flux of water for z > ' num2str(factor) ' H_{sb}'];
color = handles.icons(1).Color;
maximize; %drawnow;
% mark timesteps
axes(handles.hax(1));
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
for ii=[1 2 3 6 8]
    handles.icons(ii).Color = color;
end
hpt = findobj(handles.icons, 'Type', 'patch');
hpt.FaceColor = color;
hpt.FaceAlpha = 0.2;

axes(handles.hax(2)); % make visible
export_fig -r150 -a2 -pdf -png images/paper2/flux-diags

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface map of eddye for sampling schematic

if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface map and instantaneous flux

folders = {'ew-4341', 'ew-34', 'ew-2341_wider'};
if ~exist('ewflux', 'var'), ewflux = runArray(folders); end
N = length(folders);
times = [180 300 275];
opt.rhocontourplot = 0;
opt.csfluxplot = 1;
opt.nocolorbar = 1;
opt.csdcontourplot = 0;
opt.addvelquiver = 0;
opt.csfluxIsobath = 1;
clear handles
figure;
for ii=1:3
    run = ewflux.array(ii);
    [~,~,tres] = run.locate_resistance;
    hax(ii) = subplot(5,3,[1 4 7] + (ii-1));
    hax(N+ii) = subplot(5,3,[10 13] + (ii-1));
    opt.dy = 0; run.bathy.xsb/1000;

    handles(ii) = run.animate_field('csdye', hax([ii N+ii]), times(ii), 1, opt);
    handles(ii).hfield.CData = handles(ii).hfield.CData - run.bathy.xsb/1000;

    if opt.csdcontourplot
        handles(ii).hcsd.ZData = handles(ii).hcsd.ZData - run.bathy.xsb;
        handles(ii).hcsd.LevelList = handles(ii).hcsd.LevelList - run.bathy.xsb;
    end

    axes(hax(ii)); caxis([-30 200]);
    handles(ii).htitle.String = ['\lambda = H_{sb}/L_z = ', ...
                        num2str(run.bathy.hsb/run.eddy.Lgauss(tres), '%.2f')];
    handles(ii).hbathy{1}.ShowText = 'off';
    axes(hax(ii));
    hcen2(ii) = linex(handles(ii).hcen.XData);
    hcen2(ii).YData(2) = hcen2(ii).YData(2) * 3/4;
    xlabel(''); hax(ii).XTickLabel = []; colorbar('off');
    axes(hax(N+ii));
    title(''); pbaspect([1 0.5 0.5]); hax(N+ii).Position(2) = 0.18;

    if ii ~= 1
        axes(hax(ii)); ylabel('');
        axes(hax(N+ii)); ylabel('');
    end
end

hax(5).YLim = [-1 1]*0.4;
hax(6).YLim = [-1 1]*5;
correct_ticks('y', [], '100', hax(1));
correct_ticks('y', [], {'50'; '100'}, hax(2));
correct_ticks('y', [], {'150'; '200'}, hax(3));
correct_ticks('x', [], '400', hax(4));
correct_ticks('x', [], '200', hax(5));
correct_ticks('x', [], {'300'; '400'}, hax(6));

axes(hax(4));
ylabel('\int v(x,z) dz (m^2/s)');
htxt(1) = text(0.70, 0.15, ' Onshore', 'Color', [1 1 1]*0.55, ...
               'Units', 'Normalized');
htxt(2) = text(0.70, 0.85, 'Offshore', 'Units', 'Normalized');
linkprop(htxt, 'Color');

axes(hax(N));
hcb = colorbar('southoutside');
pos = hcb.Position;
hcb.Position(1) = 0.5 - pos(3)/2;
hcb.Position(2) = 0.5 + pos(4)/2;
hcb.Label.String = 'Cross shelf dye - Y_{sb} (km)';

axes(hax(2))
htxt(3) = text(0.7, 0.7, '(Figure 1e)', 'Units', 'normalized', ...
               'Color', 'k', 'FontSize', 16);

supax = suplabel('Along-shelf profile of cross-isobath velocity at shelfbreak', 't');
supax.Position(end) = 0.79;

export_fig -painters -r120 -a2 images/paper2/inst-flux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x-z sections

handles = ew34.plot_xzsection(isobath, 225);
correct_ticks('y', [], '-300', handles.hax);
delete(handles.hrunname);
for ii=3:4
    for jj=1:2
        delete(handles.hline(ii).hl{jj});
        delete(handles.hline(ii).htxt{jj});
    end
end

for ii=1:2
    handles.hline(2).htxt{ii}.Units = 'normalized';
    handles.hline(2).htxt{ii}.Position(1) = 0.3;
end

handles.hax(1).XLim = [-170 150];
handles.hcb(1).Ticks = sort(unique([handles.hcb(1).Ticks -0.09 0.09]));

export_fig -r120 -opengl -a2 images/paper2/ew-34-xzsection.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% secondary eddy - 2360

if ~exist('ew2360', 'var') | ~strcmpi(ew2360.name, 'ew-2360_wider')
    ew2360 = runs('../topoeddy/runew-2360_wider/');
    xx = [343 346 349]';
    yy = [231 237 231]';
end
opt.addvelquiver = 0;
opt.csdcontourplot = 0;
opt.rhocontourplot = 0;
handles = ew2360.secondary_vortices(95, [xx yy], opt);
handles.hfield.htrack.delete;
correct_ticks('y', [], '198', handles.hax(1));
handles.hax(1).Title.String = 'Cross shelf dye - X_{sb} (km)';

export_fig('-r120', '-opengl', '-a2', ...
           'images/paper2/ew-2360-secondary-cyclone.png');

handles.hax(3).YLabel.String = 'Z (m)';
export_fig(handles.hax(3), '-painters', '-a2', ...
           'images/paper2/ew-2360-secondary-cyclone-rho.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% secondary eddy - ew-34

if ~exist('ew34', 'var')
    ew34 = runs('../topoeddy/runew-34/');
end
xx = [394 371 297 282 219 204]';
yy = [70.5 70.5 46.5 66 91 108]';
opt.addvelquiver = 0;
opt.csdcontourplot = 0;
opt.rhocontourplot = 0;
handles = ew34.secondary_vortices(320, [xx yy], opt);
hbathy = handles.hfield.hbathy;
clabel(hbathy{4}, hbathy{1}, 'LabelSpacing', 108*12);

% displace eddy water plots
dx = 230; dpv = 4e-11;
for ii=[1 3 5 7]
    handles.hcsd(ii).XData = handles.hcsd(ii).XData + dx;
end

axes(handles.hax(1));
correct_ticks('y', [], '198');
handles.hax(1).Title.String = 'Cross shelf dye - X_{sb} (km)';

axes(handles.hax(2)); linex(dx);
xticks = sort(unique([0 max(handles.hcsd(1).XData) ...
                    dx max(handles.hcsd(4).XData)]));
handles.hax(2).XTick = xticks;
xlim([-10 2*dx]);
handles.hax(2).XTickLabels = {'Sh', 'Edd', 'Sh', 'Edd'};

axes(handles.hax(3));
xlim([-0.07 0.05]);

axes(handles.hax(4));
for ii=[1 3 5 7]
    handles.hpv(ii).XData = handles.hpv(ii).XData + dpv;
end
linex(dpv);

axes(handles.hax(5));
xlim([-1 1]*0.25);

% extend horizontal lines
for ii=1:4
    axes(handles.hax(ii+1));
    handles.hl(ii,1).XData = xlim;
    handles.hl(ii,2).XData = xlim;
end

export_fig -r96 -a2 images/paper2/ew-34-secondary-cyclone.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ew-34 z-slice

opt.csdcontours = ew34.csflux.x(isobath);
handles = ew34.mosaic_zslice('dye_03', 75, [225 230 235 240 245 260], opt);
for ii=1:length(handles)
    handles(ii).eddsurf.delete;
    handles(1).htext(2).delete;
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
handles(1).htext(1).Position(2) = 0.18;
handles(1).htext(3).Position(2) = 0.11;
xlim([200 400]);
ylim([0 150]);
handles(1).htitle.String = 'Eddy dye at z = -75 m | H_{sb} = 50m | Ro = 0.10';
handles(1).htitle.FontSize = 22;
handles(1).htitle.FontWeight = 'normal';
handles(1).supax.Position(4) = 0.87;
axes(handles(4).hax); correct_ticks('x', '', '400');
axes(handles(5).hax); correct_ticks('x', '', '400');
export_fig -r96 -a2 images/paper2/ew-34-mosaic-zslice.png

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

export_fig -r120 -a2 images/paper2/ew-2360-mosaic-zslice.png

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

export_fig -r96 -a2 images/paper2/ew-2360-bottomrhopv.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% secondary cyclone x-z section

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

export_fig -r120 -a2 images/paper2/sb-flux-summary.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg streamer profiles - shelfbreak

handles = ew34.plotAvgStreamer(1);
correct_ticks('y', [], '-50');
handles.ax(1).Title.String = 'Mean cross-isobath velocity (m/s)';
export_fig -r120 -a2 images/paper2/ew34-avgstreamer-sb.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg streamer profiles - offshore

handles = ew34.plotAvgStreamer(4);
handles.ax(1).Title.String = 'Mean cross-isobath velocity (m/s)';
correct_ticks('y', [], {'-50'; '-400'}, handles.ax(1));
handles.ax(3).XLim = [0 0.6];
export_fig -r120 -a2 images/paper2/ew34-avgstreamer-sl.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flux vertical profiles

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
export_fig -r96 -a3 images/paper2/fluxvertprofile.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shelfbreak flow surface snapshot

if ~exist('ew04', 'var')
    ew04 = runArray({'runew-04', 'runew-8041'});
end

opt.addvelquiver = 1;
opt.rhocontourplot = 0;
opt.csdcontourplot = 0;
opt.dxi = 8; opt.dyi = 5;

handles = ew04.mosaic_field('csdye', {'169'; '169'}, opt);
xlim([170 400]);
ylim([0 120]);
for ii=1:2
    delete(handles(ii).hrunname)
    handles(ii).hbathy{1}.ShowText = 'off';
    handles(ii).hquiv.Color = [1 1 1]*0;
    handles(ii).hquiv.LineWidth = 1.5;
    handles(ii).hquiv.AutoScaleFactor = 4.5;
    handles(ii).htlabel.Position(2) = 0.1;
    handles(ii).htrack.delete;
    handles(ii).hbathy{2}.Color = [1 1 1]*0.9;
end
correct_ticks('x', [], '450', handles(1).hax);
handles(1).hax.Title.String = 'Flat shelf | S_{sh} = 0';
handles(2).hax.Title.String = 'Sloping shelf | S_{sh} = 0.05';

handles(1).supax.Position(4) = 0.73;
handles(1).supax.Title.String = 'Surface cross-shelf dye (km)';
handles(1).supax.Title.FontSize = 20;

axes(handles.hax(2))
hl = liney(37.5-12);
hl.Color = [1 1 1]*0.9;
hanno = annotation('doublearrow', [0.6 0.6], [0.405 0.445]);
hanno.Head1Style = 'cback3';
hanno.Head2Style = 'cback3';
hanno.Head1Length = 7;
hanno.Head2Length = 7;
hanno.Head1Width = 7;
hanno.Head2Width = 7;
htxt = text(202, 33, 'L_\beta', 'Color', [1 1 1]*0.9, 'FontSize', 20);

export_fig -opengl -r120 -a2 images/paper2/sbsnapshot.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg flux

hax = csf.plot_fluxparam('avg flux');
axes(hax(3));
hleg = legend;
hleg.Position(1) = 0.82;

export_fig -r96 -a2 images/paper2/avgflux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max flux

hax = csf.plot_fluxparam('max flux');
axes(hax(3));
hleg = legend;
hleg.Position(1) = 0.82;
ylim([-1 1]*55);

export_fig -r96 -a2 images/paper2/maxflux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3d schematics
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

mergename = 'images/paper2/3dmerged.png';
multifig = 1; % multiple figure windows and stitch later?
day = [220 230 240 260];
depth = 75;
MoveToZLevel = -990;

annocolor = 'k'; %[1 1 1]*0.5;
annofs = 23;
annolw = 0.5;
annoheadstyle = 'cback3';
annofontname = 'Futura Bk BT';

opt.eddthresh = 0.8;
opt.csdcontours = ew34.bathy.xsb+5000;
opt.eddreducepatch = 0.3;
opt.csdreducepatch = 0.3;
opt.finalize = 1;
opt.linefilter = 1;

clear handles hslice hplane hax
if ~multifig
    figure; maximize; set(gcf, 'Renderer', 'opengl');
end

for ii=1:length(day)
    opt.x = [200, 410] * 1000;
    opt.y = [300, 0] * 1000;
    if multifig
        figure; maximize; hax(ii) = gca;
        set(gcf, 'Renderer', 'opengl');
    else
        hax(ii) = subplot(2,2,ii);
    end

    if ii == 2
        opt.sect = 'x';
    else
        opt.sect = 'y';
    end

    handles(ii) = ew34.animate_3d(num2str(day(ii)), opt, hax(ii));
    handles(ii).hcsd.FaceAlpha = 0.8;

    hax(ii).Title.String = '';
    hax(ii).XAxis.Visible = 'off';
    hax(ii).YAxis.Visible = 'off';
    hax(ii).ZAxis.Visible = 'off';
    hax(ii).Box = 'off';
    hax(ii).DataAspectRatio = [1 1 5];

    bathy = handles(ii).hbathy.ZData;
    handles(ii).hbathy.ZData(bathy <= -850) = -850;
    handles(ii).hsectoutline.ZData(handles(ii).hsectoutline.ZData <= -850) = -850;
    handles(ii).hplane.ZData(handles(ii).hplane.ZData <= -850) = -850;
    handles(ii).hcsd.FaceColor = [107 174 214]/255;
    handles(ii).hedd.FaceColor = brighten([215 48 31]/255, 0.1);

    % hslice(ii) = ew34.animate_zslice('dye_03', depth, day(ii), hax(ii), opt);
    % hslice(ii).bathy{1}.delete;
    % hslice(ii).bathy{2}.delete;
    % hslice(ii).bathy{3}.delete;
    % hslice(ii).eddsurf.delete;
    % hslice(ii).htime.delete;
    % colorbar('delete');
    % linkprop([hslice(ii).hvar hslice(ii).csdsurf hslice(ii).rhocont], 'ContourZLevel');
    % hslice(ii).hvar.ContourZLevel = MoveToZLevel;
    % hslice(ii).hvar.ZData = fillnan(double(hslice(ii).hvar.ZData > 0.4), 0);

    % [xmat, ymat] = meshgrid(xlim, ylim);
    % hplane(ii) = surf(xmat, ymat, -1*abs(depth) * ones(size(xmat)), ...
    %                   'FaceColor', 'k', 'FaceAlpha', 0.05, 'EdgeColor', [1 1 1]*0.3);

    % houtline(ii) = plot3([xmat(:); xmat(1)], [ymat(:,1)' flip(ymat(:,2)') ymat(1)], ...
    %                      MoveToZLevel * ones([1 5]), ...
    %                      '-', 'LineWidth', 1, 'Color', [1 1 1]*0.3, 'Tag', 'dcline');

    % linkprop([hplane(ii) handles(ii).hyplane], 'FaceAlpha');

end

%linkprop(houtline, {'LineWidth', 'Color'});
%linkprop(hplane, 'EdgeColor');

linkprop(hax, 'DataAspectRatio');
linkprop([handles.hbathy], {'FaceColor', 'FaceAlpha'});
linkprop([handles.hedd], {'FaceColor', 'AmbientStrength'});
linkprop([handles.hcsd], {'FaceColor', 'AmbientStrength'});

zlim([-850 1.1]);
linkprop(hax, 'ZLim');

% change views and lighting as needed
%axes(hax(1));
%view(-125, 28);

axes(hax(1));
view(-110,28);
hanno(1) = annotation('textarrow', [0.64 0.56], [0.52 0.62], ...
                      'String', {'The cyclonic';  'wave is'; 'a wrinkle'; 'in 3-D.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'top');
hlabel(1) = annotation('textbox', [0.4 0.18 0.05 0.05], ...
                       'String', '(a)', 'EdgeColor', 'none');

axes(hax(2));
view(-122, 42);
handles(2).hlight.Position = [-1400 1400 500];
hanno(2) = annotation('textarrow', [0.645 0.53], [0.49 0.59], ...
                      'String', {'The wave';  'propagates'; ...
                    'around the eddy'; 'with the streamer.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hlabel(2) = annotation('textbox', [0.4 0.18 0.05 0.05], ...
                       'String', '(b)', 'EdgeColor', 'none');

axes(hax(3));
view(-130, 28);
hanno(3) = annotation('textarrow', [0.40 0.47], [0.78 0.7], ...
                      'String', {'The wave rolls up into'; ...
                    'a cyclone trapping'; ...
                    'the streamer water'; 'above it.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hanno(4) = annotation('textarrow', [0.66 0.57], [0.53 0.61], ...
                      'String', {'There is a';  'persistent bulge'; ...
                    'in the eddy'; ...
                    'below'; 'shelfbreak depth'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hlabel(3) = annotation('textbox', [0.4 0.18 0.05 0.05], ...
                       'String', '(c)', 'EdgeColor', 'none');

axes(hax(4));
view(-128, 40);
hanno(5) = annotation('textarrow', [0.42 0.46], [1 1]*0.75, ...
                      'String', {'Cyclone propagates'; ...
                    'away with trapped';  'streamer water.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hanno(6) = annotation('textarrow', [0.65 0.57], [0.53 0.61], ...
                      'String', {'The process'; 'repeats.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs);
hlabel(4) = annotation('textbox', [0.4 0.18 0.05 0.05], ...
                       'String', '(d)', 'EdgeColor', 'none');

for ii=1:length(hanno)
    try
        hanno(ii).HorizontalAlignment = 'center';
        hanno(ii).TextMargin = 0.05;
        hanno(ii).Text.FontName = annofontname;
        hanno(ii).FontSize = annofs;
        hlabel(ii).Text.FontName = annofontname;
        hlabel(ii).FontSize = annofs;
    catch ME
    end
end

tic;
if multifig
    for ii=1:4
        axes(hax(ii));
        export_fig('-r120', '-a4', '-p0.01', '-opengl', ...
                   ['images/paper2/3d-schem-' num2str(ii) '.png']);
    end
    hash = githash;
    system(['montage images/paper2/3d-schem-[1-4].png ' ...
            '-geometry +0.05+0.05 ' mergename]);
    system(['exiftool -overwrite_original -Producer=' hash ' ' mergename]);
else
    export_fig('-opengl', '-r120', '-a4', '-p0.02', mergename);
end
toc;

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