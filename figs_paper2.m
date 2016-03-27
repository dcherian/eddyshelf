%% params
isobath = 4;
factor = 1;
timesteps = [1 150 200 250 300 380];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% field map

if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

opt.rhocontourplot = 1;
opt.csdcontourplot = 1;
opt.csdcontours = ew34.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.addzeta = 1;
handles = ew34.mosaic_field('csdye', timesteps, opt);
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

handles.supax.Position(4) = 0.715;
handles.htitle.String = ...
    'Surface cross-shelf dye (km) | Ro = 0.1 | Eddy scales = (25 km, 400m)';

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

export_fig -r150 -a2 -painters images/paper2/ew-34-surface-csdye.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% center tracks

colors = [51,34,136; ...
          136,204,238; ...
          221,204,119; ...
          ... % 68,170,153; ...
          17,119,51; ...
          204,102,119; ...
          170,68,153]/255;

set(groot,'DefaultAxesColorOrder', colors);
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
xlim([-10 0.5]); ylim([0 3]);

hanno = annotation('textarrow', [0.65 0.72], [0.42 0.46], ...
                   'String', 'eddy splits here', ...
                   'LineWidth', 1, 'HeadStyle', 'none', ...
                   'HorizontalAlignment', 'center');

export_fig -r150 -a2 -png -pdf images/paper2/centracks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flux diagnostics

handles = ew34.plot_fluxts(factor, 3, 3);
handles.htitle.String = ['Flux of water for z > -H_{sb}'];
handles.htitle.FontSize = 28;
maximize; %drawnow;
% mark timesteps
axes(handles.hax(1));
ylim([0 25]*1e3);
handles.hax(1).YTick(handles.hax(1).YTick == 15e3) = [];
handles.hax(1).YTick = sort(unique([handles.hax(1).YTick 0 10e3 25e3]));
handles.hax(1).YTickLabelMode = 'auto';
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

export_fig -r150 -a2 -eps -png images/paper2/flux-diags

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
lab = 'abc';

figure;
for ii=1:3
    run = ewflux.array(ii);
    [~,~,tres] = run.locate_resistance;
    hax(ii) = subplot(5,3,[1 4 7] + (ii-1));
    hax(N+ii) = subplot(5,3,[10 13] + (ii-1));
    opt.dy = 0; run.bathy.xsb/1000;

    handles(ii) = run.animate_field('csdye', hax([ii N+ii]), times(ii), 1, opt);
    handles(ii).hfield.CData = handles(ii).hfield.CData - run.bathy.xsb/1000;
    handles(ii).htext2.delete;

    if opt.csdcontourplot
        handles(ii).hcsd.ZData = handles(ii).hcsd.ZData - run.bathy.xsb;
        handles(ii).hcsd.LevelList = handles(ii).hcsd.LevelList - run.bathy.xsb;
    end

    axes(hax(ii)); caxis([-30 200]);
    handles(ii).htitle.String = ['(' lab(ii) ') \lambda = H_{sb}/L_z = ', ...
                        num2str(run.bathy.hsb/run.eddy.Lgauss(tres), '%.2f')];
    handles(ii).hbathy{1}.ShowText = 'off';
    handles(ii).hbathy{2}.Color = 'w';
    handles(ii).hbathy{3}.Color = 'w';

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

hax(1).XLim = [200 600]; hax(1).YLim = [0 350];
hax(2).XLim = [120 460]; hax(2).YLim = [0 200];
hax(3).XLim = [150 480]; hax(3).YLim = [0 300];

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
hcb.Label.String = 'Cross shelf dye - Y_{sb} (km)';

axes(hax(2))
htxt(3) = text(0.7, 0.9, '(Figure 1e)', 'Units', 'normalized', ...
               'Color', 'k', 'FontSize', 16);
hax(2).Position(2) = 0.51;

hcb.Position(1) = hax(2).Position(1); 0.5 - pos(3)/2;
hcb.Position(2) = 0.59; 0.5 + pos(4)/2;

supax = suplabel('Along-shelf profile of cross-isobath transport at shelfbreak', 't');
supax.Position(end) = 0.77;

for ii=4:6
    hax(ii).Position(2) = 0.26;
end

export_fig -r150 -a2 images/paper2/inst-flux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x-z sections

handles = ew34.plot_xzsection(isobath, 225);
correct_ticks('y', [], {'-50','-300'}, handles.hax);
delete(handles.hrunname);
for ii=3:4
    for jj=1:2
        delete(handles.hline(ii).hl{jj});
        delete(handles.hline(ii).htxt{jj});
    end
end

for ii=1:2
    handles.hline(ii).hl{1}.delete;
    handles.hline(ii).htxt{1}.delete;
    % handles.hline(2).htxt{ii}.Units = 'normalized';
    %    handles.hline(2).htxt{ii}.Position(1) = 0.3;
end

handles.hline(1).htxt{2}.Position(1) = -150;
handles.hax(1).XLim = [-80 50];
linkprop(handles.hax, 'YLim');
handles.hax(1).YLim = [-200 0];
handles.hax(2).YTickLabelMode = 'auto';
handles.hcb(1).Ticks = sort(unique([handles.hcb(1).Ticks -0.09 0.09]));
handles.hline(1).htxt{2}.Position(1) = -70;

hanno(1) = annotation('textarrow', [0.19 0.23], [0.78 0.85], 'String', 'intrusion');
hanno(2) = annotation('textarrow', [0.63 0.67], [0.343 0.39], 'String', 'intrusion', ...
                      'Color', 'w');

export_fig -r200 -opengl -png -pdf -a2 images/paper2/ew-34-xzsection

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

export_fig('-r150', '-opengl', '-a2', ...
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
clear opt handles;
opt.addvelquiver = 0;
opt.csdcontourplot = 0;
opt.rhocontourplot = 0;
handles = ew34.secondary_vortices(320, [xx yy], opt);

handles.hax(1).XLim = [120 420];
handles.hax(1).YLim = [20 120];
handles.hfield.htlabel.Position(2) = 0.1;
hbathy = handles.hfield.hbathy;
clabel(hbathy{4}, hbathy{1}, 'LabelSpacing', 108*12);
caxis([0 280]);
handles.hfield.htrack.Visible = 'off';
handles.hfield.hcen.Color = [1 1 1]*0.8;
handles.hfield.hcb.Position(1) = 0.7;

correct_ticks('y', [], '0', handles.hax(2:end));

% displace eddy water plots
dx = 230; dpv = 4e-11;
for ii=[1 3 5 7]
    handles.hcsd(ii).XData = handles.hcsd(ii).XData + dx;
end

axes(handles.hax(1));
correct_ticks('y', [], '198');
handles.hax(1).Title.String = '(a) Cross shelf dye - X_{sb} (km)';

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

% make labels stick when I change y limits
for ii=1:length(handles.htxt(:))
    handles.htxt(ii).Units = 'data';
end

handles.hax(2).YLim = [-310 0];

export_fig -r200 -png -a2 images/paper2/ew-34-secondary-cyclone

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

%export_fig -r150 -a2 -painters -pdf -png images/paper2/sb-flux-summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg streamer profiles - shelfbreak

handles = ew34.plotAvgStreamer(1);
correct_ticks('y', [], '-50');
handles.ax(1).Title.String = '(b) Mean cross-isobath velocity (m/s)';
export_fig -r150 -a2 images/paper2/ew34-avgstreamer-sb.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg streamer profiles - offshore

handles = ew34.plotAvgStreamer(4);
handles.ax(1).Title.String = '(b) Mean cross-isobath velocity (m/s)';
correct_ticks('y', [], {'-50'; '-400'}, handles.ax(1));
handles.ax(3).XLim = [0 0.6];
export_fig -r150 -a2 images/paper2/ew34-avgstreamer-sl.png

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
export_fig -r150 -a2 images/paper2/fluxvertprofile.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shelfbreak flow surface snapshot

if ~exist('ew34', 'var')
    ew34 = runs('../topoeddy/runew-34/');
end

opt.addvelquiver = 1;
opt.rhocontourplot = 0;
opt.csdcontourplot = 0;
opt.dxi = 8; opt.dyi = 5;
opt.normquiver = 1;

tindex = 295;

ew34.read_zeta(tindex);

figure; maximize;
ax = packfig(2,1);
handles = ew34.animate_field('csdye', ax(1), tindex, 1, opt);

ylim([0 120]);
hline = linex(ew34.eddy.mx(tindex)/1000, [], 'w');

handles.hbathy{1}.ShowText = 'off';
handles.hquiv.Color = [1 1 1]*1;
handles.hquiv.LineWidth = 1.5;
handles.hquiv.AutoScaleFactor = 0.5;
handles.htlabel.Position(2) = 0.9;
handles.htlabel.Position(1) = 0.75;
handles.htlabel.Color = 'k';
handles.htlabel.FontSize = 20;
handles.htrack.delete;
handles.hbathy{3}.Color = [1 1 1]*0;
handles.hbathy{2}.Color = [1 1 1]*0;
handles.hbathy{1}.Color = [1 1 1]*0;
htitle = title('Cross shelf dye and velocity direction');
htitle.FontSize = 26;
ax(1).XTickLabel = {};
handles.hcb.Position(1) = 0.82;
handles.hcb.Position(2) = 0.55;
handles.hcb.Position(4) = 0.34;

htxt = text(ew34.eddy.mx(tindex)/1000, ew34.eddy.my(tindex)/1000 - 5, ...
            '  eddy center', 'Color', 'k', 'FontSize', 18, ...
            'HorizontalAlignment', 'center');

axes(ax(2));
dy = 0.2;
ax(2).Position(2) = ax(2).Position(2) + dy;
ax(2).Position(4) = ax(2).Position(4) - dy;

[yyax,hzeta,hcsd] = plotyy(ew34.rgrid.x_rho(1,:)/1000, ew34.zeta(:,ew34.bathy.isb,tindex)*100, ...
                           ew34.rgrid.x_rho(1,:)/1000, ew34.csdsurf(:,ew34.bathy.isb,tindex)/1000);
linkaxes([ax yyax], 'x');
xlim([150 450]);
ylabel('SSH (cm)');
xlabel('X (km)');
ax(2).YTick(end) = [];

hline2 = linex(ew34.eddy.mx(tindex)/1000, [], 'k');
correct_ticks('x', [], '213');
ylim([-0.1 0.2]);
ax(2).YTick = [-0.1 0 0.1];

yyax(2).YTick = [37.5 100 150 200];
yyax(2).YTickLabel = {'Sh'; '100'; 'Edd'; '200'};
yyax(2).YLabel.String = 'Cross-shelf dye';

yyax(1).YLabel.Color = yyax(1).YAxis.Color;
yyax(2).YLabel.Color = yyax(2).YAxis.Color;

axes(yyax(1)); beautify;
axes(yyax(2)); beautify;

% get the damned axes to line up after 'axis image'
pba = ax(1).PlotBoxAspectRatio;
pos1 = ax(1).Position;
pos2 = ax(2).Position;
pbnew = [pba(1)/pba(2)*pos1(4)/pos2(4) 1 1]
ax(2).PlotBoxAspectRatio = pbnew;
yyax(2).PlotBoxAspectRatio = pbnew;

export_fig -r150 -a2 images/paper2/sbsnapshot.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg flux - flat shelf only

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
hax(2).XTick = [0:100:600];
for ii=[2 4:9]
    hax(ii).XLim = [0 600];
    hax(ii).XTick = [100:100:600];
end
hax(7).XTick = [0 hax(7).XTick];
hax(1).XColor = hax(1).YColor;
hax(1).XTickLabelMode = 'auto';
hax(1).XAxisLocation = 'top';
hax(1).YLim = [0 120];
hax(1).YTick = [0:20:100];
hax(2).YTick = [0:20:100];

axes(hax(3));
hleg = legend;
hleg.Position(1) = 0.82;
%ylim([-1 5]);
%hax(3).YTick = [-10 0 20 40];
correct_ticks('y', '%.2f', '0.20');
hax(3).Position(1) = 0.68;

export_fig -r150 -a2 -png -pdf -opengl images/paper2/avgflux

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
%% 3d schematics
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

mergename = 'images/paper2/3dmerged.png';
multifig = 1; % multiple figure windows and stitch later?
day = [220 230 240 260];
depth = 75;

annocolor = 'k'; %[1 1 1]*0.5;
annofs = 23;
annolw = 0.5;
annoheadstyle = 'none'; 'cback3';
annofontname = 'Futura Bk BT';

zoomin = 1; % for presentation
opt.MoveToZLevel = -850;
opt.eddthresh = 0.8;
opt.csdcontours = ew34.bathy.xsb+5000;
opt.eddreducepatch = 0.3;
opt.csdreducepatch = 0.3;
opt.finalize = 1;
opt.linefilter = 1;
opt.nolabels = 1;

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

    hax(ii).DataAspectRatio = [1 1 5];
end

linkprop(hax, 'DataAspectRatio');
linkprop([handles.hbathy], {'FaceColor', 'FaceAlpha'});
linkprop([handles.hedd], {'FaceColor', 'AmbientStrength'});
linkprop([handles.hcsd], {'FaceColor', 'AmbientStrength'});
linkprop(hax, 'ZLim');

% change views and lighting as needed
%axes(hax(1));
%view(-125, 28);

axes(hax(1));
view(-110,28);
hanno(1) = annotation(hax(1).Parent, 'textarrow', [0.64 0.56], [0.52 0.62], ...
                      'String', {'The cyclonic';  'wave is'; 'a wrinkle'; 'in 3-D.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'top');
hlabel(1) = annotation(hax(1).Parent, 'textbox', [0.35 0.18 0.05 0.05], ...
                       'String', '(a)', 'EdgeColor', 'none');

axes(hax(2));
view(-122, 42);
handles(2).hlight.Position = [-1400 1400 500];
hanno(2) = annotation(hax(2).Parent, 'textarrow', [0.645 0.53], [0.49 0.59], ...
                      'String', {'The wave';  'propagates'; ...
                    'around the eddy'; 'with the streamer.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hlabel(2) = annotation(hax(2).Parent, 'textbox', [0.4 0.18 0.05 0.05], ...
                       'String', '(b)', 'EdgeColor', 'none');

axes(hax(3));
view(-130, 28);
hanno(3) = annotation(hax(3).Parent, 'textarrow', [0.395 0.47], [0.84 0.75], ...
                      'String', {'The wave rolls up into'; ...
                    'a cyclone trapping'; 'the shelf-slope water'; 'above it.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hanno(4) = annotation(hax(3).Parent, 'textarrow', [0.67 0.57], [0.53 0.61], ...
                      'String', {'There is a';  'persistent bulge'; ...
                    'in the eddy'; 'below'; 'shelfbreak depth'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hline = annotation(hax(3).Parent, 'arrow', [1 1]*0.48, [0.93 0.8], 'LineWidth', annolw, ...
                   'HeadStyle', annoheadstyle);
hanno(7) = annotation(hax(3).Parent, 'textarrow', [0.6 hline.X(1)], [1 1]*hline.Y(1), ...
                      'String', {'Phase shift between'; 'cyclonic and'; 'anti-cyclonic waves'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', 'none', 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');

hlabel(3) = annotation(hax(3).Parent, 'textbox', [0.4 0.13 0.05 0.05], ...
                       'String', '(c)', 'EdgeColor', 'none');

axes(hax(4));
view(-128, 40);
hanno(5) = annotation(hax(4).Parent, 'textarrow', [0.41 0.45], [1 1]*0.80, ...
                      'String', {'Cyclone propagates'; ...
                    'away with trapped';  'streamer water.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs, ...
                      'VerticalAlignment', 'middle');
hanno(6) = annotation(hax(4).Parent, 'textarrow', [0.66 0.575], [0.53 0.61], ...
                      'String', {'The process'; 'repeats.'}, ...
                      'LineWidth', annolw, 'Color', annocolor, ...
                      'HeadStyle', annoheadstyle, 'FontSize', annofs);
hlabel(4) = annotation(hax(4).Parent, 'textbox', [0.4 0.18 0.05 0.05], ...
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
        export_fig('-r150', '-a4', '-p0.01', '-opengl', ...
                   ['images/paper2/3d-schem-' num2str(ii) '.png']);
    end
    hash = githash;
    system(['montage images/paper2/3d-schem-[1-4].png ' ...
            '-geometry +0.05+0.05 ' mergename]);
    system(['exiftool -overwrite_original -Producer=' hash ' ' mergename]);
else
    export_fig('-opengl', '-r150', '-a4', '-p0.02', mergename);
end

if zoomin
    axes(hax(3));
    hax(3).XLim = [200 350];
    hax(3).YLim = [0 150];
    hax(3).ZLim = [-300 1.1];
    hanno(3).X = [0.3210 0.3784];
    hanno(3).Y = [0.7302 0.6604];
    hanno(4).delete;
    hanno(7).X = [0.60 0.4];
    hanno(7).Y = [1 1]*0.90;
    hline.X = [1 1]*hanno(7).X(2);
    hline.Y = [hanno(7).Y(2) 0.73];
    handles(3).hsect.delete
    handles(3).hsectoutline.delete;
    handles(3).hplane.delete;
    export_fig('-r150', '-a4', '-p0.01', '-opengl', ...
               ['images/paper2/3d-schem-' num2str(3) '-zoom.png']);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% churchill figure
if ~exist('ew36', 'var')
    ew36 = runs('../topoeddy/runew-36/');
end

figure; maximize;
hax(1) = subplot(1,5,[1 2]);
hax(2) = [];
hax(3) = subplot(1,5,[3]);
hax(4) = subplot(1,5,[4 5]);
handles = PlotOleanderSection(114, [17 11], [], hax(1:3));
hax(2) = handles.hax(2);
hax(2).Title.String = '(a) Oleander XBT Temperature (C)';
hax(3).XLabel.String = '';
handles.hleg.delete;

axes(hax(1));
hax(1).XLim = [550 905];
linkaxes(hax(1:2));
handles.hax(1).YLabel.String = 'Depth (m)';
correct_ticks('x', [], {'650'; '750'; '700'}, handles.hax(1));
htxt = text(0.07, 0.1, {'Oleander XBT'; '19-20 Oct, 1983'}, 'Units', 'normalized');
handles.hax(2).Title.Position(1) = 600;
cmap = cbrewer('seq','Reds',50);
colormap(hax(1), cmap(1:end-12,:))
caxis([10 20]);

axes(hax(4));
handles = ew36.PlotSingleYZSection('csdye', '158.5', [], hax(4));
handles.hed.LevelList = 135;
handles.hed.Color = 'k';
handles.htlabel.Position = [0.02 0.07 0];
hax(4).DataAspectRatio = [1 3 1];
ylabel('');
ylim([-450 0]);
xlim([0 150]);
colorbar('delete');
hcb = colorbar('SouthOutside');
title('(b) Cross-shelf dye (km)');
hax(4).YTick = sort(unique([hax(4).YTick -1*(0:50:400)]));
hax(4).YTickLabelMode = 'auto';
hax(4).XTick = sort(unique([hax(4).XTick 0 50 100 150]));
hax(4).XTickLabelMode = 'auto';
correct_ticks('y', [], '-200');
beautify;

export_fig -r150 -a2 images/paper2/eddy-intrusion.png