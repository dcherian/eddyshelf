%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% field map

if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end
factor = 2;
isobath = 4;
timesteps = [1 150 200 250 300 380];

opt.csdcontourplot = 1;
opt.csdcontours = ew34.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.zetaOnFirstPlot = 1;
handles = ew34.mosaic_field('csdye', timesteps, opt);
handles.hfield{1}.hzeta.LevelList = handles.hfield{1}.hzeta.LevelList(1:2:end);
handles.hcb.delete;

ylim([0 250]);
correct_ticks('y', [], {'50', '100'}, handles.hax([1 4]));

handles.supax.Position(4) = 0.70;
handles.htitle.String = 'Surface cross-shelf dye (km)';

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

export_fig -a1 images/paper2/ew-34-surface-csdye.png

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
export_fig -a1 -png -pdf images/paper2/centracks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% field map

opt.csdcontourplot = 1;
opt.csdcontours = ew34.csflux.x([1 4 8]);
opt.addvelquiver = 0;
opt.zetaOnFirstPlot = 1;
handles = ew34.mosaic_field('csdye', [1 150 200 250 300 380], opt);
ylim([0 250]);
handles.supax.Position(4) = 0.70;
handles.htitle.String = 'Surface cross-shelf dye (km)';
handles.hcb.delete;
correct_ticks('y', [], {'50', '100'}, handles.hax([1 4]));
export_fig -a1 images/paper2/ew-34-surface-csdye.png

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
handles.icons(10).FaceColor = color;
handles.icons(10).FaceAlpha = 0.2;

axes(handles.hax(2)); % make visible
export_fig -a3 images/paper2/flux-diags.png

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

supax = suplabel('Along-shelf profile of cross-isobath velocity at shelfbreak', 't');
supax.Position(end) = 0.79;

export_fig -a1 images/paper2/inst-flux.png

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
drawnow;

export_fig -a1 images/paper2/ew-34-xzsection.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% secondary eddy - 2360

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

export_fig -a1 images/paper2/ew-34-secondary-cyclone.png

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
export_fig -a1 images/paper2/ew-34-mosaic-zslice.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ew-2360 z-slice mosaic

handles = ew2360.mosaic_zslice('dye_03', 200, [63 70 77 95]);
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
handles(1).supax.Position(4) = 0.87;
axes(handles(1).hax); correct_ticks('y', '', {'150' '200'});
axes(handles(3).hax); correct_ticks('y', '', {'150' '200'});

export_fig -a1 images/paper2/ew-2360-mosaic-zslice.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom rho PV

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

export_fig -a1 images/paper2/ew-2360-bottomrhopv.png

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
        'S_{sh} = 0 | r = 0'; ...
        'S_{sh} = 0 | r = 5e-3'; ...
        'S_{sh} = 0.05 | r = 0'; ...
        'S_{sh} = 0.05 | r = 5e-4'; ...
        'S_{sh} = 0.05 | r = 5e-3'; ...
            };
    shfric = runArray(folders, names);
end

handles = shfric.PlotFluxSummary(1);
axes(handles.hax(1))
ylim([0 15]);
delete(handles.hleg);
handles.hleg(1) = legend(handles.hflux(1:2), shfric.name{1:2});
handles.hleg(2) = legend(handles.henv(3:5), shfric.name{3:5});
handles.hleg(1).Box = 'off';
handles.hleg(2).Box = 'off';
handles.hleg(1).Position(1) = 0.62;
pos1 = handles.hleg(1).Position;
handles.hleg(2).Position(1) = 0.62 + pos1(3);
handles.hleg(2).Position(2) = pos1(2) + pos1(4) - handles.hleg(2).Position(4);
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

dy = dy*0;
for ii= 5
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

export_fig -a3 images/paper2/sb-flux-summary.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg streamer profiles - shelfbreak

handles = ew34.plotAvgStreamer(1);
correct_ticks('y', [], '-50');
handles.ax(1).Title.String = 'Mean cross-isobath velocity (m/s)';
export_fig -a1 images/paper2/ew34-avgstreamer-sb.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg streamer profiles - offshore

handles = ew34.plotAvgStreamer(4);
handles.ax(1).Title.String = 'Mean cross-isobath velocity (m/s)';
correct_ticks('y', [], {'-50'; '-400'}, handles.ax([1 3]));
export_fig -a1 images/paper2/ew34-avgstreamer-sl.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% flux vertical profiles

clear handles;
if ~exist('fluxvert', 'var')
    fluxvert = runArray({'runew-2360_wider', 'runew-34-cyc'});
end

day = {'119'; []; '75'};
names = {'Anticyclone'; 'Cyclone'}; %; 'Anticyclone - flat bottom'};
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
        fluxvert.array(ii).PlotSingleXZSection('v', 4, day{ii}, hax(ii));
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
export_fig -a3 images/paper2/fluxvertprofile.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shelfbreak flow surface snapshot

if ~exist('ew04', 'var')
    ew04 = runArray({'runew-04', 'runew-8041'});
end

opt.addvelquiver = 1;
opt.rhocontourplot = 1;
opt.csdcontourplot = 0;

handles = ew04.mosaic_field('csdye', {'169'; '169'}, opt);
xlim([200 450]);
ylim([0 150]);
for ii=1:2
    delete(handles(ii).hrunname)
    handles(ii).hbathy{1}.ShowText = 'off';
    handles(ii).hquiv.Color = [1 1 1]*0.65;
    handles(ii).hquiv.LineWidth = 2;
end
correct_ticks('x', [], '450', handles(1).hax);
handles(1).hax.Title.String = 'Flat shelf | S_{sh} = 0';
handles(2).hax.Title.String = 'Flat shelf | S_{sh} = 0.05';

handles(1).supax.Position(4) = 0.76;
handles(1).supax.Title.String = 'Surface cross-shelf dye (km)';
handles(1).supax.Title.FontSize = 20;

export_fig -a1 images/paper2/sbsnapshot.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% avg flux

hax = csf.plot_fluxparam('avg flux');
axes(hax(3));
hleg = legend;
hleg.Position(1) = 0.82;

export_fig -a2 images/paper2/avgflux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max flux

hax = csf.plot_fluxparam('max flux');
axes(hax(3));
hleg = legend;
hleg.Position(1) = 0.82;
ylim([-1 1]*55);

export_fig -a2 images/paper2/maxflux.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3d schematics
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end

multifig = 1; % multiple figure windows and stitch later?
day = [210 220 240 260];
depth = 75;
MoveToZLevel = -1100;

opt.eddthresh = 0.8;
opt.csdcontours = ew34.bathy.xsb+5000;
opt.eddreducepatch = 0.3;
opt.csdreducepatch = 0.3;
opt.finalize = 1;
opt.linefilter = 1;

clear handles hslice hplane
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

    hslice(ii) = ew34.animate_zslice('dye_03', depth, day(ii), hax(ii), opt);
    hslice(ii).bathy{1}.delete;
    hslice(ii).bathy{2}.delete;
    hslice(ii).bathy{3}.delete;
    hslice(ii).eddsurf.delete;
    hslice(ii).htime.delete;
    colorbar('delete');
    linkprop([hslice(ii).hvar hslice(ii).csdsurf hslice(ii).rhocont], 'ContourZLevel');
    hslice(ii).hvar.ContourZLevel = MoveToZLevel;
    hslice(ii).hvar.ZData = fillnan(double(hslice(ii).hvar.ZData > 0.4), 0);

    handles(ii) = ew34.animate_3d(num2str(day(ii)), opt, hax(ii));
    handles(ii).hcsd.FaceAlpha = 1;

    [xmat, ymat] = meshgrid(xlim, ylim);
    hplane(ii) = surf(xmat, ymat, -1*abs(depth) * ones(size(xmat)), ...
                      'FaceColor', 'k', 'FaceAlpha', 0.05, 'EdgeColor', [1 1 1]*0.3);

    houtline(ii) = plot3([xmat(:); xmat(1)], [ymat(:,1)' flip(ymat(:,2)') ymat(1)], ...
                         MoveToZLevel * ones([1 5]), ...
                         '-', 'LineWidth', 1, 'Color', [1 1 1]*0.3, 'Tag', 'dcline');

    hax(ii).Title.String = '';
    hax(ii).XAxis.Visible = 'off';
    hax(ii).YAxis.Visible = 'off';
    hax(ii).ZAxis.Visible = 'off';
    hax(ii).Box = 'off';
    hax(ii).DataAspectRatio = [1 2.5 8];

    %bathy = handles(ii).hbathy.ZData;
    %handles(ii).hbathy.ZData(bathy < -700) = -700;
    handles(ii).hcsd.FaceColor = [107 174 214]/255;
    handles(ii).hedd.FaceColor = brighten([215 48 31]/255, 0.1);
end

linkprop(houtline, {'LineWidth', 'Color'});
linkprop(hplane, {'FaceAlpha', 'EdgeColor'});

linkprop(hax, 'DataAspectRatio');
linkprop([handles.hbathy], {'FaceColor', 'FaceAlpha'});
linkprop([handles.hedd], 'FaceColor');
linkprop([handles.hcsd], 'FaceColor');

axes(hax(1));
zlim([-1250 0]);
linkprop(hax, 'ZLim');

% change views and lighting as needed
axes(hax(1));
view(-125, 28);

axes(hax(2));
view(-110,28);

axes(hax(3));
view(-130, 28);

axes(hax(4));
view(-130, 28);

tic;
if multifig
    for ii=1:4
        axes(hax(ii));
        export_fig('-r96', '-a4', '-p0.02', '-opengl', ...
                   ['images/paper2/3d-schem-' num2str(ii) '.png']);
    end
    system(['montage images/paper2/3d-schem-[1-4].png ' ...
            '-geometry +10+10 ' ...
            'images/paper2/3dmerged.png']);
else
    export_fig -opengl -r96 -a4 -p0.02 images/paper2/3d-schem.png;
end
toc;
