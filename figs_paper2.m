% describe fluxes
if ~exist('ew34', 'var') | ~strcmpi(ew34.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end
factor = 2;
isobath = 4;
timesteps = [1 150 200 250 300 380];

%% field map
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
set(gcf, 'Position', [675 175 1080 924]);
columnlegend(3, names, 'Location', 'NorthWest');
export_fig -a1 -png -pdf images/paper2/centracks

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

%% surface map and instantaneous flux
folders = {'ew-4341', 'ew-34', 'ew-2341_wider'};
if ~exist('ewflux', 'var'), ewflux = runArray(folders); end
N = length(folders);
times = [180 200 250];
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
    axes(hax(ii));
    xlabel(''); set(gca, 'XTickLabel', []); colorbar('off');
    axes(hax(N+ii));
    title('');

    if ii ~= 1
        axes(hax(ii)); ylabel('');
        axes(hax(N+ii)); ylabel('');
    end
end

axes(hax(4));
ylabel('\int v(x,z) dz');
axes(hax(N));
hcb = colorbar('southoutside');
pos = hcb.Position;
hcb.Position(1) = 0.5 - pos(3)/2;
hcb.Position(2) = 0.5 + pos(4)/2;
hcb.Label.String = 'Cross shelf dye - X_{sb} (km)';

export_fig -a1 images/paper2/inst-flux.png

%% x-z sections
handles = ew34.plot_xzsection(isobath, 225);
correct_ticks('y', [], 3, handles.hax([1 3 4]));
correct_ticks('y', [], 4, handles.hax([2]));
drawnow;
export_fig -a1 images/paper2/ew-34-xzsection.png

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

%% z-slice
opt.csdcontours = ew34.csflux.x(isobath);
handles = ew34.mosaic_zslice('dye_03', 100, [225 230 235 240 245 260], opt);
for ii=1:length(handles)
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([200 400]);
ylim([0 150]);
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.10';
handles(1).htitle.FontSize = 22;
handles(1).htitle.FontWeight = 'normal';
handles(1).supax.Position(4) = 0.87;
axes(handles(4).hax); correct_ticks('x', '', '400');
axes(handles(5).hax); correct_ticks('x', '', '400');
export_fig -a1 images/paper2/ew-34-mosaic-zslice.png

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

%% bottom rho PV

hplt = ew2360.plotBottomRhoPV(70);
title('Bottom PV with \rho contours');
xlabel('X (km)');
ylabel('Y (km)');
text(330, 149, 'shelfbreak')
ylim([140 180]);
set(gca, 'YTickMode', 'auto');

export_fig -a1 images/paper2/ew-2360-bottomrhopv.png

%% secondary cyclone x-z section

%% shelf + friction flux time series

%% shfric
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

%% flux summary
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

%% avg streamer profiles - shelfbreak
handles = ew34.plotAvgStreamer(1);
correct_ticks('y', [], '-50');
handles.ax(1).Title.String = 'Mean cross-isobath velocity (m/s)';
export_fig -a1 images/paper2/ew34-avgstreamer-sb.png

%% avg streamer profiles - offshore
handles = ew34.plotAvgStreamer(4);
handles.ax(1).Title.String = 'Mean cross-isobath velocity (m/s)';
correct_ticks('y', [], {'-50'; '-400'}, handles.ax([1 3]));
export_fig -a1 images/paper2/ew34-avgstreamer-sl.png