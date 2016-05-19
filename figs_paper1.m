%% Figures for paper 1:

%% ew-64461-5 dye
if ~exist('ew','var') | ~strcmpi(ew.name, 'ew-64461-5')
    ew = runs('../topoeddy/runew-64461-5/');
end
opt = [];
opt.addcsdye = 1;
opt.addzeta = 1;
opt.rhocontourplot = 1;
tsteps = [1 40 230];

figure;
for ii=1:6, hax(ii) = subplot(2,3,ii); end
maximize;

handles = ew.mosaic_field('eddye', tsteps, opt, hax(1:3));
handles.hcb.delete;
handles.hfield{2}.hshelflabel.delete;
handles.hfield{3}.hshelflabel.delete;
handles.supax.Position(4) = 0.65;
handles.htitle.String = 'Surface dyes';

correct_ticks('y', [], {'200'; '50'}, handles.hax(1));
for ii=1:3
    str = handles.hfield{ii}.htlabel.String;
    handles.hfield{ii}.htlabel.String = str(1:2);
    hax(ii).Title.String = str(4:end);
    handles.hfield{ii}.hbathy{2}.Color = [1 1 1]*0.77;
    if ii > 1
        handles.hfield{ii}.hzetaneg.LevelList = ...
            [-6 -4 -3 -2]*1e-3;
    end
end

labels = 'def';
for ii=1:3
    handlesyz{ii} = ew.PlotSingleYZSection('rho', tsteps(ii), [], hax(3+ii));
    title('');
    if ii~=3
        colorbar('delete');
        correct_ticks('x', [], {'200'}, hax(3+ii));
    else
        hcb = colorbar;
        hcb.Position(1) = 0.91;
        hcb.Label.String = '\rho (kg/m^3)';
    end
    handlesyz{ii}.htlabel.String{1} = [labels(ii) ') '];
    handlesyz{ii}.htlabel.String{2} = [];

    handlesyz{ii}.hlzlabel.Position(2) = ...
        handlesyz{ii}.hlzlabel.Position(2) - 75;
    handlesyz{ii}.hlzlabel.Position(1) = ...
        handlesyz{ii}.hlzlabel.Position(1) - 20;

    caxis([21.6 22.4]);

    % change colormap
    colormap(hax(3+ii), cbrewer('div', 'BrBG', 20));
    if ii > 1
        hax(3+ii).YTickLabel = [];
        hax(3+ii).YLabel.String = '';
    end
end

correct_ticks('y', [], {'-300'}, hax(4));
% move panels closer to each other
hax(1).Position(1) = hax(1).Position(1) + 0.05;
hax(4).Position(1) = hax(4).Position(1) + 0.05;
hax(3).Position(1) = hax(3).Position(1) - 0.05;
hax(6).Position(1) = hax(6).Position(1) - 0.05;
hcb.Position(1) = hcb.Position(1) - 0.05;
for ii=4:6
    hax(ii).Position(2) = hax(ii).Position(2) + 0.06;
end
hcb.Position(2) = hcb.Position(2) + 0.06;

export_fig -opengl -r150 -a2 images/paper1/xyzmap.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ew-64461-5-eddye video
if ~exist('ew','var') | ~strcmpi(ew.name, 'ew-64461-5')
    ew = runs('../topoeddy/runew-64461-5/');
end
opt = [];
opt.addcsdye = 1;
opt.addzeta = 1;
ew.makeVideo = 1;
ew.animate_field('eddye', [], 1, [], opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ew-64461-5-eddye eddydiag
if ~exist('ew','var') | ~strcmpi(ew.name, 'ew-64461-5')
    ew = runs('../topoeddy/runew-64461-5/');
end

handles = ew.EddyDiags(gcf);

correct_ticks('y', [], {'0.5'}, handles.hax(2));

export_fig -r150 -a2 images/paper1/eddydiag.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EW isobaths
ewall = runArray({ ...
    'runew-4341/', ...
    'runew-36-sym/', ...
    'runew-64361/', ...
                });
ewall.name = {'R/L_{sl} > 1', 'R/L_{sl} ~ 1', ...
              'R/L_{sl} < 1'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EW, NS isobaths
clear ew
ew = runArray({ ...
    'runew-64361-shallow', ...
    'runew-6341', ...
    'runew-6362-2', ...
    'runew-6441' });

% NS isobaths
folders = { ...
    'runns-64361', 'runns-6341', 'runns-6362-2',...
    'runns-6441', ...
          };
ns = runArray(folders);

ns.print_diag('params table', 'images/paper1/sl-params-table-ns.org');

ew.name{1} = 'ew-64361';
ew.sort(ew.print_params('params.nondim.eddy.Rh'))
ns.sort(ns.print_params('params.nondim.eddy.Rh'))

for ii=1:ew.len
    ns.name{ii} = [num2str(ns.array(ii).params.nondim.eddy.Rh', '%2.0f') ...
                   ' | ' ns.name{ii}];
    ew.name{ii} = [num2str(ew.array(ii).params.nondim.eddy.Rh', '%2.0f') ...
                   ' | ' ew.name{ii}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shelfbreak depth - xy maps
sb = runArray({ ...
    'runew-36', 'runew-2360_wider', 'runew-2361_wider', ...
    'runew-2363_wider', 'runew-2362_wider', ...
              });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wide slope runs
% folders = { ...
%     'runew-6341', 'runew-56341', 'runew-56341-2', ...
%     'runew-6341-2', 'runew-6341-3', 'runew-6342', ...
%     'runew-6362', ... %'runew-6362-1', 'runew-6362-2', 'runew-6371', ...
%     'runew-6372', ... % 'runew-6373' ...
%     ...                     %'runew-6452', 'runew-6352', ...
%     'runew-64361', 'runew-64361-2', ...
%     ... %    'runew-564361', 'runew-564361-1', 'runew-564361-2', 'runew-564361-3', ...
%     'runew-64461-3', ...%'runew-64462', ...
%     'runew-64461-4', 'runew-64461-5',...
%     'runew-64461-6', 'runew-64461-7', 'runew-64463-1', ...
%     ... %    'runew-64351-cyc', 'runew-6341-fneg', ...
%           };

folders = { ...
    'runew-6341', 'runew-6341-2', 'runew-6341-3', 'runew-6341-4', ...
    'runew-56341', 'runew-56341-2', ...
    'runew-6342', 'runew-6352', ...
    'runew-6362', 'runew-6362-1', 'runew-6362-2', ...
    'runew-6371', 'runew-6372', ...% 'runew-6373' ...
    ... % 'runew-6441', ... % high Rh!
    'runew-64462', ...
    'runew-64461-4', 'runew-64461-5', ...
    'runew-64461-6', 'runew-64461-7', ...
    'runew-64463-1', ...
    ... %'runew-b4361', ...
    ... %'runew-64351-cyc',
    'runew-6341-fneg', ...
    ... % 'runew-64361', moved to external
    'runew-64361-2', ...
    ... %'runew-564361', 'runew-564361-2', 'runew-564361-3', ...
    ... %'runew-64361-deep', 'runew-64361-shallow', 'runew-64361-sl' ...
    'runew-64461-3', 'runew-5644613-1/', ...
    'runew-64461-9-shallow', 'runew-64461-9-deep', ...
    'runew-64461-8', 'runew-64461-8-deep', 'runew-64461-10', 'runew-64461-11', ...
    ... %'runew-b371', ...
          };
sl = runArray(folders);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom friction
folders = { ...
    'runew-6341', 'runew-56341', 'runew-56341-2', ...
          };
bfrics = runArray(folders);
for ii=1:bfrics.len
    run = bfrics.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    bfrics.name{ii} = ['r = ' ...
                       num2str(bfrics.array(ii).params.misc.rdrg, '%1.0e') ' m/s'];
end
bfrics.name{1} = 'r = 0 m/s';

folders = { ...
    'runew-64361-shallow', 'runew-564361', 'runew-564361-2', ...
    'runew-564361-3', 'runew-564361-4', ...
          };
bfrics60 = runArray(folders);
for ii=1:bfrics60.len
    run = bfrics60.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    bfrics60.name{ii} = ['r = ' ...
                       num2str(run.params.misc.rdrg, '%1.0e') ' m/s'];
end
bfrics60.name{1} = 'r = 0 m/s';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% y-z cross section
name = 'ew-34';
if ~exist('run', 'var') || ~strcmpi(run.name, name)
    run = runs(['../topoeddy/run' name '/']);
end

tt = [1 250];
run.plot_eddye(tt);
subplot(121); correct_ticks('y',[],5);
suplabel('', 't');

export_fig('-r450', 'images/paper1/yzsection.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EW - center tracks - all,
fs  = 18;
ewall.plot_penetration(gca,'all');
maximize();
ax = gca;
title([]); pbaspect([1.618 1 1]);
text(-0.3, 1.5, 'R > L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(1,:));
text(0.2, 3.8, 'R ~ L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(2,:));
%text(-9, 2.1, 'R < L_{sl}', ...
%     'FontSize', fs, 'Color', ax.ColorOrder(3,:));
text(-6, 3.8, 'R < L_{sl}', ...
     'FontSize', fs, 'Color', ax.ColorOrder(3,:));
%legend('off');
export_fig('images/paper1/centrack.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EW, NS Center-tracks - wide slope
ew.filter = [];
ns.filter = [];

figure; maximize();
ax1 = subplot(121);
hew = ew.plot_penetration(ax1, 'all');
hew.hgplt2(1).Color = 'k';
linkprop([hew.hgplt2 hew.hfit hew.hslbreak], 'Color');
hew.hgplt2(1).LineStyle = '--';
hew.hgplt2(2).LineStyle = '--';
ax1.Title.String = ['a) ' ax1.Title.String];
beautify

ax2 = subplot(122);
hns = ns.plot_penetration(ax2, 'all'); drawnow;
beautify
hns.hgplt2(1).Color = 'k';
linkprop([hns.hgplt2 hns.hfit hns.hslbreak], 'Color');
hns.hgplt2(1).LineStyle = '--';
hns.hgplt2(2).LineStyle = '--';
ax2.Title.String = ['b) ' ax2.Title.String];

ax2.XTick = unique([ax2.XTick 1]);
ax1.XLim = [-15 0];
hline = findall(ax1,'Tag','dcline');
hline.XData = ax1.XLim;
export_fig('-r150', '-a2', 'images/paper1/sl-centrack.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom friction

bfrics.sort(bfrics.print_params('params.misc.rdrg'));
bfrics60.sort(bfrics60.print_params('params.misc.rdrg'));

figure; maximize
hax = packfig(1,2);
hfric = bfrics.plot_penetration(hax(1),'all'); maximize();
pbaspect([1.618 1 1])
hfric.hgplt2(1).Color = [1 1 1] * 0.7;
hfric.hgplt2(2).Color = [1 1 1] * 0.4;
hfric.hgplt2(3).Color = [1 1 1] * 0.2;

hfric60 = bfrics60.plot_penetration(hax(2),'all',1); maximize();
hax(2) = gca;
pbaspect([1.618 1 1])
hfric60.hgplt2(1).Color = [1 1 1] * 0.8;
hfric60.hgplt2(2).Color = [1 1 1] * 0.6;
hfric60.hgplt2(3).Color = [1 1 1] * 0.4;
hfric60.hgplt2(4).Color = [1 1 1] * 0.2;

linkaxes(hax, 'xy');
hax(1).YLim = [0 8];
hax(1).YTick = [0:1:8];
hax(2).YLabel.String = [];
hax(2).YTickLabel = {};
hax(1).XTickLabel{end} = '';

hax(1).Title.String = 'a) Rh = 12';
hax(2).Title.String = 'b) Rh = 60';

export_fig -r300 -a2 images/paper1/bfrics-centrack.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameterization
[~,~,~,~,~,~,~,~,handles] = ...
    sl.print_diag('bottom torque', 'regression', [], 'no_name_points');
title([]);

hax = gca; hleg = legend;
hax.FontSize = 30
handles.htext.FontSize = 28;
hleg.FontSize = handles.htext.FontSize;
%hax.XAxis.Axle.VertexData(1,:) = single([0.0051 0.1139]);
%hax.YAxis.Axle.VertexData(2,:) = single([0.0132 0.1768]);

export_fig -r200 -a2 images/paper1/penetration-param.png

sl.print_diag('params table', 'images/paper1/sl-params-table.org');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functional form
[~,~,~,~,~,~,~,~,handles] = ...
    sl.print_diag('bottom torque', 'functional', [], 'no_name_points');
title([]);

hax = gca; hleg = legend;
hax.FontSize = 30
hax.XTick(1) = 5e-3;
hax.YTick(1) = [];
hax.YTick(end) = 1.92;
%hax.XAxis.Axle.VertexData(1,:) = single([0.0051 0.1139]);
%hax.YAxis.Axle.VertexData(2,:) = single([0.0132 0.1768]);

export_fig -r200 -a2 images/paper1/functional.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% energy decay
hlines = sl.plot_dEdt; maximize(); pause(1);
prop = linkprop([hlines.hke hlines.hpe], 'Color');
hlines.hke(1).Color = [1 1 1]*0.55;
a = strcmp(sl.name, 'ew-64461-5');
prop.delete;
hlines.hke(a).Color = [0 0 0];
hlines.hpe(a).Color = [0 0 0];
uistack(hlines.hke(a), 'top');
uistack(hlines.hpe(a), 'top');
subplot(121); title([]);
pbaspect([1.618 1 1]);
xlim([0 200]);
legend('off');
subplot(122);
xlim([0 200]);
pbaspect([1.618 1 1]);
export_fig('-r150', '-a2', 'images/paper1/energy-decay.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shallow and deep

if ~exist('ewshdp', 'var')
    ewshdp = runArray({'ew-64461-9-shallow', 'ew-64461-9-deep', ...
                       'ew-64361-shallow', 'ew-64361-deep'});
end

fontsize = [22 24 28];
figure; maximize;
ax(1) = subplot(121);
ewshdp.filter = [1 2];
handles = ewshdp.plot_ts('eddy.hcen', ax(1))
handles.hplt(1).Color = 'k';
handles.hplt(2).Color = 'k';
handles.hplt(2).LineStyle = '--';
pbaspect([1.618 1 1]);
title('$$a) \frac{U}{\beta L^2} \sim 20 $$');
ax(1).Title.Interpreter = 'latex';
legend('off');
beautify(fontsize);

ax(2) = subplot(122);
ewshdp.filter = [3 4];
handles = ewshdp.plot_ts('eddy.hcen', ax(2));
handles.hplt(1).Color = 'k';
handles.hplt(2).Color = 'k';
handles.hplt(2).LineStyle = '--';
pbaspect([1.618 1 1]);
title('$$b) \frac{U}{\beta L^2} \sim 60 $$');
ax(2).Title.Interpreter = 'latex';
legend('off');
beautify(fontsize);

linkaxes(ax, 'y');
packfig(1,2, 'columns');
ax(2).XTick(1) = [];

[ax(3), htitle] = suplabel('Water depth at eddy center', 't');
ax(3).Position(4) = 0.80;
ax(3).FontSize = fontsize(end);
beautify(fontsize);

export_fig -r150 -a2 images/paper1/shallow-deep-hcen.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% high, low Rh

rh = runArray({'ew-64361-shallow', 'ew-64361-deep', ...
               'ew-6441', 'ew-564361-2', 'runew-564361-3', ...
               'ew-6040', 'ew-6041', 'ew-6042', ...
               'ew-6050', 'ew-6051'})
rh.print_diag('params table', 'images/paper1/sl-params-table-Rh.org');