%% Figures for paper 1:
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
%% EW isobaths
ewall = runArray({ ...
    'runew-4341/', ...
    'runew-36-sym/', ...
    'runew-64361/', ...
                });
ewall.name = {'R/L_{sl} > 1', 'R/L_{sl} ~ 1', ...
              'R/L_{sl} < 1'};

ew = runArray({ ...
    'runew-64361-shallow', ...
    'runew-6341', ...
    'runew-6362-2', ...
    'runew-6441' });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NS isobaths
folders = { ...
    'runns-64361', 'runns-6341', 'runns-6362-2',...
    'runns-6441', ...
          };
ns = runArray(folders);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shelfbreak depth - xy maps
sb = runArray({ ...
    'runew-36', 'runew-2360_wider', 'runew-2361_wider', ...
    'runew-2363_wider', 'runew-2362_wider', ...
              });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wide slope runs
folders = { ...
    'runew-6341-2', 'runew-6341-3', 'runew-6342', ...
    'runew-6362', ... %'runew-6362-1', 'runew-6362-2', 'runew-6371', ...
    'runew-6372', ... % 'runew-6373' ...
    ...                     %'runew-6452', 'runew-6352', ...
    'runew-64361', 'runew-64361-2', ...
... %    'runew-564361', 'runew-564361-1', 'runew-564361-2', 'runew-564361-3', ...
    'runew-64461-3', ...%'runew-64462', ...
    'runew-64461-4', 'runew-64461-5',...
    'runew-64461-6', 'runew-64461-7', 'runew-64463-1', ...
... %    'runew-64351-cyc', 'runew-6341-fneg', ...
          };
sl = runArray(folders);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom friction
folders = { ...
    ... %'runew-34', 'runew-5341', 'runew-5343'
    'runew-6341', 'runew-56341', 'runew-56341-2', ...
          };
bfrics = runArray(folders);
for ii=1:bfrics.len
    run = bfrics.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    bfrics.name{ii} = ['r = ' ...
                       num2str(bfrics.array(ii).params.misc.rdrg) ' m/s'];
end

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

fontSize = [22 22 30]
figure; maximize();
ax1 = subplot(121);
ew.plot_penetration(ax1, 'all');
beautify(fontSize)
legend('off');
subplot(122);
ns.plot_penetration(gca, 'all'); drawnow;
beautify(fontSize)
legend('off');
ax1 = gca; ax1.XTick = unique([ax1.XTick 1])
export_fig('-r150','images/paper1/sl-centrack.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bottom friction

bfrics.plot_penetration([],'all'); maximize();
pbaspect([1.618 1 1]);
export_fig('images/paper1/bfrics-centrack.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameterization
sl.print_diag('bottom torque');
title([]); pause(1);
export_fig('images/paper1/penetration-param.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% energy decay
sl.filter = [3 4 5 7 8 12]
sl.plot_dEdt; maximize(); pause(1);
subplot(121); title([]);
pbaspect([1.618 1 1]);
xlim([0 200]);
legend('off');
subplot(122);
xlim([0 200]);
pbaspect([1.618 1 1]);
export_fig('images/paper1/energy-decay.pdf');
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
