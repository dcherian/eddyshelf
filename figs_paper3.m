%% shelf slopes
folders = { ...
    'ew-8041', 'ew-8042', ...
    ... %'ew-8150', 'ew-8151', ...
    'ew-82342', 'ew-82343', ... %'ew-82344' ...
    'ew-8341', 'ew-8361', ...
    'ew-8352', 'ew-8352-2', ...
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
sh.filter = [1:12];
sh.plot_fluxes(1,1,[],[0 1]);
hleg = findobj('type','legend');
hleg.delete;

export_fig -r150 images/paper3/sb-vert-profiles.png

%% eddy water on shelf scale
figure;
hax(1) = subplot(121);
hax(2) = subplot(122);

sh.print_diag('supply', [], hax(1), 'no_name_points');
sh.print_diag('eddyonshelf', [], hax(2), 'no_name_points');

hax(1).Title.String = '';
hax(2).Title.String = '';

export_fig -r150 images/paper3/parameterizations.png
%%
for ii=1:sh.len
    % sbssh = sh.array(ii).sbssh;
    % supply = sh.array(ii).supply;

    % figure; maximize;
    % plot(sbssh.xvec/1000, sbssh.ssh - min(sbssh.ssh)); hold on;
    % plot(supply.ymat(:,1)/1000 - sh.array(ii).bathy.xsb/1000, ...
    %      supply.zeta.zetamean - min(supply.zeta.zetamean));
    % title(sh.array(ii).name);
    %sh.array(ii).avgSupplyJet;
    %sh.array(ii).plotAvgSupplyJet;
    %export_fig('-r96',  ['images/supply-' sh.array(ii).name '.png']);

    % sh.array(ii).VolumeBudget;
    % export_fig(['images/volume-budget-' sh.array(ii).name '.png']);
end

%% volume budget
handles = sh.array(4).VolumeBudget;
xlim([0 450]);
handles.hax.XLabel.Position(1) = 400;
handles.hisponge.DisplayName = 'Along-shelf: supply';

export_fig -r150 -a2 images/paper3/ew-8341-volume-budget.png

%% compare flat and sloping
if ~exist('ew04', 'var')
    ew04 = runArray({'ew-04', 'ew-8041'});
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

export_fig -opengl -r150 -a2 images/paper3/sbsnapshot.png

%%
