%% shelf slopes
folders = { ...
    'runew-8041', 'runew-8042','runew-8150', 'runew-8151', ...
    ... %'runew-82342', 'runew-82343', 'runew-82344' ...
    'runew-8341', 'runew-8352', ...
    'runew-8383', 'ew-8384', 'ew-8385' ...
          };
sh = runArray(folders);

%% volume budget
handles = sh.array(4).VolumeBudget;
xlim([0 450]);
handles.hax.XLabel.Position(1) = 400;
handles.hisponge.DisplayName = 'Along-shelf: supply';

export_fig -r150 -a2 images/paper3/ew-8341-volume-budget.png

%% compare flat and sloping
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

export_fig -opengl -r150 -a2 images/paper3/sbsnapshot.png
