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

%% flux diagnostics
handles = ew34.plot_fluxts(factor, 3, 3);
handles.htitle.String = ['Flux of water above ' num2str(factor) 'H_{sb}'];
maximize; drawnow;
% mark timesteps
export_fig -a1 images/paper2/flux-diags.png

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

%% secondary eddy
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

%% z-slice
opt.csdcontours = ew34.csflux.x(isobath);
handles = ew34.mosaic_zslice('dye_03', 100, [225 235 245 255], opt);
for ii=1:length(handles)
    handles(ii).eddsurf.Visible = 'off';
    handles(1).htext(2).Visible = 'off';
    handles(ii).rhocont.LineWidth = 2.5;
    handles(ii).csdsurf.LineWidth = 2.5;
end
xlim([200 400]);
ylim([0 150]);
handles(3).hax.XTickLabel{end} = '';
handles(1).htitle.String = 'Eddy dye at z = -200 m | H_{sb} = 50m | Ro = 0.10';
handles(1).htitle.FontSize = 22;
handles(1).htitle.FontWeight = 'normal';
handles(1).supax.Position(4) = 0.87;
axes(handles(1).hax); correct_ticks('y', '', '0');
axes(handles(3).hax); correct_ticks('y', '', '0');
export_fig -a1 images/paper2/ew-34-mosaic-zslice.png