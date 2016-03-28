% do this as an array
rootdir = '';
folders = { ...
    'runew-03', 'runew-04', 'runew-05', 'runew-06', ...
    'runew-13', 'runew-15', ...
    'runew-4040', 'runew-4050', ...
    'runew-5040', 'runew-5041', 'runew-5043', ...
    'runew-8040', 'runew-8041','runew-8042', ...
    'runew-8150', 'runew-8151', ...
    ... %'runew-2040-Bu05', ...
    ... %'runew-2340', 'runew-2341_wider', 'runew-2342_wider', ...
    ... %'runew-2360', 'runew-2360-20km', 'runew-2361_wider', ...
    ... %'runew-2362_wider', 'runew-2363_wider', ...
    ... %'runew-33', 'runew-34', 'runew-35', 'runew-36','runew-37', ...
    ... %'runew-36-sym', ...
    ... %'runew-3340', 'runew-3341-2', ...
    ... %'runew-4341', 'runew-4343', ...
    ... %'runew-5341', 'runew-5343', ...
          };
all = runArray(folders);

for ii=1:all.len
    if isempty(all.array(ii).csflux) || ...
        ~isfield(all.array(ii).csflux, 'westmask')
        all.array(ii).csfluxes;
    end
    %    if isempty(all.array(ii).asflux)
    %    all.array(ii).asfluxes;
    %end
end

%% Ro
%folders = { ...
%     'runew-03', 'runew-04/', 'runew-05', 'runew-06'};
%    'runew-13', 'runew-15/', ...
%    'runew-33', 'runew-34', 'runew-35', ...
%    'runew-4040', 'runew-4050', ...
%    'runew-8041', 'runew-8051', ...
%          };

folders = {'runew-03', 'runew-04', 'runew-05', 'runew-06'};
Rh3 = runArray(folders);

for ii=1:Ro3.len
    tind = Ro3.array(ii).eddy.tscaleind;
    Ro3.name{ii} = ['Ro = ' num2str(Ro3.array(ii).eddy.Ro(tind), '%.2f') ' | Rh = ' ...
                    num2str(Ro3.array(ii).params.nondim.eddy.Rh, '%.2f')];
end

folders = {'runew-33', 'runew-34', 'runew-35', 'runew-36', 'runew-37'};
Rh12 = runArray(folders);

for ii=1:Ro12.len
    tind = Ro12.array(ii).eddy.tscaleind;
    Ro12.name{ii} = ['Ro = ' num2str(Ro12.array(ii).eddy.Ro(tind), '%.2f') ' | Rh = ' ...
                     num2str(Ro12.array(ii).params.nondim.eddy.Rh, '%.2f')];
end

folders = { ...
    'runew-04', 'runew-34', ...
    'runew-05', 'runew-35', ...
    'runew-06' ,'runew-36', ...
          };
Ro = runArray(folders);

ew09 = runArray({'runew-09-dx750', 'runew-09-teles4-2'});

% ew
folders = { ...
    'runew-34', 'runew-35', 'runew-36', ...
    'runew-2340', 'runew-2341_wider', ...
    'runew-2361_wider', ...
    'runew-3340', 'runew-23340', ...
          };
ew = runArray(folders);

%% gaussian pressure
folders = { ...
    'runew-a330', 'runew-a340', 'runew-a350'
          }'
ewa = runArray(folders);

%% 642
folders = { ...
    'runew-6042-new', ...
    'runew-642-fb', ...
          };
ew642 = runArray(folders);

%% process
folders = { ...
'runew-6341', 'runew-6341-2', 'runew-6341-3', 'runew-6341-4', ...
          };
process = runArray(folders);;

for ii=1:process.len
    process.array(ii).csfluxes;
end

%%%
ew = runArray({'runew-64361', 'runew-64361-shallow', 'runew-64361-deep'});

%% asfluxes
folders = { ...
    'runew-04', 'runew-34' ...
          };
beta = runArray(folders);

%% NS isobaths
folders = { ...
    ... %'runns-64361', 'runns-6341', 'runns-6362-2',...
    ... %'runns-6441', ...
    'runns-35', 'runnw-36', ...
    'runns-2360', 'runns-2360_wider', ...
          };
ns = runArray(folders);

for ii=3:ns.len
    if ~isfield(ns.array(ii).eddy, 'rhovor')
        track_eddy(ns.array(ii));
    end
    ns.array(ii).csfluxes;
end

% look at deep shelfbreaks
folders = { ...
    'runew-34', 'runew-35', 'runew-36', 'runew-37', ...
    'runew-6370', ... %'runew-6341', 'runew-6371', ...
    'runew-2041', 'runew-2043', ...
    'runew-2340', 'runew-2345', ...
    'runew-2341_wider','runew-2361_wider', 'runew-2363_wider',  ...
    'runew-2360_wider', 'runew-2365-75km', ...
          };
sb = runArray(folders);

for ii=1:sb.len
    sb.array(ii).avgStreamerVelSection;
end

%% shelf slopes
folders = { ...
    'runew-8041', 'runew-8042',' runew-8150', 'runew-8151', ...
    ... %'runew-82342', 'runew-82343', 'runew-82344' ...
    'runew-8341', 'runew-8352', ...
    'runew-8383', 'ew-8384', 'ew-8385' ...
          };
sh = runArray(folders);

for ii=1:sh.len
    sh.array(ii).avgSupplyJet;
end

%% csflux
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
        %'runew-8040', ...
    'runew-8041', 'runew-8042',' runew-8150', 'runew-8151', ...
    ... %'runew-82342', 'runew-82343', 'runew-82344' ...
    'runew-8341', 'runew-8352', ...
    'runew-583411', 'runew-583413', ...
    'runew-8383', ...
    ... %'runew-b4361', ...
          };
csf = runArray(folders);

comm6filter = [4:20 30];
didntwork = [];
for ii=1:csf.len
    % if isempty(strfind(csf.array(ii).eddy.hash, '31d29d2416ec1199'))
    %     track_eddy(csf.array(ii));
    %     csf.array(ii).usurf = []; csf.array(ii).vsurf = [];
    %     csf.array(ii).zeta = []; csf.array(ii).rhosurf = [];
    %     csf.array(ii).vorsurf = [];
    % end
    try
        csf.array(ii).csfluxes;
        %csf.array(ii).avgStreamerVelSection;
    catch ME
        disp(ME.message);
        didntwork = [didntwork ii];
    end
end

%% shfric
folders = { ...
    'runew-34', 'runew-8341', ...
    'runew-583411', 'runew-583413', ...
    'runew-583414', 'runew-583415', ...
    };
shfric = runArray(folders);

%% process
process = runArray({'runew-4341', 'runew-4344', 'runew-8340', 'runew-583413'});
didntwork = [];
for ii=1:process.len
    try
        process.array(ii).csfluxes;
        process.array(ii).avgStreamerVelSection;
    catch ME
        disp(ME.message)
        didntwork = [didntwork ii];
    end
end
didntwork

%% slope
folders = { ...
    'runew-6341', 'runew-6341-2', 'runew-6341-3', 'runew-6341-4', ...
    'runew-56341', 'runew-56341-2', ...
    'runew-6342', 'runew-6352', ...
    'runew-6362', 'runew-6362-1', 'runew-6362-2', ...
    'runew-6371', 'runew-6372', ...% 'runew-6373' ...
    'runew-6441', ...
    'runew-64461-3', 'runew-64462', ...
    'runew-64461-4', 'runew-64461-5', ...
    'runew-64461-6', 'runew-64461-7', ...
    'runew-64463-1', ...
    'runew-b4361', ...
    'runew-64351-cyc', 'runew-6341-fneg', ...
    ... % 'runew-64361', moved to external
    'runew-64361-2', ...
    'runew-564361', 'runew-564361-2', 'runew-564361-3', ...
    'runew-64361-deep', 'runew-64361-shallow', 'runew-64361-sl' ...
    'runew-64461-3', 'runew-5644613-1/', ...
    'runew-64461-9-shallow', 'runew-64461-9-deep', ...
    ... %'runew-05', 'runew-6052', ...
    ... %'runew-06', 'runew-6062', ...
          };
sl = runArray(folders);

for ii=1:sl.len
    tind = sl.array(ii).eddy.tscaleind;
    run = sl.array(ii);

    slSa(ii) = run.bathy.S_sl;
    slRo(ii) = run.eddy.Ro(1);
    slRh(ii) = run.params.nondim.eddy.Rh;
    beta(ii) = run.params.phys.beta;
    beta_t(ii) = run.params.phys.f0 ./ run.eddy.Lgauss(1) ...
        .* run.bathy.sl_slope;
    blf(ii) = beta(ii) .* run.eddy.vor.dia(1)/2 ...
              ./ run.eddy.fcen(1);
    %slRo = run.eddy.V(1)./run.params.phys.f0/run.bathy.L_slope;
    %sl.name{ii} = ['S_\alpha = ' num2str(slSa(ii), '%.2f') ' | Ro = ' ...
    %               num2str(slRo(ii), '%.2f') ' | Ro/S_\alpha = ' ...
    %               num2str(slRo(ii)./slSa(ii), '%.2f')];
    %sl.name{ii} = ['Ro/S_\alpha = ' ...
    %               num2str(slRo(ii)./slSa(ii), '%.2f') ...
    %               ' | ' run.name];
    %sl.name{ii} = ['Rh = ' ...
    %               num2str(slRh(ii), '%.2f') ...
    %               ' | ' run.name];
    sl.name{ii} = [run.name];
end
sl.sort(blf);

for ii=1:3
    sl.array(ii).asfluxes;
end

%%
folders = { ...
    ... %    'runew-332/', ...
    'runew-04', ...
    'runew-15', ...
    'runew-34/', ...
    'runew-35/', ...
    'runew-36/', ...
          };
beta = runArray(folders);

folders = {'runew-04', 'runew-34'};
temp = runArray(folders);

%% shelfbreak depth
folders = {...
    'runew-34', 'runew-2340', ... %'runew-2341_wider', 'runew-2342_wider', ...
    'runew-2345', ...
    'runew-36', ... %'runew-2360', ...
    'runew-2360_wider', 'runew-2365-75km', ...
    ... %'runew-2360-20km', 'runew-2361_wider', 'runew-2362_wider', ...
    ... %'runew-2363_wider', ...
    'runew-3340-80km', 'runew-23340', ...
    'runew-35', 'runew-37', ... %'runew-9340', 'runew-9341', ...
    };
sb = runArray(folders);

for ii=1:sb.len
    tind = sb.array(ii).eddy.tscaleind;
    lambda(ii) = sb.array(ii).bathy.hsb ./ ...
                           sb.array(ii).eddy.Lgauss(1);
    sb.name{ii} = ['\lambda = ' ...
                   num2str(lambda(ii), '%.2f') ...
                           ' | ' sb.array(ii).name];
end
sb.sort(lambda);

%%ns
folders = { ...
    'runns-35/', 'runns-36/' ...
    'runns-2360_wider', ...
          };
ns = runArray(folders);

%%ew-4
folders = { ...
    'runew-34/', ...
    'runew-3340/', ...
%    'runew-4040/', ...
%    'runew-5040', 'runew-5041', 'runew-5043', ...
%    'runew-6040/', 'runew-6041/', 'runew-6042', ...
%    'runew-34/', ...
%    'runew-6343', ...
%    'runew-8040/', 'runew-8041/', 'runew-8042/', ...
          };
ew4 = runArray(folders);
for ii=1:ew4.len
    ew4.name{ii} = [ew4.array(ii).name ' | ' ...
                    num2str(ew4.array(ii).eddy.Ro(ew4.array(ii).tscaleind))];
end

%%ew-5
folders = { ...
    'runew-05/', ...
    'runew-15/', ...
    'runew-4050/', ...
    'runew-6150-closer/', ...
    'runew-6151/', 'runew-6152/', ...
    'runew-35/', ...
    'runew-8150/', 'runew-8151/', ...
          };
ew5 = runArray(folders);
for ii=1:ew5.len
    ew5.name{ii} = [ew5.array(ii).name ' | ' ...
                    num2str(ew5.array(ii).eddy.Ro(1)) ' | ' ...
                    num2str(ew5.array(ii).params.nondim.eddy.Rh)];
end

%% VERTICAL SCALES - FLAT BOTTOM
folders = { ...
    'runew-741-flat/', ...
    'runew-742-flat/', ...
    'runew-743-flat/', ...
    'runew-745-fb/', ...
    'runew-746-flat/', ...
          };
vscalesflat = runArray(folders);

for ii=2:vscalesflat.len
    track_eddy(vscalesflat.array(ii));
end

%% 745
folders = { ...
    'runew-745-nobg/', ...
    'runew-745-fb/', ...
    'runew-745-bg/', ...
    'runew-745-flat-bg/', ...
          };
ew745 = runArray(folders);
ew745.array(2).eddy.tscale = ew745.array(1).eddy.tscale;
ew745.array(4).eddy.tscale = ew745.array(1).eddy.tscale;

%% 742
folders = { ...
    'runew-742-nobg/', ...
    'runew-742-flat/', ...
    'runew-742-a22/', ...
    'runew-742-nobg-a22/', ...
          };
ew742 = runArray(folders);
ew742.array(2).eddy.tscale = ew742.array(1).eddy.tscale;
ew742.array(3).eddy.tscale = ew742.array(1).eddy.tscale;

%% VERTICAL SCALES
folders = { ...
    'runew-04-nobg-2/', ...
    'runew-741-nobg/', ...
    'runew-742-nobg/', ...
    'runew-743-nobg/', ...
    'runew-745-nobg/', ...
          };
vscales = runArray(folders);

for ii=1:vscales.len
    run = vscales.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    %vscales.name{ii} = num2str(run.eddy.Lgauss(3));
    vscales.name{ii} = ...
        [num2str(run.eddy.Lgauss(3)) ' m | ' ...
         num2str(run.eddy.V(tind)/run.eddy.Lgauss(tind)/sqrt(1e-5))];
end

for ii=1:vscales.len

end

vscales.plot_fluxes;

%% BOTTOM FRICTION
folders = { ...
    ... %'runew-34', 'runew-5341', 'runew-5343'
    ... %'runew-6341', 'runew-56341', 'runew-56341-2', ...
    'runew-64361', ...
    'runew-564361', 'runew-564361-2', 'runew-564361-3', ...
    'runew-64361-sl', 'runew-64361-sl-fb', ...
    'runew-64361-fb', 'runew-564361-fb', ...
    'runew-64461-3', 'runew-5644613-1/', ...
          };
bfrics = runArray(folders);
bfrics.filter = [];
diags = bfrics.print_diag('bfric');
for ii=1:bfrics.len
    run = bfrics.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    %bfrics.name{ii} = num2str(bfrics.array(ii).params.misc.rdrg);
    bfrics.name{ii} = num2str(diags(ii));
end

%%
tind = [];
figure;
hax = subplot(122);
bfrics.array(2).animate_field('pbot', hax, tind, 1); clim = caxis;
hax = subplot(121);
bfrics.array(1).animate_field('pbot', hax, tind, 1); caxis(clim);

%%
var = 'rho'; yy = 160e3;
bfrics.filter = 1;
bfrics.plot_sections(var, find_approx(bfrics.fitler);
title(bfrics.array(bfrics.filter).name);
ax(1) = gca;
bfrics.filter = 2; bfrics.plot_sections(var, 142)
title(bfrics.array(bfrics.filter).name);
ax(2) = gca;
linkaxes(ax, 'xy')

%% MISC
for ii=2:vscales.len
    try
        roms_pv(vscales.array(ii).dir, [], 'his');
    catch ME
        disp(ME);
        disp(vscales.array(ii).name)
    end
end
%folders = { ...
%    'runew-630-nobg/', ...
%    'runew-631-nobg/', ...
%    'runew-632-nobg/', ...
%    'runew-640-nobg/', ...
%    'runew-641-nobg/', ...
%    'runew-642-nobg/', ...
%    'runew-650-nobg/', ...
%    'runew-651-nobg/', ...
%    };
kk = 1;
for ii = 1:length(folders)
    warning off;
    try
        array(kk) = runs([rootdir folders{ii}]);
        array(kk).name
        kk = kk + 1
    catch ME
        disp(ME)
        continue;
    end
end

%% flat v/s topo
ew05 = runs('../topoeddy/runew-05-nobg-2/');
ew05flat = runs('../topoeddy/runew-05-flat/');

topo = ew05;
flat = ew05flat;

topo = ew744;
flat = ew744flat;

figure; hold all
plot(flat.eddy.t/topo.eddy.tscale * 86400, flat.eddy.Lgauss);
plot(topo.eddy.t/topo.eddy.tscale * 86400, topo.eddy.Lgauss);
legend('flat', 'topo');


figure(7);
hold all
for ii=1:length(array)
    ndtime = array(ii).eddy.t/array(ii).eddy.tscale * 86400;
    hgplt = plot(ndtime, array(ii).eddy.Lgauss);
    addlegend(hgplt, array(ii).name);
end


for ii=1:length(array)
    disp(array(ii).params.eddy.dia);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% COMMITTEE MEETING II

%% ew-04 runs
fontSizes = [];
run03 = [1 5 7 8];
run04 = [2 6 9 10 13];
colors = distinguishable_colors(15);

indices = [1 2 3 4];
names = { ...
    'Ro = 0.06', ...
    'Ro = 0.10', ...
    'Ro = 0.20', ...
    'Ro = 0.40'};

indices = run04;
names = { ...
    'Base case', ...
    'H_{sb} = 75m', ...
    'S = 1.25', ...
    'S = 0.96', ...
    'linear drag = 5e-4 m/s'};


indices=  [1 2]
names = { ...
    'rdrg = 5e-4 m/s', ...
    'rdrg = 5e-3 m/s', ...
        };
figure;
subplot(2,1,1); hold all
subplot(2,1,2); hold all
for jj=1:length(indices)
    %if jj == 1
    %    alphaval = 1
    %else
    %    alphaval = 0.7
    %end
    ii = indices(jj);
    ndtime = array(ii).csflux.time(1:end-2) / array(ii).eddy.tscale;
    subplot(2,1,1)
    hgplt = plot(ndtime, array(ii).csflux.west.shelf(1:end-2)/1e6, ...
                 'Color', colors(jj,:));
    addlegend(hgplt, names{jj}, 'NorthWest');
    plot(ndtime, array(ii).csflux.east.eddy(1:end-2)/1e6, 'LineStyle', ...
         '--', 'Color', colors(jj,:));

    subplot(2,1,2)
    plot((array(ii).csflux.west.shelfwater.bins/1000 - ...
          array(ii).bathy.xsb/1000) / (array(ii).rrshelf/1000), ...
         array(ii).csflux.west.shelfwater.itrans, 'Color', colors(jj,:));
end
subplot(2,1,1)
liney(0, [], [1 1 1]*0.7);
xlabel('Non-dimensional time');
ylabel('Transport (Sv)');
beautify(fontSizes);
subplot(2,1,2)
xlabel('Distance from shelfbreak / Shelfbreak Rossby Radius');
xlim([-30 0]);
ylabel('Volume (m^3)');
beautify(fontSizes);

%% variation of slope burger number

figure;
subplot(2,1,1)
hold all
subplot(2,1,2)
hold all

indices = run04;
names = { ...
    'Base case', ...
    'S = 1.25', ...
    'S = 0.96', ...
        };

for jj=1:length(indices)
    %if jj == 1
    %    alphaval = 1
    %else
    %    alphaval = 0.7
    %end
    ii = indices(jj);
    run = array(ii);
    ndtime = run.eddy.t ./ run.eddy.tscale * 86400;
    subplot(2,1,1)
    hgplt = plot(ndtime, array(ii).eddy.hcen);
    addlegend(hgplt, names{jj}, 'NorthWest');

    subplot(2,1,2)
    plot(ndtime, array(ii).eddy.V./array(ii).eddy.Ls/array(ii).params.phys.f0);
    %    plot(array(ii).csflux.west.shelfwater.bins/1000 - ...
    %     array(ii).bathy.xsb/1000, ...
    %     array(ii).csflux.west.shelfwater.itrans);
end

%% vertical scales
indices = run04;
names = { ...
    'Base case', ...
    'H_{sb} = 75m', ...
    'S = 1.25', ...
    'S = 0.96', ...
    'linear drag = 5e-4'};

indices = [1 2 3 4];
names = { ...
    'Ro = 0.06', ...
    'Ro = 0.10', ...
    'Ro = 0.20', ...
    'Ro = 0.40'};

figure;
subplot(1,2,1); hold all
subplot(1,2,2); hold all
%subplot(2,1,3); hold all
%for jj=1:length(indices)
%    ii = indices(jj);
for ii=1:length(array)
    jj = ii
    run = array(ii);
    ndtime = run.eddy.t / run.eddy.tscale * 86400;

    subplot(1,2,1)
    hgplt = plot(ndtime, run.eddy.Lgauss, 'color', colors(jj,:));
    addlegend(hgplt, run.name);
    plot(ndtime, run.eddy.hcen, 'color', colors(jj,:), ...
         'LineStyle', '--');

    subplot(1,2,2)
    hold all
    plot(ndtime, run.eddy.prox/1000,'Color', colors(jj,:));
    plot(ndtime, run.eddy.my/1000 - run.bathy.xsb/1000, 'Color', ...
         colors(jj,:), 'LineStyle', '--');
end
subplot(1,2,1)
ylabel('(m)');
title(['Dashed = water depth at eddy center | Solid = vertical scale ' ...
       'of eddy']);
beautify(fontSizes);
subplot(1,2,2)
ylabel('(km) from shelfbreak');
xlabel('non-dim time');
liney(0, [], [1 1 1]*0.5);
beautify(fontSizes);
title('Dashed = center | solid = southern edge');

%%

figure;
hold all
plot(ew744.csflux.time/ew744.eddy.tscale, ew744.csflux.west.shelf);
plot(ew745nobg.csflux.time/ew745nobg.eddy.tscale, ew745nobg.csflux.west.shelf);
legend('744','745');

%%

ew543.plot_shelfvorbudget;
subplot(2,1,1)
hold all;
plot(ew04.vorbudget.time/86400, ew04.vorbudget.shelf.rv, 'Color', [0.68, ...
                    0.85, 0.90]);
legend('linear drag = 5e-4 m/s', '', 'inviscid', 'Location', 'NorthWest');

%%
ew543.animate_vorbudget(120,0);
ew04.animate_vorbudget(120,0);
array(9).animate_vorbudget(120,0);
array(10).animate_vorbudget(120,0);
ew24.animate_vorbudget(120,0);