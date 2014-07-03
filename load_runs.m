%{
% inviscid - vary Ro
ew03 = runs('../topoeddy/runew-03-nobg-2');
ew04 = runs('../topoeddy/runew-04-nobg-2');
ew05 = runs('../topoeddy/runew-05-nobg-2');
ew06 = runs('../topoeddy/runew-06-nobg');

% deeper shelfbreak
ew23 = runs('../topoeddy/runew-23-nobg/');
ew24 = runs('../topoeddy/runew-24-nobg/');
%ew25 = runs('../topoeddy/runew-25-nobg/');
%ew23.fluxes; ew24.fluxes; ew25.fluxes;
%ew23.water_census; ew24.water_census; ew25.water_census;
%ew23.jetdetect; ew24.jetdetect; ew25.jetdetect;

% with bottom friction
ew540 = runs('../topoeddy/runew-540');
ew541 = runs('../topoeddy/runew-541');
ew542 = runs('../topoeddy/runew-542');
%}

%% do this as an array
rootdir = '../topoeddy/';
folders = { ...
    'runew-03-nobg-2/', ...
    'runew-04-nobg-2/', ...
    'runew-05-nobg-2/', ...
    'runew-06-nobg/', ...
    'runew-23-nobg/', ...
    'runew-24-nobg/', ...
    'runew-630-nobg/', ...
    'runew-631-nobg/', ...
%'runew-632-nobg/', ...
    'runew-640-nobg/', ...
    'runew-641-nobg/', ...
    'runew-650-nobg/', ...
    'runew-651-nobg/', ...
    'runew-543-nobstress/', ...
    ...%   'runew-540/', ...
    ...%   'runew-541/', ...
    ...%   'runew-542/', ...
    ...%   'runew-544-nobg/', ...
          };

%% VERTICAL SCALES
folders = { ...
    '../topoeddy/runew-04-nobg-2/', ...
    '../topoeddy/runew-741-nobg/', ...
    '../topoeddy/runew-742-nobg/', ...
    '../topoeddy/runew-743-nobg/', ...
    '../topoeddy/runew-745-nobg/', ...
          };
vscales = runArray(folders);

for ii=1:vscales.len
    run = vscales.array(ii);
    tind = find_approx(run.eddy.t/run.eddy.tscale*86400, 1);
    vscales.name{ii} = num2str(run.eddy.Lgauss(3));
    %vscales.name{ii} = ...
    %    num2str(run.eddy.V(tind)/run.eddy.Lgauss(tind)/sqrt(1e-5))
end
vscales.plot_fluxes;

%% BOTTOM FRICTION
folders = { ...
    '../topoeddy/runew-04-nobg-2/', ...
    '../topoeddy/runew-540/', ...
    '../topoeddy/runew-541-nobstress/', ...
    '../topoeddy/runew-542/', ...
    '../topoeddy/runew-543-nobg/', ...
    '../topoeddy/runew-544-nobg/', ...
    '../topoeddy/runew-546-nobstress/', ...
          };
bfrics = runArray(folders);
for ii=1:bfrics.len
    bfrics.name{ii} = num2str(bfrics.array(ii).params.misc.rdrg);
end

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
%   'runew-632-nobg/', ...
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