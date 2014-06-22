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

%folders = { ...
%    '../topoeddy/runew-04-nobg-2/', ...
%    '~/scratch/topoeddy/runew-741-nobg/', ...
%    '~/scratch/topoeddy/runew-742-nobg/', ...
%    '~/scratch/topoeddy/runew-743-nobg/', ...
%          };

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
        continue;
    end
end

%% flat v/s topo
figure; hold all
plot(ew05flat.eddy.t/ew05flat.eddy.tscale * 86400, ew05flat.eddy.Lgauss);
plot(ew05.eddy.t/ew05flat.eddy.tscale * 86400, ew05.eddy.Lgauss);
legend('flat', 'topo');

%%
figure(1)
clf; hold all
figure(2)
clf; hold all
for ii=1:length(array)
    %    hgplt = plot(array(ii).eddy.t/array(ii).eddy.tscale*86400, ...
    %             array(ii).eddy.prox);    
    eddy_ndtime = array(ii).eddy.t/array(ii).eddy.tscale*86400;
    csflx_ndtime = array(ii).csflux.time/array(ii).eddy.tscale * 86400;
    etind = find_approx(eddy_ndtime, 1.5, 1)
    cstind = find_approx(csflx_ndtime, 1.5, 1);

    meanprox(ii) = nanmean(array(ii).eddy.hcen(etind:end));
    meanflux(ii) = nanmean(array(ii).csflux.west.shelf(cstind:end));

    param(ii) = array(ii).params.nondim.eddy.Ro/ ...
        array(ii).params.nondim.S_sl .* ...
        array(ii).params.nondim.eddy.Bu;

    figure(1);
    hgplt = plot(param(ii), meanprox(ii), '.');
    addlegend(hgplt, array(ii).name);
    disp(['run = ', array(ii).name , ' | mean prox. = ', ...
          num2str(meanprox(ii))]);
    %    pause;

    figure(2);
    hgplt = plot(param(ii), meanflux(ii), '.');
    addlegend(hgplt, array(ii).name);
    disp(['run = ', array(ii).name , ' | mean prox. = ', ...
          num2str(meanflux(ii))]); 
end

g
% look at frictional runs
figure(5);
hold all
figure(6);
hold all;
for ii=1:length(array)
    if ii == 2 || ...
           array(ii).params.misc.rdrg ~= 0

        ndtime = array(ii).eddy.t/array(ii).eddy.tscale * 86400;
        figure(5);
        subplot(2,1,1)
        hold all
        hgplt = plot(ndtime, array(ii).eddy.V);
        addlegend(hgplt, array(ii).name);

        subplot(2,1,2)
        hold all
        plot(ndtime, array(ii).eddy.Ls);

        figure(6);
        hgplt = plot(array(ii).csflux.time/array(ii).eddy.tscale * ...
                     86400, smooth(array(ii).csflux.west.shelf/1e6, 3));
        addlegend(hgplt, array(ii).name);
    end
end
subplot(2,1,1)
title('Eddy velocity scale (m/s)');
subplot(2,1,2)
title('Eddy Length scale');
xlabel('Non-dim. time');
figure(6);
title('Fluxes');

figure(7);
hold all
for ii=1:length(array)
    ndtime = array(ii).eddy.t/array(ii).eddy.tscale * 86400;
    hgplt = plot(ndtime, array(ii).eddy.Lgauss);
    addlegend(hgplt, array(ii).name);
end

figure;
hold all
for ii=1:length(array)
    ndtime = array(ii).eddy.t; %/array(ii).eddy.tscale * 86400;
    hgplt = plot(ndtime, array(ii).eddy.hcen);
    addlegend(hgplt, array(ii).name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% COMMITTEE MEETING II

%% ew-04 runs
fontSizes = [];
figure;
subplot(2,1,1); hold all
subplot(2,1,2); hold all
run03 = [1 5 7 8];
run04 = [2 6 9 10 13];

indices = run04;
names = { ...
    'Base case', ...
    'H_{sb} = 75m', ...
    'S = 1.25', ...
    'S = 0.96', ...
    'linear drag = 5e-4'};

for jj=1:length(indices)
    %if jj == 1
    %    alphaval = 1
    %else
    %    alphaval = 0.7
    %end
    ii = indices(jj);
    ndtime = array(ii).csflux.time(1:end-2) / array(ii).eddy.tscale;
    subplot(2,1,1)
    hgplt = plot(ndtime, array(ii).csflux.west.shelf(1:end-2)/1e6);
    addlegend(hgplt, names{jj}, 'NorthWest');
    
    subplot(2,1,2)
    plot(array(ii).csflux.west.shelfwater.bins/1000 - ...
         array(ii).bathy.xsb/1000, ...
         array(ii).csflux.west.shelfwater.itrans);
end
subplot(2,1,1)
liney(0, [], [1 1 1]*0.7);
title(['Preserve Ro = 0.05, \alpha = 0.02,  Eddy diameter/L_{slope} = 2 | Base case = ' ...
       'inviscid, H_{sb} = 50m, S = 1.5'], 'interpreter', 'tex');
xlabel('Non-dimensional time');
ylabel('Transport (Sv)');
beautify(fontSizes);
subplot(2,1,2)
xlabel('Distance from shelfbreak');
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

%% eddy water on shelf
indices = run03;
names = { ...
    'Base case', ...
    'H_{sb} = 75m', ...
    'S = 1.25', ...
    'S = 0.96', ...
    'linear drag = 5e-4'};

figure;
hold all

for jj=1:length(indices)
    ii = indices(jj);
    run = array(ii);
    ndtime = run.csflux.time ./ run.eddy.tscale;
    hgplt = plot(ndtime, run.csflux.east.eddy/1e6);
end
ylabel('Transport (Sv)');