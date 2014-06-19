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
    'runew-540/', ...
    'runew-541/', ...
    'runew-542/', ...
    'runew-543-nobg/', ...
    'runew-544-nobg/', ...
    'runew-630-nobg/', ...
    'runew-631-nobg/', ...
    'runew-640-nobg/', ...
    'runew-641-nobg/', ...
          };

kk = 1;
for ii = 7:length(folders)
    warning off;
    try
        array(kk) = runs([rootdir folders{ii}]);
        array(kk).name
        kk = kk + 1
    catch ME
        continue;
    end
end

%%


for ii=1:4%length(array)
    figure;
    plot(array(ii).csflux.time/86400, smooth(array(ii).csflux.west.shelf, ...
                                             2));
    hold all
    plot(array(ii).csflux.time(1:2:end)/86400, ...
         smooth(array(ii).csflux.west.shelf(1:2:end), 1));
    title(array(ii).name);
end

%% do plots
figure(1);
clf; hold all
%figure(2);
%clf; hold all;
%figure(3);
%clf; hold all;
figure(4);
clf; hold all;
for ii = 1:length(array)
    ndtime = array(ii).csflux.time/array(ii).eddy.tscale;
    figure(1);
    subplot(2,1,1)
    hold all
    hgplt = plot(ndtime, smooth(array(ii).csflux.west.shelf/1e6, 3));
    addlegend(hgplt, array(ii).name);
    subplot(2,1,2)
    hold all
    plot(ndtime, cumtrapz(array(ii).csflux.time, ...
                          repnan(array(ii).csflux.west.shelf, 0)));

    %    figure(2);
    %hgplt = plot(array(ii).time/array(ii).eddy.tscale, ...
    %             array(ii).water.sh.slope + array(ii).water.sh.deep);
    %addlegend(hgplt, array(ii).name);

    %figure(3);
    %hgplt = plot(array(ii).time/array(ii).eddy.tscale, ...
    %             array(ii).water.sl.shelf);
    %addlegend(hgplt, array(ii).name);

    figure(4);
    subplot(2,1,1)
    hold all
    %    if array(ii).params.misc.rdrg ~= 0
    hgplt = plot(ndtime, ...
                 smooth(array(ii).csflux.west.shelfwater.envelope, ...
                        4));
    addlegend(hgplt, array(ii).name);
    %end
    subplot(2,1,2)
    hold all
    ndx = (array(ii).csflux.west.shelfwater.bins / 1000 - ...
           array(ii).bathy.xsb/1000);
    plot(ndx, array(ii).csflux.west.shelfwater.itrans);
end
figure(1);
subplot(2,1,1)
title('Flux of shelfwater (Sv)');
subplot(2,1,2)
title('Total volume transported (m^3)');
figure(2);
title('Volume of exported shelf water');
figure(3);
title('Volume of slope water on shelf');
figure(4);
subplot(2,1,1)
title('Farthest location from where shelfwater is exported');
subplot(2,1,2)
title('Transport binned by location on shelf');
xlabel('Location (km)');
ylabel('Volume (m^3)');

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
    hgplt = plot(ndtime, array(ii).eddy.prox);
    addlegend(hgplt, array(ii).name);
end