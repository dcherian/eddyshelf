fname = '../data/seep1.mat';

seep = load(fname);
datevec = datetime(seep.mday, 'ConvertFrom', 'datenum');

% velocity is stored in complex form

%%
figure
plot(datevec, seep.wtsp6(:,1));
hold on
plot(datevec, seep.wtsp4(:,2));
plot(datevec, seep.wtsp7(:,1));
ylim([10 25]);
xlim([datenum('1983-09-01') datenum('1983-10-31')]);
linex([datenum('1983-10-16') datenum('1983-10-21')]);

%%
figure; insertAnnotation('PlotSEEP1.m'); maximize;
ax = plotyy(datevec, seep.wtsp7(:,1), ...
            datevec, seep.saltsp7(:,1));
ax(1).YLabel.String = 'Temp';
ax(2).YLabel.String = 'Salinity';
linex([datenum('1983-09-25') datenum('1983-10-07')],  ...
       {'Cyclone'; 'Anticyclone'});
xlim([datenum('1983-09-01') datenum('1983-10-31')])
linkprop(ax, {'XLim','XTick'});
ax(1).XTickMode = 'auto';
axes(ax(1)); beautify;
axes(ax(2)); beautify;
title('Mooring 7 at 400m');