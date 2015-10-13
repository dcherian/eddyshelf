%% surface dye field - anticyclone and cyclone
%ew34 = runArray({'runew-34', 'runew-34-cyc'});
figure;
ax(1) = subplot(211);
ew34.array(1).animate_field('dye_03', ax(1), 320, 1);
title('Anticyclone moving southwards');
ax(2) = subplot(212);
ew34.array(2).animate_field('dye_03', ax(2), 290, 1);
title('Cyclone moving northwards');
export_fig images/ew-34-cyc-acyc-eddye.png

%% eddye y-z section
%ew36 = runs('../topoeddy/runew-36/');
tind = 320; limx = [0 150];
ew36.plot_yzsection(tind)
ylim([-400 0]);
xlim(limx);
hax = gca;
hax.YTick = sort([hax.YTick -floor(ew36.bathy.hsb)]);
suplabel('', 't'); title('Eddy dye');
text(0.85*limx(2), -1 * ew36.eddy.Lgauss(tind), ...
         {'vertical','scale'}, 'VerticalAlignment', 'Bottom', ...
         'HorizontalAlignment','Center');
export_fig images/ew-36-eddye-yz-320.png