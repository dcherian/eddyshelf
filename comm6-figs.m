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