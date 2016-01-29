fname = '../data/All_Oleander_3D.mat';

ole = load(fname);

%% churchill (1986) section
% Oct 19-20, 1983
% stored in reverse order, that's why I reverse axis direction.
tt = 114;
zint = [-800:1:0];
clear Tint

figure;
insertAnnotation('PlotOleander.m');
ax(1) = subplot(1,3,[1 2]);
for aa=1:size(ole.Temp_Fix, 2) % along-shelf direction
    if isempty(cut_nan(squeeze(ole.Depth(tt,aa,:))))
        break;
    end
    Tint(aa,:) = interp1(-1*cut_nan(squeeze(ole.Depth(tt,aa,:))), ...
                         cut_nan(squeeze(ole.Temp_Fix(tt,aa,:))), ...
                         zint);
end
contourf(ole.Dist(tt,1:size(Tint,1))/1000, zint, Tint', 40);
set(gca, 'xdir', 'reverse');
linex(ole.Dist(tt,[11 17])/1000);
colorbar;
ylabel('Z (m)'); xlabel('Along-track Distance (km)');
title('XBT Data | Oleander | 19-20 Oct, 1986');
beautify;

ax(2) = subplot(1,3,3);
hh(1) = plot(squeeze(ole.Temp_Fix(tt,17,:)), -1*squeeze(ole.Depth(tt,17,:)), 'o', ...
             'LineWidth', 2);
hold on;
plot(squeeze(ole.Temp_Fix(tt,11,:)), -1*squeeze(ole.Depth(tt,11,:)), 'o', ...
     'LineWidth', 2);
hh(2) = plot(Tint(17,:), zint, '--');
plot(Tint(11,:), zint, '--');
ax(2).XAxisLocation = 'top';
xlabel('Temp (C)');
legend(hh, 'XBT data', 'Interpolated', 'Location', 'SouthEast');
beautify;

linkaxes(ax, 'y');
ylim([-500 0]);

export_fig -r150 -a2 images/oleander-oct1986.png

%contourf(max(ole.Dist(tt,:)) - repmat(ole.Dist(tt,:)', [1 500]), ...
%         -1*squeeze(ole.Depth(tt,:,:)), squeeze(ole.Temp_Fix(tt,:,:)), 40);
%colorbar;
%xlim([nanmin(ole.Dist(:)) nanmax(ole.Dist(:))/8]);

%%
tt = 480;
for tt=480:500
    cla
    contourf(repmat(ole.Dist(tt,:)', [1 500]), ...
             -1*squeeze(ole.Depth(tt,:,:)), squeeze(ole.Temp_Fix(tt,:,:)), 40);
    colorbar;
    xlim([nanmin(ole.Dist(:)) nanmax(ole.Dist(:))/8]);
    pause
end