fname = '../data/All_Oleander_3D.mat';

ole = load(fname);

%% churchill (1986) section
% Oct 19-20, 1983
% stored in reverse order, that's why I reverse axis direction.
PlotOleanderSection(114, [11 17]);
export_fig -r150 -a2 images/oleander-oct1983.png

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