% overlay bathymetry image for JHU APL AVHRR SST

load 'data/etopo2_extract.mat'

%%
ix = find_approx(topo.x, -76, 1):find_approx(topo.x,-65,1);
iy = find_approx(topo.y, 35, 1):find_approx(topo.y, 44, 1);

hf = figure; hax = gca; maximize;
insertAnnotation('OverlayBathy.m');
%m_proj('mercator','lon', [min(topo.x(ix)) max(topo.x(ix))], ...
%       'lat',[min(topo.y(iy)) max(topo.y(iy))]);
contour(topo.x(ix), topo.y(iy), topo.z(ix,iy)',[0 -100], 'k', 'LineWidth', 1.5);
hold on
[~,hc] = contour(topo.x(ix), topo.y(iy), topo.z(ix,iy)',-[50 80 2000 4000 5000], ...
                 'Color', [1 1 1]*0.6 , 'LineWidth', 1.5);
%clabel(cc,hh);
hax.Units = 'pixels';
hax.Position = [50 50 747 792] - [0 0 6 5];
hax.XTickLabels = {};
hax.YTickLabels = {};

%export_fig -painters images/overlay-bathy-non-transparent.png
hax.Visible = 'off';
export_fig -a4 -opengl -transparent images/overlay-bathy.png

% png't work as well
%system('convert images/overlay-bathy.png -transparent white images/overlay-bathy.png');

% off left = (80,30)
% bottom right = (825, 820);
% composite overlay-bathy.png wcr-avhrr-label.png -geometry +78+28 avhrr.png
% non-transparent version is actually smaller by 2 pixels?! +80+30 works ok