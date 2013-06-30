function [] = plot_coast(mm,flag)
    load topo_useast
    
    if ~exist('mm','var'), mm = 1; end
    if ~exist('flag','var'), flag = 1; end
    
    [cc,hh] = contour(Xuseast(1:mm:end,1:mm:end), ...
                      Yuseast(1:mm:end,1:mm:end), ...
                      Zuseast(1:mm:end,1:mm:end), ...
                      [100 200 500 1000 2000 3000 4000],'k');
    if flag, clabel(cc,hh); end
    Z_dar
    hold on
    plot(coast(1:mm:end,1),coast(1:mm:end,2),'k')
    ylim([34.5 43]);
end