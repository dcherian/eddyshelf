function [] = plot_chelton(data,mask)

    data = apply_mask(data,mask);
    
    % figure out unique tracks
     utrack = unique(data.track);
    nunique = length(utrack);
  
    colors = distinguishable_colors(nunique);
    
    figure; hold on
    m_proj('mercator','longitudes',[-76 -60], 'latitudes',[32 46]);
    m_coord('geographic');
    
    % plot coast
    %m_coast('patch',[0 0 0],'edgecolor','none');
    m_grid
    
    % plot bathy
    load 'E:\Work\eddyshelf\data\etopo2_extract.mat';
    levels = [-50 -100 -1000 -2000 -4000];
    [c,h] = m_contour(topo.x,topo.y,topo.z',levels,'k-');
    clabel(c,h,'LabelSpacing',500);
    
    % plot eddy tracks
    for i=1:nunique
        maski = fillnan(double(data.track == utrack(i)),0);
        m_plot(data.lon .* maski,data.lat .* maski,'Color',colors(i,:));
    end
    