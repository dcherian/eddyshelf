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
    m_coast('patch',[0 0 0],'edgecolor','none');
    m_grid
    
    % plot bathy
    load 'E:\Work\eddyshelf\data\etopo2_extract.mat';
    levels = [-50 -100 -1000 -2000 -3000 -4000];
    [c,h] = m_contour(topo.x,topo.y,topo.z',levels,'k-');
    clabel(c,h,'LabelSpacing',500);
    
    % plot eddy tracks
    for i=1:nunique
        maski = fillnan(double(data.track == utrack(i)),0);
        index_first = find(maski == 1,1,'first');
        index_last = find(maski == 1,1,'last');
        
        m_plot(data.lon .* maski,data.lat .* maski,'Color',colors(i,:));
        
        % isolate first and last points to plot circle
        maskt = nan(size(maski));        
        maskt(index_first) = 1; maskt(index_last) = 1; 
        maski = maski .* maskt;
        
        xCenter = data.lon .* maski;
        yCenter = data.lat .* maski;
        theta = 0 : 0.01 : 2*pi;
        theta = repmat(theta,length(maski),1);
        radius = (data.L .* maski)./100/1000;
        x = bsxfun(@plus,bsxfun(@times,radius,cos(theta)),xCenter);
        y = bsxfun(@plus,bsxfun(@times,radius,sin(theta)),yCenter);
        %m_plot(x, y,'.','Color',colors(i,:),'MarkerSize',6);
        
        % plot 'x' at initial position
        m_plot(data.lon(index_first),data.lat(index_first),'*','Color',colors(i,:));
    end
    
    beautify([14 14 16]);
    