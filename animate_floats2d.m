function [] = animate_floats2d(rgrid,zeta,floats)

    % need to color floats based on starting position
   % colors = distinguishable_colors(size(floats.x,2));

    cmap = flipud(cbrewer('div', 'BrBG', 32));

    figure;
    % plot bathymetry
%     ax(2) = subplot(132);
%     plot(rgrid.x_rho(1,:)/1000,-rgrid.h(1,:));
%     liney(-114,'-114');
%     liney(-130,'-130');
%     beautify
%     ax(3) = subplot(133);
%     plot3(floats.x/1000,floats.y/1000,floats.z);
%     hold on;
    
    ind0 = find(rgrid.ocean_time == floats.time(1));

    for i=1:2:size(floats.x,1)
        %ax(1) = subplot(131);
        cla
        contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,ind0+i)');
        shading flat; axis image
        colormap(cmap); 
        caxis([min(zeta(:)) max(zeta(:))]);
        hold on
        [C,hc] = contour(rgrid.x_rho./1000,rgrid.y_rho./1000,rgrid.h,[114 500 750 1100],'k');
        clabel(C,hc);
        %plot(floats.init(:,1)/1000,floats.init(:,2)/1000,'x','MarkerSize',12);
        for j = 1:size(floats.x,2)
              plot(floats.x(1:floats.fac*i,j)/1000,floats.y(1:floats.fac*i,j)/1000, ...
                     'Color','k','MarkerSize',10);
    %     scatter3(floats.x(1:fac*i,j),floats.y(1:fac*i,j),floats.z(1:fac*i,j), ...
    %                 10,floats.z(1:fac*i,j));
        end
        %linkaxes(ax,'x')
        title(['t = ' num2str(floats.time(1)+i) ' days']);
        pause(0.01)
    end
end