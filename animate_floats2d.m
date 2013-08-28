function [] = animate_floats2d(rgrid,zeta,floats)

    % need to color floats based on starting position
   % colors = distinguishable_colors(size(floats.x,2));

    cmap = flipud(cbrewer('div', 'BrBG', 32));

    figure;
    
    tfilt = cut_nan(fillnan(floats.init(:,4),0));
    
    ind0 = find_approx(rgrid.ocean_time, tfilt(1),1);
        
    for i=ind0:size(zeta,3)
        cla
        contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,i)');
        shading flat; axis image
        colormap(cmap); 
        caxis([min(zeta(:)) max(zeta(:))]);
        hold on
        [C,hc] = contour(rgrid.x_rho./1000,rgrid.y_rho./1000,rgrid.h,[114 500 750 1100],'k');
        clabel(C,hc);
        plot(floats.init(:,1)/1000,floats.init(:,2)/1000,'x','MarkerSize',12);
        floats.fac = floor(floats.fac);
        nn = find_approx(floats.time,rgrid.ocean_time(i),1);
        plot(floats.x(nn,:)/1000,floats.y(nn,:)/1000,'k.','MarkerSize',10);
        title(['t = ' num2str(floats.time(1)+i) ' days']);
        pause(0.01)
    end
end