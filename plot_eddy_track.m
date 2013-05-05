function [] = plot_eddy_track(eddy,index)

    hold on
    plot(eddy.mx(1:index)/1000,eddy.my(1:index)/1000,'r');
    plot(eddy.mx(index)/1000,eddy.my(index)/1000,'r.','MarkerSize',16);
    contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,index),1,'k');