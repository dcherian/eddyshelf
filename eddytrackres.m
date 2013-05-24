% Compares track variability with resolution changes

function [] = eddytrackres()

    dir1 = 'runs/topoeddy/';
    str  = 'runteb-04-*';
    
    names = ls([dir1 str]);
    
    colors = distinguishable_colors(size(names,1));
    
    figure; hold on
    for ii=1:size(names,1)
        fname = [dir1 strtrim(names(ii,:)) '/eddytrack.mat'];
        load(fname,'eddy');
        legstr{ii} = make_legend(eddy);
        plot(eddy.mx/1000,eddy.my/1000,'Color',colors(ii,:),'LineWidth',2);
    end
    
    xlim([0 180]);
    axis image;
    legend(legstr);
    

function [legend] = make_legend(track)

    dx = (track.xr(2,1)-track.xr(1,1))/1000;
    dy = (track.yr(1,2)-track.yr(1,1))/1000;
    legend = [num2str(dx) ' km x ' num2str(dy) ' km'];
    
% old code
%     lores = load('runs/topoeddy/runteb-04-lores/eddytrack.mat','eddy');
%     lores2 = load('runs/topoeddy/runteb-04-lores-2/eddytrack.mat','eddy');
%     hires = load('runs/topoeddy/runteb-04/eddytrack.mat','eddy');
%     hires2 = load('runs/topoeddy/runteb-04-hires-2/eddytrack.mat','eddy');
%     hires3 = load('runs/topoeddy/runteb-04-hires-2/eddytrack.mat','eddy');
%     
%     lores = lores.eddy; hires = hires.eddy; lores2 = lores2.eddy; hires2 = hires2.eddy;
% 
%     lores.legend = make_legend(lores); lores2.legend = make_legend(lores2);
%     hires.legend = make_legend(hires); hires2.legend = make_legend(hires2);
    
%     %dx = hires.xr(2,1)-hires.xr(1,1);
%     %dy = hires.yr(1,2)-hires.yr(1,1);
%     figure;
%     %subplot(2,2,[1 3])
%     hold on
%     plot(lores.mx/1000,lores.my/1000,'b',hires.mx/1000,hires.my/1000,'r',...
%          lores2.mx/1000,lores2.my/1000,'c',hires2.mx/1000,hires2.my/1000,'m');
%     %plot(lores.cx/1000,lores.cy/1000,'b--',hires.cx/1000,hires.cy/1000,'r--',...
%     %     lores2.cx/1000,lores2.cy/1000,'c--',hires2.cx/1000,hires2.cy/1000,'m--');
%     contour(hires.xr/1000,hires.yr/1000,hires.h(2:end-1,2:end-1),10,'k');
%     legend(lores.legend,hires.legend,lores2.legend,hires2.legend);
%     xlim([0 180]); axis image
%     subplot(222)
%     plot(lores.t,(lores.mx-hires.mx)/dx,'b',lores.t,(lores2.mx-hires.mx)/dx,'r');
%     legend('lores','lores-2');
%     title('diff. in x - number of points - hires');
%     subplot(224)
%     plot(lores.t,(lores.my-hires.my)/dy,'b',lores.t,(lores2.my-hires.my)/dy,'r');
%     legend('lores','lores-2');
%     title('diff in y - number of points - hires');