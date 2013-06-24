% Compares track variability with resolution changes

function [] = eddytrackres()

    dir1 = 'runs/topoeddy/';
    str  = 'runteb-04-*';
    
    names = ls([dir1 str]);
    
    plots = 0;
    
    colors = distinguishable_colors(size(names,1));
    % generate gray colors
    %arr = linspace(0,200,size(names,1))/255;
    %colors = repmat(arr',1,3);
    
    figure; hold on
    for ii=1:size(names,1)
        files{ii} = [dir1 strtrim(names(ii,:))];
        fname = [files{ii} '/eddytrack.mat'];
        load(fname,'eddy');
        eddies{ii} = eddy;
        legstr{ii} = make_legend(eddy);
    end
    [sorted,sort_ind] = sort(legstr);
    
    for ii=1:length(sort_ind)
        kk = sort_ind(ii);
        eddy = eddies{kk};
        aa = 5; bb = aa * 2;
        %if ii == 2, eddy.cy = eddy.cy-50000; end
        subplot(aa,2,[1:2:bb]); hold on
        plot(eddy.cx/1000,eddy.cy/1000,'Color',colors(ii,:),'LineWidth',2);
        subplot(aa,2,2); hold on
        plot(eddy.t,eddy.amp,'Color',colors(ii,:));
        subplot(aa,2,4); hold on
        plot(eddy.t,eddy.dia/1000,'Color',colors(ii,:));
        subplot(aa,2,6); hold on
        plot(eddy.t,eddy.cx/1000,'Color',colors(ii,:));
        subplot(aa,2,8); hold on
        plot(eddy.t,eddy.cy/1000,'Color',colors(ii,:));
        subplot(aa,2,10); hold on
        plot(eddy.t,eddy.Lz2,'Color',colors(ii,:));
    end
    
    xtick = 0:25:300;
    subplot(aa,2,[1:2:bb]);
    hold on;
    [cc,hh] = contour(eddy.xr/1000,eddy.yr/1000,eddy.h(2:end-1,2:end-1),[500 1000 2000 2300],'k');
    clabel(cc,hh);
    axis image;
    legend(sorted,'Location','NorthWestOutside');
    title('center track | red crosses at t=150 days');
    xlim([0 180]);
    xlabel('x (km)'); ylabel('y (km)');
    for kk = 1:length(sort_ind)
        eddy = eddies{kk};
        subplot(aa,2,[1:2:bb]);
        index = find_approx(eddy.t,150);
        plot(eddy.cx(index)/1000,eddy.cy(index)/1000,'rx','MarkerSize',12);
    end
    beautify
    
    ax(1) = subplot(aa,2,2);
    ylabel('amplitude (m)');
    xlabel('time (days)');
    set(gca,'XTick',xtick);
    beautify
    
    ax(2) = subplot(aa,2,4);
    ylabel('diameter (km)');
    xlabel('time (days)');
    set(gca,'XTick',xtick);
    beautify
    
    ax(3) = subplot(aa,2,6);
    ylabel('x - center (km)');
    xlabel('time (days)');
    set(gca,'XTick',xtick);
    beautify
    
    ax(4) = subplot(aa,2,8);
    ylabel('y - center (km)');
    xlabel('time (days)');
    set(gca,'XTick',xtick);
    beautify
    
    ax(5) = subplot(aa,2,10);
    ylabel('vertical scale (m)');
    xlabel('time (days)');
    set(gca,'XTick',xtick);
    beautify
    
    linkaxes(ax,'x');
    if plots
        set(gcf,'position',[0 0 1600 900]);
        export_fig eddytrackres.png
    end

    % compare both 0.5km runs 
    figure; hold on
    II = 2;
    dy = -(eddies{sort_ind(II)}.cy(II) - eddies{sort_ind(II+1)}.cy(II))/1000;
    plot(eddies{sort_ind(II)}.cx/1000,eddies{sort_ind(II)}.cy/1000,'b');
    plot(eddies{sort_ind(II+1)}.cx/1000,eddies{sort_ind(II+1)}.cy/1000 - dy, 'r');
    plot(eddies{sort_ind(II+1)}.cx/1000,eddies{sort_ind(II+1)}.cy/1000, 'k');
    [cc,hh] = contour(eddy.xr/1000,eddy.yr/1000,eddy.h(2:end-1,2:end-1),[500 1000 2000 2300],'k');
    clabel(cc,hh); 
    axis image;xlim([0 220]);
    xlabel('x (km)'); ylabel('y (km)');
    title('comparing both 0.5km x 3km runs');
    legend('center of domain','higher-displaced','higher');
    
    if plots
        set(gcf,'position',[0 0 1600 900]);
        export_fig compare500.png
    end
    
    %analyze_cyclones(eddies,legstr,sort_ind,files,70.5);
    
    if plots
        set(gcf,'position',[0 0 1600 900]);
        export_fig cyclones.png
    end
    
function [] = analyze_cyclones(eddies,legstr,sort_ind,files,t0)
    % x-locations of centers seen to diverge then after inspecting image 
    % created by eddytrackres()
    N = length(eddies);
    a = factor(N-1);
    clim = linspace(-0.02,0.05,35);
    clim2 = linspace(-2e-5,2e-5,10);
    clim3 = linspace(-5e-6,5e-6,8);
    
    h1 = figure;
    h2 = figure;
    h3 = figure;
    for ii=1:N-1
        % read data
        kk = sort_ind(ii);
        tind = find(eddies{kk}.t == t0);
        if isempty(tind), continue; end
        rgrid = roms_get_grid([files{kk} '/ocean_avg.nc'],[files{kk} '/ocean_avg.nc']);
        zeta = double(ncread([files{kk} '/ocean_avg.nc'],'zeta',[1 1 tind],[Inf Inf 1]));
        u = double(ncread([files{kk} '/ocean_avg.nc'],'u',[1 1 rgrid.N tind],[Inf Inf 1 1]));
        v = double(ncread([files{kk} '/ocean_avg.nc'],'v',[1 1 rgrid.N tind],[Inf Inf 1 1]));
        w = double(ncread([files{kk} '/ocean_avg.nc'],'w',[1 1 rgrid.N tind],[Inf Inf 1 1]));
        rgrid.xv = rgrid.x_v'; rgrid.yv = rgrid.y_v';
        rgrid.yu = rgrid.y_u'; rgrid.zv = rgrid.z_r;
        [vor,xvor,yvor,~] = vorticity_cgrid(rgrid,u,v);
        
        % compare zeta
        figure(h1)
        ax1(ii) = subplot(a(1),a(2),ii);
        kk = sort_ind(ii);
        [cc,hh] = contourf(eddies{kk}.xr/1000,eddies{kk}.yr/1000,zeta(2:end-1,2:end-1),clim);
        shading flat;
        %clabel(cc,hh);
        title(legstr{kk})
        xlabel('x (km)'); ylabel('y (km)');
        hold on
        plot(eddies{kk}.cx/1000,eddies{kk}.cy/1000,'w');
        plot(eddies{kk}.cx(tind)/1000,eddies{kk}.cy(tind)/1000,'b.','MarkerSize',16);
        %axis image; 
        colorbar
        beautify
        caxis([clim(1) clim(end)]);
        
        % compare surface vorticity
        figure(h2)
        ax2(ii) = subplot(a(1),a(2),ii);
        [cc,hh] = contourf(xvor/1000,yvor/1000,vor,clim2); shading flat
        %clabel(cc,hh);
        title(legstr{kk})
        xlabel('x (km)'); ylabel('y (km)');
        hold on
        plot(eddies{kk}.cx/1000,eddies{kk}.cy/1000,'w');
        plot(eddies{kk}.cx(tind)/1000,eddies{kk}.cy(tind)/1000,'b.','MarkerSize',16);
        %axis image; 
        colorbar
        beautify
        caxis([clim2(1) clim2(end)]);
        
        % compare w
        figure(h3)
        ax3(ii) = subplot(a(1),a(2),ii);
        kk = sort_ind(ii);
        [cc,hh] = contourf(eddies{kk}.xr/1000,eddies{kk}.yr/1000,w(2:end-1,2:end-1),clim3);
        shading flat
        %clabel(cc,hh);
        title(legstr{kk})
        xlabel('x (km)'); ylabel('y (km)');
        hold on
        plot(eddies{kk}.cx/1000,eddies{kk}.cy/1000,'w');
        plot(eddies{kk}.cx(tind)/1000,eddies{kk}.cy(tind)/1000,'b.','MarkerSize',16);
        %axis image; 
        colorbar
        beautify
        caxis([clim3(1) clim3(end)]);
    end
    figure(h1)
    [~,hh1] = suplabel([' zeta | t_0 = ' num2str(t0) ' days'],'t');
    set(hh1,'FontSize',12);
    linkaxes(ax1,'xy');
    
    figure(h2)
    [~,hh1] = suplabel([' surface vor. | t_0 = ' num2str(t0) ' days'],'t');
    set(hh1,'FontSize',12);
    linkaxes(ax2,'xy');
    
    figure(h2)
    [~,hh1] = suplabel([' w | t_0 = ' num2str(t0) ' days'],'t');
    set(hh1,'FontSize',12);
    linkaxes(ax3,'xy');

function [legend] = make_legend(track)

    dx = (track.xr(2,1)-track.xr(1,1))/1000;
    dy = (track.yr(1,2)-track.yr(1,1))/1000;
    legend = sprintf('%.2f km x %.2f km',dx,dy);
    
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