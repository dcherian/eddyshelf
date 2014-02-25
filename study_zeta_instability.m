% dir0 = '/mnt/scylla/runew-06-scylla/';
% dirbig = '/mnt/scylla/runew-06-big/';
% dirbig2 = '/mnt/scylla/runew-06-big2/';
% dirsp = '/mnt/scylla/runew-06-bigsp/';
% dirper = '/mnt/scylla/spen/runew-03-per/';
% 
% rgrid0 = roms_get_grid([dir0 '/ocean_his.nc'],[]);
% rgridbig = roms_get_grid([dirbig '/ocean_his.nc'],[]);
% rgridbig2 = roms_get_grid([dirbig2 '/ocean_his.nc'],[]);
% rgridsp = roms_get_grid([dirsp '/ocean_his.nc'],[]);
% rgridper = roms_get_grid([dirper '/ocean_his.nc'],[]);
% 
% IX0 = size(rgrid0.x_rho,2);IY0 = size(rgrid0.y_rho,1);
% IXbig = size(rgridbig.x_rho,2);IYbig = size(rgridbig.y_rho,1);
% IXper = size(rgridper.x_rho,2);IYper = size(rgridper.y_rho,1);
% 
% %% read zeta
% zeta0 = dc_roms_read_data(dir0,'zeta',[],{'x' IX0-80 ...
%                     IX0-80; 'y' IY0-80 IY0-80},[],rgrid0);
% zetabig = dc_roms_read_data(dirbig,'zeta',[],{'x' IXbig-80 ...
%                     IXbig-80; 'y' IYbig-80 IYbig-80},[],rgridbig);
% zetabig2 = dc_roms_read_data(dirbig2,'zeta',[],{'x' IXbig-80 ...
%                     IXbig-80; 'y' IYbig-80 IYbig-80},[],rgridbig2);
% zetasp = dc_roms_read_data(dirsp,'zeta',[],{'x' IX0-80 ...
%                     IX0-80; 'y' IY0-80 IY0-80},[],rgridsp);
% %zetaper = dc_roms_read_data(dirper,'zeta',[],{'x' IXper-80 ...
% %                    IXper-80; 'y' IYper-80 IYper-80},[],rgridper);
%                 
%                 
% %% plot
% 
% figure;
% hold on
% plot(zeta0 - zeta0(1),'b');
% plot(zetasp - zetasp(1),'r');
% plot(zetabig - zetabig(1),'k');
% plot(zetabig2 - zetabig2(1),'c');
% %plot(zetaper - zetaper(1),'m');
% legend('original','bigsp','big','big2');


dirs = {'runs/topoeddy/spen/runew-03/';
       'runs/topoeddy/runew-03/'};
clear grd dz
for ii=1:length(dirs)
   grd(ii) = roms_get_grid([dirs{ii} '/ocean_avg.nc'],[]);
   dz{ii} = permute(diff(grd(ii).z_uw,1,1),[3 2 1]);
end

%%
figure; maximize;
suplabel('t', 'color = ubar (m/s) | contours = zeta ');
for tt=1:120    
    for ii=1:2%length(dirs)
        u = dc_roms_read_data(dirs{ii},'u',tt,{},[],grd(ii));
        zeta = dc_roms_read_data(dirs{ii},'zeta',tt,{},[],grd(ii));
        
        U = sum( u .* dz{ii},3) ./ avg1(grd(ii).h',1);
        
        subplot(1,2,ii);
        if tt == 1
            if ii == 1
                zmin = min(zeta(:));
                zmax = max(zeta(:));
            end
            hu(ii) = pcolorcen(U'); hold on
            caxis([-1 1] * 0.07); colorbar;
            [~,hc(ii)] = contour(avg1(zeta,1)', ...
                    linspace(zmin,zmax,15),'k');
            xlim([0 400]);
            ylim([0 212]);
            linex(400-40,[],'b');liney(180-40,[],'b');
            linex(400-40,[],'k'); liney(210-40,[],'k');
            liney(find_shelfbreak([dirs{ii} '/ocean_avg.nc'])/1500,[],'g');
        else
            set(hu(ii),'CData',U');
            set(hc(ii),'ZData',zeta');
        end
        title([dirs{ii} '| ' num2str(tt)]);
    end
    export_fig(sprintf('videos/zeta-instability/%03d.png',tt));
    pause(0.01);
end

%% 
zeta = run3.zeta;

figure;
hold on;
plot(squeeze(run3.zeta(60,20,:)),'r');
plot(squeeze(run3.zeta(320,20,:)),'r--');

plot(squeeze(run3.zeta(60,40,:)),'g');
plot(squeeze(run3.zeta(320,40,:)),'g*-');


plot(squeeze(run3.zeta(5,160,:)),'k');
plot(squeeze(run3.zeta(390,160,:)),'k-.');

plot(squeeze(run3.zeta(320,120,:)),'b--');
plot(squeeze(run3.zeta(60,120,:)),'b');
legend('(60,20)','(320,20)', '(60,40)','(320,40)',...
       '(5,16k0)','(390,160)','(320,120)','(60,120)','Location','SouthWest');
