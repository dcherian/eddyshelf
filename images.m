%% summary plot for thesis proposal

dirname = 'E:\Work\eddyshelf\runs\runte-03\';
fnames = ls([dirname '*_his*']);
tstep = [1 21 10 20];
clim = [-0.2 0.2];

fnames = {'ocean_his_0001.nc','ocean_his_0001.nc','ocean_his_0002.nc','ocean_his_0003.nc'};
vars = {'zeta'};
titles = {'\zeta','T'};
count = {[Inf Inf 1]; [Inf Inf 1 1]};

figure;
     
for jj=1:length(vars)
    for ii=1:4
        file = [dirname fnames{ii}];
        start = {[1 1 tstep(ii)]; [1 1 40 tstep(ii)]};
        zeta = ncread(file,vars{jj},start{jj},count{jj});
        time = ncread(file,'ocean_time');

        if ii == 1
            h = ncread(file,'h');
            grid1 = roms_get_grid(file,file,0,0);
        end

        subplot(2,2,ii)
        contourf(grid1.x_rho/1000,grid1.y_rho/1000,zeta',10,'LineWidth',1.5);
        caxis(clim);
        hold on;
        [c,hh] = contour(grid1.x_rho/1000,grid1.y_rho/1000,h',[200 750 1500],'b--','LineWidth',1.5);
        clabel(c,hh,'FontSize',12,'Color','b');
        xlabel(' X (km) ');
        ylabel(' Y (km) ');
        title([titles{jj} ' at t = ' num2str(time(tstep(ii))/86400) ' days']);
        colorbar
        beautify([14 16 18]);
    end
    pause
    %spaceplots([0 0 0 0], [0.02 0.02]);
    export_fig(['E:\Work\eddyshelf\images\' vars{jj} '.png']);
end

