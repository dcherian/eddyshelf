
mm_instance = mm_setup;
mm_instance.pixelSize = [1600 900];
mm_instance.outputFile = 'mm_initial-pt-05.avi';
mm_instance.ffmpegArgs = '-q:v 1 -g 1';
mm_instance.InputFrameRate = 3;
mm_instance.frameRate = 3;

dirname = 'runs/runte-05-rst/';
varname = 'dye_02';
fnames  = ls([dirname '/*his*.nc']);

caxis_zeta = [-0.15 0.15];
caxis_dye_02 = [0 400000];

contour_zeta = [-0.15:0.015:0.15];
contour_dye_02 = linspace(0,400000,20);

grid1 = roms_get_grid([dirname fnames(1,:)]);

for ff = 1:size(fnames,1)
    file = [dirname fnames(ff,:)];
    if ff == 1 && strcmp(varname,'dye_02');
        read_start = [1 1 39 61];
        read_count = [Inf Inf 1 Inf];
    else
        read_start = [1 1 39 1];
        read_count = [Inf Inf 1 Inf];
    end
    time = ncread(file,'ocean_time',read_start(end),read_count(end));
    var  = ncread(file,varname,read_start,read_count);
    zeta = ncread(file,'zeta',[read_start(1) read_start(2) read_start(end)], ...
                        [read_count(1) read_count(2) read_count(end)]);
    
    t = length(time);
    for i=1:t
        %mod_movie(file,'zeta',[i i],{},'s',39,'caxis([-0.15 0.15])');
        clf
        eval(['contourf(grid1.x_rho/1000,grid1.y_rho/1000,var(:,:,i)'',contour_' varname ');']);
        set(gcf,'renderer','zbuffer'); shading flat;
        colorbar;
        eval(['caxis(caxis_' varname ');']);
        hold on;
        [c,h] = contour(grid1.x_rho/1000,grid1.y_rho/1000,grid1.h,[200 500 1000 1500],'k','LineWidth',2);
        clabel(c,h,'LabelSpacing',108*2);
        contour(grid1.x_rho/1000,grid1.y_rho/1000,zeta(:,:,i)',10,'w','LineWidth',2);
        
        xlabel('X (km)');
        ylabel('Y (km)');
        
        title([varname ' | t = ' num2str(time(i)/86400) ' days']);
        beautify([14 14 16]);
        mm_addFrame(mm_instance,gcf);
    end
end

mm_render(mm_instance)