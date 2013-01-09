
mm_instance = mm_setup;
mm_instance.pixelSize = [1600 900];
mm_instance.outputFile = 'mm_initial-y-04.avi';
mm_instance.ffmpegArgs = '-q:v 1 -g 1';
mm_instance.InputFrameRate = 3;
mm_instance.frameRate = 3;

dirname = 'runs/runte-05-rst/';
fnames  = ls([dirname '/*his*.nc']);

grid1 = roms_get_grid([dirname fnames(1,:)]);

for ff = 1:size(fnames,1)
    file = [dirname fnames(ff,:)];
    t = length(ncread(file,'ocean_time'));
    for i=1:t
        mod_movie(file,'dye_02',[i i],{},'s',39,'caxis([1e3 4e5])');
        hold on;
        [c,h] = contour(grid1.x_rho/1000,grid1.y_rho/1000,grid1.h,'k');
        clabel(c,h,108);
        mm_addFrame(mm_instance,gcf);
    end
end

mm_render(mm_instance)