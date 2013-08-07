%% compare N-S & E-W isobaths

runns = runs('runs/topoeddy/runteb-04-hires-9');
%runew = runs('runs/topoeddy/runbathysouth-02/ocean_avg_001.nc');
% updated for run with closed boundaries, same resolution, and similar
% initial location - this is a much better comparision
runew = runs('runs/topoeddy/runbathysouth-03-new/');

fontsize = [14 14 16];
figure
subplot(121)
hold on
pcolorcen(runns.rgrid.xr/1000,runns.rgrid.yr/1000,runns.bathy.h);
xlabel('X (km)'); ylabel('Y (km)');
plot(runns.eddy.cx/1000,runns.eddy.cy/1000,'Color','b','LineWidth',2);
colorbar
clim = caxis;
plot(runns.eddy.we/1000,runns.eddy.cy/1000,'k');
axis image
title('Bathymetry (color), center (blue), west edge (black)');
beautify(fontsize);

subplot(122)
hold on
pcolorcen(runew.rgrid.xr/1000,runew.rgrid.yr/1000,runew.bathy.h);
xlabel('X (km)'); ylabel('Y (km)');
plot(runew.eddy.cx/1000,runew.eddy.cy/1000,'Color','b','LineWidth',2);
plot(runew.eddy.cx/1000,runew.eddy.se/1000,'k');
axis image
caxis(clim);
title('Bathymetry (color), center (blue), south edge (black)');
beautify(fontsize);

export_fig E:/Work/notes/research/eddyshelf/images/compare-ns-ew.png

%% convergence

% run eddytrackres.m

%% compare floats for runbathysouth-02

run = runs('runs/topoeddy/runbathysouth-02/ocean_avg_004.nc');
%run.eddy = load('runs/topoeddy/runbathysouth-02/eddytrack_004.mat','eddy');
%run.eddy = run.eddy.eddy;
run.ltrans = floats('ltrans','runs/topoeddy/runbathysouth-02/ltrans-compare.nc',run.rgrid);
run.ltrans.time = run.ltrans.time + 48*86400;
%%
figure
run.roms.plot_stats
run.ltrans.plot_stats
% 
% figure
% plot(run.roms.time,run.roms.N); hold on
% plot(run.ltrans.time,run.ltrans.N,'r');
% legend('roms','ltrans');

%% shelfbreak front movie

folder = ['runs/sbeddy/old expts/runsbe-04/'];
fname = [folder '/ocean_his_0001.nc'];
sst = roms_read_data(folder,'temp',[1 1 40 1],[Inf Inf 1 Inf]);
temp = roms_read_data(folder,'temp',[1 100 1 1],[Inf 1 Inf Inf]);
time = roms_read_data(folder,'ocean_time');
rgrid = roms_get_grid(fname,fname,0,1);

%%
fonts = [16 16 18];
makeVideo = 1;

if makeVideo
    aviobj = VideoWriter('sbreak1','MPEG-4');
    set(aviobj,'FrameRate',6,'Quality',100);
    open(aviobj);
end

tstr = [num2str(time(1)/86400) ' days'];
figure;
subplot(121)
[hsst] = pcolor(rgrid.x_rho/1000,rgrid.y_rho/1000,sst(:,:,1)');
shading interp;
colorbar; caxis([19 20.5]); liney(rgrid.y_rho(100,1)/1000,[]);
axis image; xlabel('X (km)'); ylabel('Y (km)');
ht1 = title(['SST | ' tstr]);
beautify(fonts);

xzr = repmat(rgrid.x_rho(1,:)',[1 rgrid.N]);
ax = subplot(122);
[~,htemp] = contourf(xzr/1000,squeeze(rgrid.z_r(:,100,:))',temp(:,:,1),24);
axis square;ylim([-800 0]); colorbar; 
ht2 = title(['Temp section | ' tstr]);
ylabel('Z (m)'); xlabel('X (km)');
beautify(fonts);
if makeVideo
    %shading(gca,'interp');
    disp('maximize!');
    pause; 
    
    spaceplots(gcf,0.1*ones([1 4]),0.1*ones([1 4]))
%               mm_addFrame(mm_instance,gcf);
    F = getframe(gcf);
    writeVideo(aviobj,F);
end
for tt=2:size(temp,3)
    tstr = [num2str(time(tt)/86400) ' days'];
    set(ht1,'String',['SST | ' tstr]);
    set(ht2,'String',['Temp section | ' tstr]);
    
    set(hsst,'CData',sst(:,:,tt)');
    set(htemp,'ZData',temp(:,:,tt));
    if makeVideo
       % shading(gca,'interp');
        %mm_addFrame(mm_instance,gcf);
        F = getframe(gcf);
        writeVideo(aviobj,F);
    end
    pause(0.01);
end

if makeVideo
   % mm_render(mm_instance);
   close(aviobj);
end