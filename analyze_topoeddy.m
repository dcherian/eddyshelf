% analyze topoeddy runs
% makes subplots of eddy tracks

%% read data
dir1 = 'runs/topoeddy/';
dir2 = {'runteb-02';'runteb-03';'runteb-04';'runteb-05';'runteb-06';'runteb-07'};
redo = 0;
for ii=1:length(dir2)
    dir3 = [dir1 dir2{ii} '/'];
    fname = [dir3 '/eddytrack.mat'];
    if ~exist(fname,'file') || redo == 1
        disp(['Calculating track for ' dir3]);
        track{ii} = track_eddy(dir3);
    else
        track{ii} = load(fname,'eddy');
        track{ii} = track{ii}.eddy;
    end
    params = read_params_from_ini([dir3 '/config/te_ini.nc']);
    track{ii}.sl_slope = params.bathy.sl_slope;
    track{ii}.L_slope  = params.bathy.L_slope;
end

%% make plot
nn = length(track);
if mod(nn,2) == 0
    a = factor(nn);
else
    a = factor(nn+1);
end

figure
for ii=1:nn
    subplot(a(1),a(2),ii);
    eddy = track{ii};
    vscale = eddy.Lz2;
    [c,h] = contour(eddy.xr/1000,eddy.yr/1000,eddy.h(2:end-1,2:end-1),10,'k');
    hold on; caxis([min(vscale(:)) max(vscale(:))]); colorbar;
    scatter(eddy.mx/1000,eddy.my/1000,10,vscale,'filled'); 
    plot(eddy.mx(1)/1000,eddy.my(1)/1000,'ko','MarkerSize',7);    
     axis image; xlim([0 220]);
    title(['runteb-0' num2str(ii+1) '| ' num2str(length(eddy.t)) ' days ' ...
            '| sl = ' num2str(eddy.sl_slope) ' | L = ' num2str(eddy.L_slope/1000)]);
    xlabel('X (km)');
    ylabel('Y (km)');
end
colormap(pmkmp);