% analyze topoeddy runs
% makes subplots of eddy tracks

%% read data
dir1 = 'runs/topoeddy/';
dir2 = {'runteb-02';'runteb-03';'runteb-04';'runteb-05';'runteb-06'};
redo = 1;
for ii=1:length(dir2)
    dir3 = [dir1 dir2{ii} '/'];
    fname = [dir3 '/eddytrack.mat'];
    if ~exist(fname,'file') || redo == 1
        disp(['Calculating track for ' dir3]);
        track{ii} = track_eddy(dir3);
    else
        track{ii} = load(fname,'eddy');
    end
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
    plot(eddy.mx/1000,eddy.my/1000,'k','LineWidth',2); hold on;
    [c,h] = contour(eddy.xr/1000,eddy.yr/1000,eddy.h(2:end-1,2:end-1),10,'k');
    clabel(c,h);
    xlim([0 220]);
    title(['runteb-0' num2str(ii+1)]);
    xlabel('X (km)');
    ylabel('Y (km)');
end