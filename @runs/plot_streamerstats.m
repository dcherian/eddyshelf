% plot streamer profiles
function [] = plot_streamerstats(runs)
    bins = runs.streamer.bins;
    figure
    subplot(121)
    cmap = brighten(cbrewer('seq','YlOrRd',runs.streamer.sz4dsp(2)),0);
    cmap = cmap(3:end,:,:); % chuck out lightest colors
    set(gca,'ColorOrder',cmap); colormap(cmap);
    line(runs.streamer.west.Vbin, repmat(avg1(bins'),[1 runs.streamer.sz4dsp(2)]));
    hold on
    zcenbin = vecfind(bins, cut_nan(runs.streamer.west.zcen));
    Vcenbin = diag(runs.streamer.west.Vbin(zcenbin,:));
    colorbar; cblabel('day');
    scatter(gca,zeros(20,1),runs.streamer.west.zcen, ...
            96,runs.streamer.time/86400,'filled');
    caxis([min(runs.streamer.time) max(runs.streamer.time)]/86400);
    xlabel('Volume (m^3)');
    ylabel(['Depth (' num2str(dbin) ' m bins)']);

    subplot(122); hold on
    plot(runs.streamer.time/86400,runs.streamer.west.zcen,'r');
    plot(runs.streamer.time/86400,runs.streamer.west.zdcen,'b');
    legend('z centroid','z-dye centroid');
    ylabel(' Depth (m) '); xlabel('day');
end
