function [handles] = plotAvgStreamer(runArray, isobath, hax)

    if ~exist('isobath', 'var'), isobath = 1; end

    if ~exist('hax', 'var')
        figure; maximize;
        hax = gca;
        hold on;
    else
        axes(hax);
        hold on
    end
    insertAnnotation('runArray.plotAvgStreamer');

    if isempty(runArray.filter)
        runArray.filter = 1:runArray.len;
    end

    ii = 1;
    for ff = runArray.filter
        run = runArray.array(ff);
        name = runArray.getname(ff);

        [start,stop] = run.flux_tindices(isobath);
        [V0,L0,Lz0] = run.EddyScalesForFlux(start);

        handles.hplt(ii) = ...
            plot(run.streamer.xivec * 1e3./run.params.eddy.dia(1)/2, ...
                 run.streamer.xprof(:,isobath), ...
                 'DisplayName', name);
        ii = ii+1;
    end

    liney(0);
    handles.hleg = legend('Location', 'NorthWest');
    handles.isolabel = runArray.array(1).add_isobathlabel(isobath);
    ylabel('Vertically integrated transport (m^2/s)');
    xlabel('(Along-isobath distance from eddy center)/(Eddy radius)');
    beautify;
end