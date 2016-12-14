% use the 2D sections at each location to partition fluxes as I
% want to post-calculation of actual fluxes

function [] = process_asfluxes(runs)

    assert(~isempty(runs.asflux));

    yvec = runs.asflux.yvec;
    sy1 = runs.asflux.iy(1);

    % correct for sy1 = isb by default
    isl = runs.bathy.isl - sy1 + 1;
    % these ranges are correct! even if you don't think so.
    % ikeflux = deep.ikeflux + topo.ikeflux
    ideep = [isl:length(yvec)];
    itopo = [1:isl];

    % rename variables
    ikefluxyt = runs.asflux.ikefluxyt;
    ipefluxyt = runs.asflux.ipefluxyt;
    ikefluxyt_edd = runs.asflux.eddy.ikefluxyt;
    ipelfuxyt_edd = runs.asflux.eddy.ipefluxyt;

    % deep water
    runs.asflux.deep.ikeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                                                   ikefluxyt(ideep,:,ii), 1));
    runs.asflux.deep.ipeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                                                   ipefluxyt(ideep,:,ii), 1));
    runs.asflux.deep.eddy.ikeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                                                      ikefluxyt_edd(ideep,:,ii), 1));
    runs.asflux.deep.eddy.ipeflux(:,ii) = squeeze(trapz(yvec(ideep), ...
                                                      ipefluxyt_edd(ideep,:,ii), 1));

    % over shelf-slope
    runs.asflux.topo.ikeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                                                   ikefluxyt(itopo,:,ii), 1));
    runs.asflux.topo.ipeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                                                   ipefluxyt(itopo,:,ii), 1));
    runs.asflux.topo.eddy.ikeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                                                      ikefluxyt_edd(itopo,:,ii), 1));
    runs.asflux.topo.eddy.ipeflux(:,ii) = squeeze(trapz(yvec(itopo), ...
                                                      ipefluxyt_edd(itopo,:,ii), 1));

end