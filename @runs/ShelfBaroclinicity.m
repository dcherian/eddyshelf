function [] = ShelfBaroclinicity(runs)

    ticstart = tic;

    disp(['ShelfBaroclinicity(' runs.name ')']);

    [start, stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1));

    if runs.params.misc.rdrg == 0
        zbot = 1;
    else
        zbot = 10;
    end

    usurf = dc_roms_read_data(runs, 'u', [start stop], ...
                              {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' 1 runs.bathy.isb; ...
                        'z' runs.rgrid.N runs.rgrid.N});

    ubot = dc_roms_read_data(runs, 'u', [start stop], ...
                             {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' 1 runs.bathy.isb; ...
                        'z' zbot zbot});

    csdsurf = dc_roms_read_data(runs, runs.csdname, [start stop], ...
                                {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' 1 runs.bathy.isb; ...
                        'z' runs.rgrid.N runs.rgrid.N});
    csdbot = dc_roms_read_data(runs, runs.csdname, [start stop], ...
                               {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' 1 runs.bathy.isb; ...
                        'z' zbot zbot});

    shelfbc.shelf = baroclinicity(usurf, ubot, csdsurf > runs.bathy.xsb);
    shelfbc.nonshelf = baroclinicity(usurf, ubot, csdsurf < runs.bathy.xsb);

    shelfbc.time = runs.csflux.time(start:stop);
    shelfbc.tind = [start stop];

    hash = githash([mfilename('fullpath') '.m']);
    shelfbc.hash = hash;

    runs.shelfbc = shelfbc;

    save([runs.dir '/shelfbc.mat'], 'shelfbc');
    toc(ticstart);
end

function [bcmn] = baroclinicity(usurf, ubot, nanmask)
    usurf(nanmask) = NaN;
    ubot(nanmask) = NaN;

    uabs = (abs(usurf) + abs(ubot));
    bc = abs(usurf - ubot)./uabs; % (us - ub)/(|us| + |ub|)
    bc(bsxfun(@lt, uabs, 0.3 * max(uabs, [], 3))) = NaN; % mask out low velocity points

    bcmn = squeeze(nanmedian(nanmedian(bc,1),2));
end