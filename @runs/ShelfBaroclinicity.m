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

    vsurf = permute( ...
        dc_roms_read_data(runs, 'v', [start stop], ...
                          {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' runs.bathy.isb runs.bathy.isb; ...
                        'z' runs.rgrid.N runs.rgrid.N}), [1 3 2]);

    vbot = permute( ...
        dc_roms_read_data(runs, 'v', [start stop], ...
                          {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' runs.bathy.isb runs.bathy.isb; ...
                        'z' zbot zbot}), [1 3 2]);

    csdsurf = dc_roms_read_data(runs, runs.csdname, [start stop], ...
                                {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' 1 runs.bathy.isb; ...
                        'z' runs.rgrid.N runs.rgrid.N});
    csdbot = dc_roms_read_data(runs, runs.csdname, [start stop], ...
                               {'x'  runs.spng.sx1 runs.spng.sx2; ...
                        'y' 1 runs.bathy.isb; ...
                        'z' zbot zbot});

    % mask of 1s is NaN-ed out
    nonshelfmask = (csdsurf > runs.bathy.xsb) | (csdbot > runs.bathy.xsb);

    % find most on-shore extent of eddy water, mask out everything north of that
    % this criterion is too restrictive, lose a lot of points
    % yind = bsxfun(@times, fillnan(double(nonshelfmask),0), [1:70]);
    % nonshelfmask(:,nanmin(yind(:)):end,:) = 1;

    % look only for negative velocities
    upos = (usurf > 0) | (ubot > 0);
    nonshelfmask = (nonshelfmask | upos);

    % us = usurf; ub = ubot;
    % nanmask = nonshelfmask;
    % us(nanmask) = NaN;
    % ub(nanmask) = NaN;

    % uabs = abs(us);
    % magmask = bsxfun(@lt, uabs, 0.2 * nanmax(uabs(:)));
    % bc = abs(us - ub)./uabs; % (us - ub)/(|us| + |ub|)
    % bc(magmask) = NaN; % mask out low velocity points
    % bcmn = squeeze(nanmedian(nanmedian(bc,1),2));
    % bc = nanmean(bcmn)

    [shelfbc.shelf, shelfbc.thresh] = baroclinicity(usurf, ubot, nonshelfmask);
    shelfbc.nonshelf = baroclinicity(usurf, ubot, ~nonshelfmask);

    shelfbc.sbreak.shelf = baroclinicity(vsurf, vbot, nonshelfmask(:,end,:));
    shelfbc.sbreak.nonshelf = baroclinicity(vsurf, vbot, ~nonshelfmask(:,end,:));

    shelfbc.time = runs.csflux.time(start:stop);
    shelfbc.tind = [start stop];

    hash = githash([mfilename('fullpath') '.m']);
    shelfbc.hash = hash;

    runs.shelfbc = shelfbc;
    save([runs.dir '/shelfbc.mat'], 'shelfbc');
    toc(ticstart);
end

function [bcmn, thresh] = baroclinicity(usurf, ubot, nanmask)
    usurf(nanmask) = NaN;
    ubot(nanmask) = NaN;

    thresh = [0.1 0.2 0.3 0.4];

    uabs = abs(usurf);
    for ii=1:length(thresh)
        bc = abs(usurf - ubot)./uabs; % (us - ub)/(|us| + |ub|)
        bc(uabs < thresh(ii) * max(uabs(:))) = NaN; %mask out low velocity points
        bcmn(:,ii) = squeeze(nanmedian(nanmedian(bc,1),2));
    end
end