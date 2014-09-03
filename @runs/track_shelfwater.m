% track where on the shelf the water exported across
% the shelfbreak is coming from
function [] = track_shelfwater(runs)
    bins = runs.csflux.west.shelfwater.bins;
    time = runs.csflux.time/86400;

    % check
    if any(avg1(bins) > runs.bathy.xsb)
        ind = find(bins > runs.bathy.xsb, 1, 'first');
    else
        ind = length(bins);
    end

    % in Sv and filtered with above check
    trans = squeeze(runs.csflux.west.shelfwater.trans(:,1,1: ...
                                                      ind)/1e6);

    % grid matrices
    bmat = repmat(bins(1:ind),[length(time) 1])./runs.rrshelf;
    tmat = repmat(time', [1 size(bmat,2)]);

    figure;
    pcolorcen(bmat, tmat, trans);
    ylabel('Time days');
    xlabel('rossby radii from shelfbreak');
    colorbar;
    cblabel('Sv');
    title(runs.name);
end
