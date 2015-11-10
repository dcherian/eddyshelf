% This function will calculate the extent to which an eddy penetrates the shelf.
% Could be done based on
%     1) Eddy dye
%     2) density
% At best of times, density envelope is the same as eddy dye envelope.
% At the worst of times, density envelope is zero.
function [] = calcShelfPenetration(runs)
    ticstart = tic;

    runs.read_eddsurf;
    runs.read_rhosurf;

    % look at values 'factor' × eddy radii behind center
    % too close to eddy and you miss a lot.
    factor = 3;

    % (along-shelf, cross-shelf, time)
    if runs.bathy.axis == 'y'
        eddsurf = runs.eddsurf;
        rhosurf = runs.rhosurf - runs.rhosurf(1,1,1);
        imx = runs.eddy.imx;
        yvec = runs.rgrid.y_rho(1:runs.bathy.isb, 1);
    else
        eddsurf = permute(runs.eddsurf, [2 1 3]);
        rhosurf = permute(runs.rhosurf, [2 1 3]);
        imx = runs.eddy.imy;
        yvec = runs.rgrid.x_rho(1, 1:runs.bathy.isb);
        factor = factor * -1 * runs.sgntamp;
    end

    if runs.bathy.loc == 'h'
        error('Not adapted for northern coast');
    end

    % number of points in one radius
    nr = ceil(runs.eddy.rhovor.dia(1)/2/1000);

    for tt = 1:runs.eddy.tend
        ix = min(imx(tt) + factor*nr, size(eddsurf,1));
        onshelf.edd.profile(:,tt) = squeeze(eddsurf(ix, 1:runs.bathy.isb, tt));
        onshelf.rho.profile(:,tt) = squeeze(rhosurf(ix, 1:runs.bathy.isb, tt));
    end

    onshelf.edd = calcdiags(onshelf.edd, runs.eddy_thresh, yvec, runs.bathy.isb);
    onshelf.rho = calcdiags(onshelf.rho, 0.4*max(abs(onshelf.rho.profile(:))), yvec, runs.bathy.isb);

    onshelf.rho.minrho = min(onshelf.rho.profile, [], 1);

    hash = githash([mfilename('fullpath') '.m']);
    onshelf.hash = hash;
    onshelf.factor = factor;
    onshelf.nr = nr;

    runs.onshelf = onshelf;
    save([runs.dir '/onshelf.mat'], 'onshelf');
    toc(ticstart);
end

function [in] = calcdiags(in, thresh, yvec, isb)

    profile = abs(in.profile) > abs(thresh); % (cross-shelf, time)
    in.thresh = thresh;

    % calculate envelope as width of inflow
    in.env = fillnan(max(bsxfun(@times, profile,abs(yvec - yvec(isb))), [], 1), 0);
end