% calculate upwelling in eddy
function [] = eddy_upwelling(runs)

% two possibilities here - use eddye as a
% 1 - use mask
% 2 - use vormask
    xr = runs.eddy.xr;
    yr = runs.eddy.yr;

    use_sshmask = 0;

    ixmin = vecfind(runs.rgrid.xr(:,1),runs.eddy.vor.we);
    ixmax = vecfind(runs.rgrid.xr(:,1),runs.eddy.vor.ee);
    iymin = vecfind(runs.rgrid.yr(1,:),runs.eddy.vor.se);
    iymax = vecfind(runs.rgrid.yr(1,:),runs.eddy.vor.ne);

    ixm = min(ixmin); ixM = max(ixmax);
    iym = min(iymin); iyM = max(iymax);

    volume = {'x' ixm ixM; 'y' iym iyM};
    zdye = dc_roms_read_data(runs.dir, runs.zdname, [], volume, [], runs.rgrid);

    % from animate_pt output it looks like vormask tracks the edge of
    % eddye contour pretty well. I don't get the stuff that spreads
    % along shelf but get the eddy pretty well.
    try
        eddye = dc_roms_read_data(runs.dir,runs.eddname, [], volume, [], runs.rgrid);
    catch ME
        warning('no eddy dye (dye_04) found');
        return;
    end

    sz4dfull = size(zdye);
    sz4dsp = [prod(sz4dfull(1:3)) sz4dfull(end)];
    sz3dsp = [sz4dsp(1) 1];

    % make my mask matrices 4d and sparse
    eddye = sparse(reshape(eddye > runs.eddy_thresh, sz4dsp));
    vormask = sparse(reshape(repmat( ...
        permute(runs.eddy.vor.mask(ixm-1:ixM-1,iym-1:iyM-1,:),[1 2 4 3]) ...
        , [1 1 runs.rgrid.N 1]), sz4dsp));

    % this is the combined mask matrix
    evormask = eddye .* vormask;
    clear eddye vormask

    disp('cleared memory');

    dV = reshape(runs.rgrid.dV(ixm:ixM,iym:iyM,:),sz3dsp);

    zdye = reshape(zdye,sz4dsp);

    % first with vormask
    zdyevor = zdye .* evormask;
    zedd = bsxfun(@times, evormask, reshape(permute( ...
        runs.rgrid.z_r(:,iym:iyM,ixm:ixM),[3 2 1]), sz3dsp));
    runs.eddy.vor.vol = full(squeeze(sum(bsxfun(@times, evormask, ...
                                                dV),1))');
    runs.eddy.vor.zdcen = runs.domain_integratesp(zdyevor, dV)' ...
        ./ runs.eddy.vor.vol;
    runs.eddy.vor.zcen = runs.domain_integratesp(zedd,dV)' ...
        ./ runs.eddy.vor.vol;

    % then with ssh mask - though really vormask is what i'm looking
    % for
    %{
    if use_sshmask
        mask = sparse(reshape(repmat( ...
        permute(runs.eddy.mask(ixm-1:ixM-1,iym-1:iyM-1,:),[1 2 4 3]) ...
        , [1 1 runs.rgrid.N 1]), sz4dsp));
        zdyessh = zdye .* mask;
        runs.eddy.vol = full(squeeze(sum(bsxfun(@times,mask,dV) ,1))');
        runs.eddy.zdcen = runs.domain_integratesp(zdyessh, reshape(dV,sz3dsp))' ...
        ./ runs.eddy.vol;
        end
        %}

        eddy = runs.eddy;
        eddy.hash = githash;
        save([runs.dir '/eddytrack.mat'], 'eddy');
    end
