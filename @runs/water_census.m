
% water mass census in full domain
function [] = water_census(runs)

    ticstart = tic;
    % if dye_04 > thresh then it is "eddy water"
    % else i classify it as "mixed"
    eddye_thresh = runs.eddy_thresh;

    % check classified vol against total volume
    debug = 0;

    % my region boundaries are based on location of shelfbreak and
    % slopebreak. Let's make it easy.
    xsl = runs.bathy.xsl;
    isl = runs.bathy.isl;
    xsb = runs.bathy.xsb;
    isb = runs.bathy.isb;

    slab = 40; % read 10 at a time

    sz4dfull = [fliplr(size(runs.rgrid.z_r)) slab];
    sz4dsp = [prod(sz4dfull(1:3)) slab];
    sz3dfull = sz4dfull(1:3);
    sz3dsp = [sz4dsp(1) 1];

    % define cross-shore grid co-ordinate
    if runs.bathy.axis == 'y'
        cs = repmat(runs.rgrid.yr,[1 1 runs.rgrid.N]);
    else
        cs = repmat(runs.rgrid.xr,[1 1 runs.rgrid.N]);
    end

    % define regions
    % deep region
    regdp = sparse(reshape(cs > xsl, sz3dsp));
    regsl = sparse(reshape(cs <= xsl & cs >= xsb, sz3dsp));
    regsh = sparse(reshape(cs < xsb, sz3dsp));

    % include sponge filtering in dV
    % i.e., set dV=0 in sponge region
    dV = reshape(bsxfun(@times, runs.rgrid.dV, ~runs.sponge), sz3dsp);

    % not sure if dz is needed
    %dz = dV ./ runs.rgrid.dx ./ runs.rgrid.dy;
    xsp = reshape(repmat(runs.rgrid.xr,[1 1 runs.rgrid.N]), sz3dsp);
    ysp = reshape(repmat(runs.rgrid.yr,[1 1 runs.rgrid.N]), sz3dsp);
    zsp = reshape(permute(runs.rgrid.z_r,[3 2 1]), sz3dsp);

    ntime = length(runs.time);

    for tt=1:slab:length(runs.time)
        tend = tt + slab -1;
        if tend > length(runs.time)
            tend = length(runs.time);
            sz4dfull(end) = tend-tt+1;
            sz4dsp(end) = tend-tt+1;
        end
        csdye = dc_roms_read_data(runs.dir,runs.csdname,[tt ...
                            tend],{},[],runs.rgrid, 'avg', 'single');
        eddye = dc_roms_read_data(runs.dir,runs.eddname,[tt ...
                            tend],{},[],runs.rgrid, 'avg', 'single');
        % define water masses
        % offshore water
        maskoff = sparse(reshape(csdye > xsl, sz4dsp));
        % slope water
        masksl  = sparse(reshape(csdye <= xsl & csdye >= xsb, sz4dsp));
        % shelf water
        masksh  = sparse(reshape(csdye < xsb, sz4dsp));
        % eddy water
        masked  = sparse(reshape(eddye > eddye_thresh, sz4dsp));
        % "mixed water"
        maskmix = sparse(reshape(eddye <= eddye_thresh & eddye > 0.01,sz4dsp));

        % shift in water parcels
        %dcsmask = reshape(bsxfun(@minus, csdye, cs), sz4dsp)/ ...
        %          1000;
        %dcsmask = dcsmask .* (abs(dcsmask)<0.5);

        %cssh =  (dcsmask .* full(masksh)) > 0;
        %masksl(cssh) = 1;


        % the eddy's velocity field mixes up csdye and makes it look
        % like slope water?
        % in any case i want all 5 to add up to total volume, so let's
        % remove the volume that's in the eddy from csdye.
        maskoff = maskoff & ~(masked | maskmix);
        masksl  = masksl  & ~(masked | maskmix);
        masksh  = masksh  & ~(masked | maskmix);

        % now census in each region
        % first deep water region
        runs.water.off.deep(tt:tend) = full(sum( ...
            bsxfun(@times, maskoff, regdp .* dV),1));
        runs.water.sl.deep(tt:tend) = full(sum( ...
            bsxfun(@times, masksl, regdp .* dV),1));
        runs.water.sh.deep(tt:tend) = full(sum( ...
            bsxfun(@times, masksh, regdp .* dV),1));
        runs.water.edd.deep(tt:tend) = full(sum( ...
            bsxfun(@times, masked, regdp .* dV),1));
        runs.water.mix.deep(tt:tend) = full(sum( ...
            bsxfun(@times, maskmix, regdp .* dV),1));

        % now slope region
        runs.water.off.slope(tt:tend) = full(sum( ...
            bsxfun(@times, maskoff, regsl .* dV),1));
        runs.water.sl.slope(tt:tend) = full(sum( ...
            bsxfun(@times, masksl, regsl .* dV),1));
        runs.water.sh.slope(tt:tend) = full(sum( ...
            bsxfun(@times, masksh, regsl .* dV),1));
        runs.water.edd.slope(tt:tend) = full(sum( ...
            bsxfun(@times, masked, regsl .* dV),1));
        runs.water.mix.slope(tt:tend) = full(sum( ...
            bsxfun(@times, maskmix, regsl .* dV),1));

        % now shelf region
        runs.water.off.shelf(tt:tend) = full(sum( ...
            bsxfun(@times, maskoff, regsh .* dV),1));
        runs.water.sl.shelf(tt:tend) = full(sum( ...
            bsxfun(@times, masksl, regsh .* dV),1));
        runs.water.sh.shelf(tt:tend) = full(sum( ...
            bsxfun(@times, masksh, regsh .* dV),1));
        runs.water.edd.shelf(tt:tend) = full(sum( ...
            bsxfun(@times, masked, regsh .* dV),1));
        runs.water.mix.shelf(tt:tend) = full(sum( ...
            bsxfun(@times, maskmix, regsh .* dV),1));

        % statistics of shelf water on the slope
        runs.water.sh.xslope(tt:tend) = full(sum( ...
            bsxfun(@times, masksh, regsl .* xsp .*dV),1)) ./ ...
            (runs.water.sh.slope(tt:tend));
        runs.water.sh.yslope(tt:tend) = full(sum( ...
            bsxfun(@times, masksh, regsl .* ysp .*dV),1)) ./ ...
            (runs.water.sh.slope(tt:tend));
        runs.water.sh.zslope(tt:tend) = full(sum( ...
            bsxfun(@times, masksh, regsl .* zsp .*dV),1)) ./ ...
            (runs.water.sh.slope(tt:tend));

        % statistics of eddy/mix water on shelf.
        % - x location, y-location, z-location
        runs.water.eddmix.xshelf(tt:tend) = full(sum( ...
            bsxfun(@times, (maskmix | masked), regsh .* xsp .*dV),1)) ./ ...
            (runs.water.edd.shelf(tt:tend) + runs.water.mix.shelf(tt:tend));

        runs.water.eddmix.yshelf(tt:tend) = full(sum( ...
            bsxfun(@times, (maskmix | masked), regsh .* ysp .*dV),1)) ./ ...
            (runs.water.edd.shelf(tt:tend) + runs.water.mix.shelf(tt:tend));

        runs.water.eddmix.zshelf(tt:tend) = full(sum( ...
            bsxfun(@times, (maskmix | masked), regsh .* zsp .*dV),1)) ./ ...
            (runs.water.edd.shelf(tt:tend) + runs.water.mix.shelf(tt:tend));

        % how uniform is the "plume" in the vertical - i.e.,
        % baroclinicity - RMS (tracer) / RMS (depth avg tracer)

    end
    toc(ticstart);

    time = runs.time/86400;
    runs.water.comment = [''];
    water = runs.water;

    water.totvol = sum(dV(:));
    water.shvol = sum(dV(:) .* full(regsh));
    water.slvol = sum(dV(:) .* full(regsl));
    water.dpvol = sum(dV(:) .* full(regdp));

    water.hash = githash;

    save([runs.dir '/watermass.mat'], 'water');

    %%
    % calculate total classified volume
    if debug
        masses = fieldnames(runs.water);
        classvol = zeros(size(runs.water.mix.deep));
        for ii=1:length(masses)
            if strcmpi(masses{ii},'eddmix'), continue; end
            try
                regions = fieldnames(runs.water.(masses{ii}));
                for jj=1:length(regions)
                    if regexp(regions{jj},'^[xyz]'), continue; end
                    classvol = classvol + runs.water.(masses{ii}).(regions{jj});
                end
            catch ME
            end
        end
        water.classvol = classvol;
        figure;
        plot(water.totvol - water.classvol);
    end
end
