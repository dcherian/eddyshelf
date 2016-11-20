% Averages 'varname' at 'axname' = 'ix'
% Average is taken over time interval returned by flux_tindices.
%   mask == 1, indicate boundary of shelf/eddy water
%      [varmean,maskmean,xivec, varmaskedmean] = avgProfile(runs, varname, axname, ix, mask)
function [varmean,maskmean,xivec, varmaskedmean] = ...
        avgProfile(runs, varname, axname, ix, mask, reference)

    if ~exist('mask', 'var'), mask = 0; end
    if ~exist('reference', 'var'), reference = 'center'; end

    varname = runs.process_varname(varname);

    [xres, yres] = runs.locate_resistance;

    if strcmpi(ix, 'sb')
        ix = runs.bathy.isb;
    end

    if axname == 'x'
        if ~exist('ix', 'var') | isempty(ix) | strcmpi(ix, 'cen')
            ix = runs.eddy.imx(tindex);
        end
        if strcmpi(ix, 'res')
            ix = num2str(xres);
        end
    else
        if ~exist('ix', 'var') | isempty(ix) | strcmpi(ix, 'cen')
            ix = runs.eddy.imy(tindex);
        end
        if strcmpi(ix, 'res')
            ix = num2str(yres);
        end
    end

    if ~strcmpi(varname, 'zeta') ...
            & ~strcmpi(varname, 'ubar') & ~strcmpi(varname, 'vbar')
        vol = {axname ix ix; 'z' runs.rgrid.N runs.rgrid.N};
    else
        vol = {axname ix ix};
    end
    volmask = {axname ix ix; 'z' runs.rgrid.N runs.rgrid.N};

    [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:,1,1));

    [var, xax, yax]  = dc_roms_read_data(runs, varname, [], vol);
    % restrict to unique time indices
    [tvec, uind] = unique(runs.eddy.t);
    var = var(:,uind);

    if mask
        % Using csdye > bathy.xsl works a *lot* better than using
        % interpolated eddye > runs.eddy_thresh. Results are not
        % sensitive to bathy.xsl or (bathy.xsb + bathy.xsl)/2.
        % bathy.xsb is not good because at the western edge some dye
        % moves inshore, screwing up mask detection (first change from
        % 0 to 1).
        maskvar = dc_roms_read_data(runs, runs.csdname, [], volmask);
        maskvar = maskvar(:,uind);
    end

    if axname == 'x'
        % for x = ix, plot against y
        xvec = yax(1,:);
        mx = runs.bathy.xsb * ones(size(runs.eddy.my)); runs.eddy.my;
    else
        % for y = iy, plot against x
        xvec = xax(:,1);
        mx = runs.eddy.mx(uind);
    end

    xivec = [-200:200]*1000;
    for tt=1:length(uind) % only as long as I can track the eddy.
        if ~strcmpi(reference, 'center')
            % find location of eddy-shelf/slope water boundary
            ref = find(maskvar(:,tt) > runs.bathy.xsl, 1, 'first');
            if isempty(ref)
                mx(tt) = nan;
                continue;
            else
                mx(tt) = xvec(ref);
            end
        end

        % interpolate to common reference frame.
        vari(:,tt) = interp1(xvec-mx(tt), var(:,tt), xivec);
        if mask
            maski(:,tt) = interp1(xvec-mx(tt), maskvar(:,tt), xivec);
            varmaski(:,tt) = interp1(xvec-mx(tt), ...
                                     var(:,tt).*(maskvar(:,tt) < runs.bathy.xsl), xivec);
        end
    end

    % only average when eddy has penetrated shelf.
    % start = find(~isnan(mx) == 1, 1, 'first') % this makes start=1 always
    varmean = nanmean(vari(:,start:stop), 2);
    if mask
        maskmean = nanmean(maski(:,start:stop), 2) > runs.bathy.xsl;
        varmaskedmean = nanmean(varmaski(:,start:stop), 2);
    else
        maskmean = ones(size(varmean));
    end
end