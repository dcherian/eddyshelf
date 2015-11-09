% Averages 'varname' at 'axname' = 'ix' after applying 'mask'
%   mask == 1, indicate boundary of shelf/eddy water
%      [varmean,maskmean,xivec] = avgProfile(runs, varname, axname, ix, mask)
function [varmean,maskmean,xivec] = avgProfile(runs, varname, axname, ix, mask)

    if ~exist('mask', 'var'), mask = 0; end

    varname = runs.process_varname(varname);

    if strcmpi(ix, 'sb')
        ix = runs.bathy.isb;
    end

    if axname == 'x'
        if ~exist('ix', 'var') | isempty(ix) | strcmpi(ix, 'cen')
            ix = runs.eddy.imx(tindex);
        end
    else
        if ~exist('ix', 'var') | isempty(ix) | strcmpi(ix, 'cen')
            ix = runs.eddy.imy(tindex);
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

    if mask
        % Using csdye > bathy.xsl works a *lot* better than using
        % interpolated eddye > runs.eddy_thresh. Results are not
        % sensitive to bathy.xsl or (bathy.xsb + bathy.xsl)/2.
        % bathy.xsb is not good because at the western edge some dye
        % moves inshore, screwing up mask detection (first change from
        % 0 to 1).
        maskvar = dc_roms_read_data(runs, runs.csdname, [], volmask);
    end

    if axname == 'x'
        % for x = ix, plot against y
        xvec = yax(1,:);
        mx = runs.eddy.my;
    else
        % for y = iy, plot against x
        xvec = xax(:,1);
        mx = runs.eddy.mx;
    end

    xivec = [-200:200]*1000;
    for tt=1:size(var,2)
        vari(:,tt) = interp1(xvec-mx(tt), var(:,tt), xivec);
        if mask
            maski(:,tt) = interp1(xvec-mx(tt), maskvar(:,tt), xivec);
        end
    end

    varmean = nanmean(vari(:,start:stop), 2);
    if mask
        maskmean = nanmean(maski(:,start:stop), 2) > runs.bathy.xsl;
    else
        maskmean = ones(size(varmean));
    end
end