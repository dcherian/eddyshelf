% calculate and save average along-shelf supply velocity (using csdye for mask) and SSH;
function [] = avgSupplyJet(runs, debug)

    ticstart = tic;

    if ~exist('debug', 'var'), debug = 0; end

    assert(runs.bathy.axis == 'y', 'Error: EW isobaths only!');

    [xx,yy,tt] = runs.locate_resistance;
    [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:, 1, 1));
    tindices = [start stop];

    xivec = -100:runs.rgrid.dx/1000:0;
    xl = length(xivec);
    N = runs.rgrid.N;
    isb = runs.bathy.isb;

    xloc = xx + 2 * runs.params.eddy.dia;
    xind = find_approx(runs.rgrid.x_rho(1,:), xloc, 1);
    assert(xind < runs.spng.sx2, 'Error: Location within sponge!');
    yvec = runs.rgrid.y_rho(2:isb, 1);
    zmat = runs.rgrid.z_r(:,2:isb,xind)';
    ymat = repmat(yvec, [1 N]);

    % ignore wall ghost point
    volr = {'x' xind xind;
            'y' 2 isb};

    volv = {'x' xind-1 xind;
            'y' 2 isb};

    % vars are (y,z,t)
    asvel = squeeze(avg1(dc_roms_read_data(runs, runs.asvelname, tindices, volv),1));
    csdye = dc_roms_read_data(runs, runs.csdname, tindices, volr);
    zeta = dc_roms_read_data(runs, 'zeta', tindices, volr);

    % restrict to unique time indices
    [tvec, uind] = unique(runs.eddy.t(start:stop));
    asvel = asvel(:,:,uind);
    csdye = csdye(:,:,uind);
    zeta = zeta(:,uind);

    % look for shelf water all throughout water column
    %csdvint = squeeze(sum(csdye < runs.bathy.xsl, 2)) >= runs.rgrid.N;
    csdvint = squeeze(csdye(:,end,:) < runs.bathy.xsl);

    % the fit doesn't seem helpful
    % for tt=1:size(zeta,2)
    %     zetamasked = cut_nan(zeta(:,tt) .* csdvint(:,tt));
    %     % Eddy's SSH is gaussian, so presumably that's how it's decaying
    %     [z0,Xzeta(tt),z1] = gauss_fit(yvec(1:length(zetamasked)), ...
    %                                   zetamasked - min(zetamasked), 1);
    % end
    zetamean = nanmean(zeta .* csdvint, 2);
    zetamean = zetamean - mean(zetamean);
    asvshmean = mean(asvel .* (csdye <= runs.bathy.xsb), 3);
    asvslmean = mean(asvel .* (csdye <= runs.bathy.xsl), 3);
    csdmean = mean(csdye, 3);
    asvmean = mean(asvel, 3);

    for ii=1:size(asvmean,1)
        asvshint(ii) = trapz(zmat(ii,:), asvshmean(ii,:));
        asvslint(ii) = trapz(zmat(ii,:), asvslmean(ii,:));
        asvint(ii) = trapz(zmat(ii,:), asvmean(ii,:));
    end

    % fit time-averaged, dye-masked SSH
    [zeta0, Xzeta, X0, zeta1,zetaconf, zfitobj] = tanh_fit(yvec, zetamean, debug);

    if debug
        title(['zeta | ' runs.name]);
    end

    asvavg   = asvint./runs.bathy.h(1,2:isb);
    asvshavg = asvshint./runs.bathy.h(1,2:isb);
    asvslavg = asvslint./runs.bathy.h(1,2:isb);

    % fit vertically averaged, time-averaged, dye-masked along-shelf velocity
    [~,imin] = min(asvshint./runs.bathy.h(1,2:isb));
    exitflag = 0;
    i0 = 1;
    while ~exitflag
        [v0,X,x0,v1,exitflag,conf] = gauss_fit(yvec(i0:imin), asvshavg(i0:imin), debug);
        if ~exitflag
            i0 = i0 + 1;
        end
        % this condition is bad. I can make any value I want.
        %if (yvec(imin) + X) > runs.bathy.xsb*1.5
        %    exitflag = 0;
        %end
    end

    % fit vertically averaged, time-averaged, *unmasked* along-shelf velocity
    [v0,Xfull,x0,v1,exitflag,confFull] = gauss_fit(yvec, asvavg, debug);
    supply.IntersectIndex = find(abs(asvavg - asvshavg) > 0.01, 1, 'first');
    supply.IntersectLocation = yvec(supply.IntersectIndex);
    supply.IntersectScale = runs.bathy.xsb - supply.IntersectLocation;

    [v0,Xsl,x0,v1,exitflag,confsl] = gauss_fit(yvec, asvavg, debug);

    if debug, title(runs.name); end
    supply.shelf.vmean = asvshmean;
    supply.shelf.vint = asvshint;
    supply.shelf.vavg = asvshavg;
    supply.shelf.xscale = X;
    supply.shelf.xmin = yvec(imin);
    supply.shelf.imin = imin;
    supply.shelf.conf = conf(:,2);
    supply.shelf.comment = 'shelf-water only';

    supply.shsl.vmean = asvslmean;
    supply.shsl.vint = asvslint;
    supply.shsl.vavg = asvslavg;
    supply.shsl.xscale = Xsl;
    supply.shsl.conf = confsl(:,2);
    supply.shsl.comment = 'shelf + slope water only';

    supply.zeta.zetamean = zetamean;
    supply.zeta.xscale = Xzeta;
    supply.zeta.conf = zetaconf(:,2);
    supply.zeta.fitobj = zfitobj;

    supply.vmean = asvmean;
    supply.vint = asvint;
    supply.vavg = asvavg;
    supply.xscale = Xfull;
    supply.conf = confFull(:,2);
    supply.csdmean = csdmean;

    supply.hash = githash([mfilename('fullpath') '.m']);
    supply.tindices = tindices;
    supply.ymat = ymat;
    supply.zmat = zmat;
    supply.x = xloc;
    supply.ix = xind;
    supply.i0 = i0;
    supply.comment = ['Along-shelf velocity averaged (flux_tindices for sb) ' ...
                      'as function of (y,z) ' ...
                      '| (v*mean, v*int, v*avg) = (mean, depth integrated, depth averaged) | ' ...
                      'i0 ≠ 0 means I had to remove points near coast | ' ...
                      'Intersect(Index,Location,Scale) = point where shelf.vint and ' ...
                      ' vint intersect: this is the last point where the jet is ' ...
                      ' all shelf water. ' ...
                      ' IntersectScale is essentially a length scale for penetration ' ...
                      'of slope-eddy water on the shelf.'];

    runs.supply = supply;

    save([runs.dir '/supplyjet.mat'], 'supply');
    toc(ticstart);
end