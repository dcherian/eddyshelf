% calculate and save average along-shelf supply velocity (using csdye for mask) and SSH;
function [] = avgSupplyJet(runs)

    ticstart = tic;

    assert(runs.bathy.axis == 'y', 'Error: EW isobaths only!');

    [xx,yy,tt] = runs.locate_resistance;
    [start,stop] = runs.flux_tindices(runs.csflux.off.slope(:, 1, 1));
    tindices = [start ceil((start+stop)/2)];

    xivec = -100:runs.rgrid.dx/1000:0;
    xl = length(xivec);
    N = runs.rgrid.N;
    isb = runs.bathy.isb;

    xloc = xx + 30e3;
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

    % look for shelf water all throughout water column
    csdvint = fillnan(double( ...
        squeeze(sum(csdye <= runs.bathy.xsb, 2)) == runs.rgrid.N), 0);

    % the fit doesn't seem helpful
    % for tt=1:size(zeta,2)
    %     zetamasked = cut_nan(zeta(:,tt) .* csdvint(:,tt));
    %     % Eddy's SSH is gaussian, so presumably that's how it's decaying
    %     [z0,Xzeta(tt),z1] = gauss_fit(yvec(1:length(zetamasked)), ...
    %                                   zetamasked - min(zetamasked), 1);
    % end

    zetamean = nanmean(zeta .* csdvint, 2);
    asvmean = mean(asvel .* (asvel < 0) .* (csdye <= runs.bathy.xsl), 3);
    csdmean = mean(csdye, 3);

    for ii=1:size(asvmean,1)
        asvint(ii) = trapz(zmat(ii,:), asvmean(ii,:));
    end

    % fit time-averaged, dye-masked SSH
    [~,imax] = nanmax(zetamean);
    [zeta0, Xzeta, zeta1] = gauss_fit(yvec(1:imax), ...
                                      zetamean(1:imax) - min(zetamean(1:imax)), 0);

    % fit vertically integrated, time-averaged, dye-masked along-shelf velocity
    [~,imin] = min(asvint./runs.bathy.h(1,2:isb));
    exitflag = 0;
    i0 = 1;
    while ~exitflag
        [v0,X,v1,exitflag] = gauss_fit(yvec(i0:imin), asvint(i0:imin), 0);
        if ~exitflag
            i0 = i0 + 1;
        end
        % this condition is bad. I can make any value I want.
        %if (yvec(imin) + X) > runs.bathy.xsb*1.5
        %    exitflag = 0;
        %end
    end
    title(runs.name);

    supply.hash = githash([mfilename('fullpath') '.m']);
    supply.tindices = tindices;
    supply.asvmean = asvmean;
    supply.asvint = asvint;
    supply.xscale = X;
    supply.xscaleZeta = Xzeta;
    supply.xmin = yvec(imin);
    supply.imin = imin;
    supply.csdmean = csdmean;
    supply.ymat = ymat;
    supply.zmat = zmat;
    supply.x = xloc;
    supply.ix = xind;
    supply.i0 = i0;
    supply.comment = ['Along-shelf velocity averaged (flux_tindices for sb) ' ...
                      'as function of (y,z) ' ...
                      '| (asvmean, asvint) = (mean, depth integrated) | ' ...
                      'i0 ≠ 0 means I had to remove points near coast'];

    runs.supply = supply;

    save([runs.dir '/supplyjet.mat'], 'supply');
    toc(ticstart);
end