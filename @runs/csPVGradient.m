% cross-shelf PV gradient
function [] = csPVGradient(runs, fname, loc)
    if ~exist('fname', 'var') | isempty(fname)
        fname = 'pv-t180.nc';
    end

    [pv,xpv,ypv,zpv] = dc_roms_read_data([runs.dir '/' fname], 'pv');
    tpv = dc_roms_read_data([runs.dir '/' fname], 'ocean_time');
    depthavgPV = squeeze(sum(avg1(pv, 3) .* diff(zpv, 1, 3), 3)) ...
        ./ (zpv(:,:,end) - zpv(:,:,1));

    if runs.bathy.axis == 'y'
        % (along-shelf, cross-shelf, z)
        %dpv = diff_cgrid()... %diff(pv, 1, 2)./diff(ypv, 1, 2);

        % (along-shelf, cross-shelf)
        ddapv = diff(depthavgPV, 1, 2)./diff(ypv(:,:,end), 1, 2);

        asname = 'x';
        csvec = avg1(ypv(1,:,1),2);
        asvec = xpv(:,1,1);
    else
        % (along-shelf, cross-shelf, z)
        % dpv = diff_cgrid
        % dpv = permute(dpv, [2 1 3]);

        % (along-shelf, cross-shelf)
        ddapv = diff(depthavgPV, 1, 1)./diff(xpv(:,:,end), 1, 1)';

        asname = 'y';
        csvec = avg1(xpv(:,1,1),1);
        asvec = ypv(1,:,1)';
    end

    if ischar(loc)
        loc = find_approx(str2double(loc), asvec, 1);
    end

    tind = find_approx(runs.time, tpv(1), 1);

    [asvel,~,csvec_asvel] = dc_roms_read_data(runs, [runs.asvelname 'bar'], tind, ...
                                              {asname loc loc; 'z' runs.rgrid.N runs.rgrid.N});

    phys = runs.params.phys;
    beta = phys.beta;
    beta_sh = phys.f0/runs.bathy.hsb * runs.bathy.sl_shelf;
    ddapv = ddapv ./ (phys.N2/phys.g); % factor needed because I'm using Ertel PV.
                                       % This renormalizes so that dq/dy ~ β

    betashvec = beta_sh * (csvec < runs.bathy.xsb);

    nsmooth = 6;

    figure;
    insertAnnotation([runs.name '.csPVGradient']);
    ax(1) = subplot(211);
    plot((csvec - runs.bathy.xsb)/1000, smooth(ddapv(loc,:,end), nsmooth));
    hold on;
    plot((csvec - runs.bathy.xsb)/1000, smooth(ddapv(loc,:,end) - betashvec, nsmooth));
    legend('dq/dy', 'dq/dy + \beta_t', 'Location', 'NorthEast');
    xlabel('Cross-shelf dist. - X_{sb} (km)');
    title([runs.name ' | Depth averaged cross-shelf PV Gradient | t = ' num2str(tpv(1)/86400)]);
    text(0.70, 0.12, {[asname ' = ' num2str(asvec(loc)/1000, '%.0f') ' km']; ...
                      ['\beta = ' num2str(runs.params.phys.beta, '%.2e')]; ...
                      ['\beta_t = ' num2str(-beta_sh, '%.2e')]}, ...
         'HorizontalAlignment', 'left', 'Units', 'normalized');
    liney(0); linex(0);
    beautify;

    ax(2) = subplot(212);
    plot((csvec_asvel - runs.bathy.xsb)/1000, asvel);
    title('Depth-averaged along-shelf velocity');
    linex(0); liney(0);
    beautify;

    linkaxes(ax, 'x');
end
