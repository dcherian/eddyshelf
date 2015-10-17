% plots average velocity seen by shelf-slope water at a given isobath
function [] = avgStreamerVelSection(runs, isobath, source)

    debug = 0;
    if ~exist('source', 'var')
        warning('Using source = isobath');
        source = isobath;
    end

    if isobath > length(runs.csflux.ix)
        error(['Choose a lower isobath number. Max = ' ...
               length(runs.csflux.ix)])
    end

    if source > isobath
        error('Source is offshore of isobath.');
    end

    [start,stop] = runs.flux_tindices(runs.csflux.west.slope(:, isobath, source));
    tindices = [start stop];
    ix = runs.csflux.ix(isobath); % WORKAROUND FOR CSFLUXES BUG
    volr = {runs.bathy.axis ix ix};
    volv = {runs.bathy.axis ix-1 ix};
    xivec = -200:runs.rgrid.dx/1000:200;
    zvec = runs.csflux.vertbins(:,isobath);

    if runs.bathy.axis == 'x'
        xvec = runs.rgrid.y_rho(2:end-1,1)/1000;
        mx = runs.eddy.my(start:stop)/1000;
        bathyax = 1;
    else
        xvec = runs.rgrid.x_rho(1,2:end-1)/1000;
        mx = runs.eddy.mx(start:stop)/1000;
        bathyax = 2;
    end

    tic; disp([runs.name ': Reading data...'])
    v = squeeze(avg1(dc_roms_read_data(runs, runs.csvelname, tindices, volv), bathyax));
    csdye = dc_roms_read_data(runs, runs.csdname, tindices, volr);
    toc;

    % interpolate onto common grid centered on eddy center
    % this is not working well!!!
    tic; disp([runs.name ': interpolating data...']);
    clear vi csdyei
    for tt=1:size(v,3)
        for zz=1:size(v,2)
            vi(:,zz,tt) = interp1(xvec-mx(tt), v(2:end-1,zz,tt), xivec);
            csdyei(:,zz,tt) = interp1(xvec-mx(tt), csdye(2:end-1,zz,tt), xivec);
        end
    end
    toc;

    vmeanwest = mean(bsxfun(@times, vi .* (csdyei < runs.csflux.x(source)), ...
                            xivec' < 0), 3);
    vmean = mean(vi .* (csdyei < runs.csflux.x(source)), 3);
    pmean = trapz(xivec, repnan(vmeanwest,0), 1); % mean profile
    pint = runs.csflux.west.slopewater.vertitrans(:,isobath,source);
    actualx = trapz(zvec, vmean, 2);

    hf = figure; maximize(); pause(1);
    insertAnnotation([runs.name '.avgStreamerVelSection']);
    ax1 = subplot(3,3,[4 5 7 8]);
    pcolorcen(xivec, zvec, vmean');
    xlabel('X - X_{eddy} (km)'); ylabel('Z (m)');
    liney(-runs.bathy.hsb);
    linex(0);
    title([runs.name ' | mean streamer velocity | y/R = ' num2str(runs.csflux.ndloc(isobath))]);
    hcb = center_colorbar;
    hcb.Position(1) = 0.5;

    ax2 = subplot(3,3,[1 2]);
    hx = plot(xivec, actualx);
    title('\int dz');
    hold on;
    linex(0); liney(0);
    ylabel('Flux (m^2/s)');

    ax3 = subplot(3,3,[6 9]);
    hz = plot(pmean, zvec);
    title('Offshore transport (\int dx)');
    xlabel('Flux (m^2/s)');

    if debug
        L = runs.eddy.rhovor.lmaj(1)/1000;
        R = runs.csflux.R/1000;

        a = 3; Ln = L/3;
        yoR = runs.csflux.ndloc(isobath); % y/R - used in csflux
        y0oL = R/L * (1 - yoR); % y0/L - used in derivation
        ideal = runs.streamer_ideal_profile(isobath);
        idealx = trapz(zvec, ideal) *  ...
                 diff(exp(-abs(xivec'/Ln).^a))./diff(xivec'/Ln) ...
                 * exp(-y0oL.^2);

        axes(ax2);
        hx.YData = hx.YData / max(actualx);
        plot(avg1(xivec), idealx./max(idealx), 'k-');
        ylabel('Flux / max flux');

        axes(ax3);
        hold on;
        hz.YData = hz.YData / max(pmean);
        plot(pint./max(pint), zvec);
        plot(ideal, zvec);
        legend('Mean', 'Integrated', 'Idealized', 'Location', 'SouthEast');
        xlabel('Flux / max flux');
    end

    linkaxes([ax1 ax2], 'x');
    linkaxes([ax1 ax3], 'y');
end