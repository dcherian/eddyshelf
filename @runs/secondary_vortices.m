% Given a time index, show surface map.
% Then, allow choosing a point and plot up z-profiles
function [] = secondary_vortices(runs, tindex, n)

    if ~exist('n', 'var'), n = 1; end

    figure; insertAnnotation([runs.name '.secondary_vortices']);
    hax = subplot(3,2,1);
    runs.animate_field('csdye', hax, tindex, 1);
    xlim([-100 100] + runs.eddy.mx(tindex)/1000);
    ylim([-100 100] + runs.eddy.my(tindex)/1000);

    [xx,yy] = ginput(n);
    colors = cbrewer('qual', 'Dark2', n);

    pvback = runs.params.phys.f0 * runs.params.phys.N2/ ...
             runs.params.phys.g;

    tic;
    disp('Reading data');
    for nn = 1:n
        x = xx(nn); y = yy(nn);
        plot(x,y,'x','Color',colors(nn,:),'MarkerSize', 16);

        x = x*1e3; y = y*1e3;
        ix = find_approx(runs.rgrid.x_rho(1,:), x, 1);
        iy = find_approx(runs.rgrid.y_rho(:,1), y, 1);

        volr = {'x' ix ix; 'y' iy iy};
        zr(:,nn) = runs.rgrid.z_r(:,iy,ix);

        csdye(:,nn) = dc_roms_read_data(runs, runs.csdname, tindex, volr)/1000;
        rho(:,nn) = dc_roms_read_data(runs, 'rho', tindex, volr);
        rback(:,nn) = dc_roms_read_data(runs, 'rho', 1, volr);
        [pv(:,nn),~,~,zpv(:,nn),~] = dc_roms_read_data(runs, 'pv', tindex, volr);
        rv(:,nn) = dc_roms_read_data(runs, 'rv', tindex, volr);
    end
    toc;

    backup - get(groot, 'DefaultAxesColorOrder);
    set(groot, 'DefaultAxesColorOrder', colors);
    ax(1) = subplot(323);
    plot(csdye - runs.bathy.xsb/1000, zr);
    linex(0);
    common(runs, tindex);
    xlabel('csdye');

    ax(2) = subplot(324);
    plot(rho-rback, zr);
    linex(0);
    common(runs, tindex);
    xlabel('\Delta\rho')

    ax(3) = subplot(325);
    plot(pv-pvback, zpv);
    linex(0);
    common(runs, tindex);
    xlabel('\Delta PV');

    ax(4) = subplot(326)
    plot(rv/runs.params.phys.f0, zpv);
    linex(0); common(runs, tindex);
    xlabel('rv/f_0');

    linkaxes(ax, 'y');
    ylim([-2 0] * runs.eddy.Lgauss(tindex));
    set(groot, 'DefaultAxesColorOrder', backup);
end

function [] = common(obj, tindex)

    liney(-1 * obj.bathy.hsb);
    liney(-1 * obj.eddy.Lgauss(tindex));
    axis tight; axis square;
    beautify;
end