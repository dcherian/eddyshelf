% Given a time index, show surface map.
% Then, allow choosing a point and plot up z-profiles
function [handles,xx,yy] = secondary_vortices(runs, tindex, n, opt)

    if ~exist('n', 'var'), n = 1; end
    if ~exist('opt', 'var'), opt = []; end

    figure; insertAnnotation([runs.name '.secondary_vortices']);
    handles.hax(1) = subplot(3,2,[1 2]);
    hfield = runs.animate_field('csdye', handles.hax(1), tindex, 1, opt);
    xlim([-1 1]*120 + runs.eddy.mx(tindex)/1000);
    ylim([-1 1]*120 + runs.eddy.my(tindex)/1000);
    hfield.hcb.LimitsMode = 'auto';
    hfield.hfield.CData = hfield.hfield.CData - runs.bathy.xsb/1000;
    caxis(hfield.hcb.Limits - runs.bathy.xsb/1000);
    center_colorbar(hfield.hcb);
    hfield.hcb.Position(1) = hfield.hcb.Position(1) + 0.02;
    index = find(hfield.hcb.Ticks == 0);
    hfield.hcb.TickLabels{index} = '0 - Shelfbreak';

    if size(n,2) ~= 2
        [xx,yy] = ginput(n);
    else
        xx = n(:,1);
        yy = n(:,2);

        n = size(n, 1);
    end
    colors = cbrewer('qual', 'Dark2', n);

    pvback = runs.params.phys.f0 * runs.params.phys.N2/ ...
             runs.params.phys.g;

    tic;
    disp('Reading data');
    for nn = 1:n
        x = xx(nn); y = yy(nn);
        plot(x,y,'x','Color',colors(nn,:),'MarkerSize', 10);

        x = x*1e3; y = y*1e3;
        ix = find_approx(runs.rgrid.x_rho(1,:), x, 1);
        iy = find_approx(runs.rgrid.y_rho(:,1), y, 1);

        volr = {'x' ix ix; 'y' iy iy};
        zr(:,nn) = runs.rgrid.z_r(:,iy,ix);

        csdye(:,nn) = dc_roms_read_data(runs, runs.csdname, tindex, volr)/1000;
        rho(:,nn) = dc_roms_read_data(runs, 'rho', tindex, volr);
        rback(:,nn) = dc_roms_read_data(runs, 'rho', 1, volr);
        try
            [pv(:,nn),~,~,zpv(:,nn),~] = dc_roms_read_data(runs, 'pv', tindex, volr);
            rv(:,nn) = dc_roms_read_data(runs, 'rv', tindex, volr);
            nopv = 0;
        catch ME
            nopv = 1;
        end
    end
    toc;

    backup = get(groot, 'DefaultAxesColorOrder');
    set(groot, 'DefaultAxesColorOrder', colors);

    handles.hax(2) = subplot(323);
    plot(csdye - runs.bathy.xsb/1000, zr);
    linex(0);
    common(runs, tindex);
    xlabel('Cross-shelf dye - X_{sb} (km)');

    handles.hax(3) = subplot(324);
    plot(rho-rback, zr);
    linex(0);
    common(runs, tindex);
    handles.hax(3).XTickLabelRotation = 45;
    xlabel('\Delta\rho (kg/m^3)')
    reduceSubplotHorizontalSpace(handles.hax(2:3));

    if ~nopv
        handles.hax(4) = subplot(325);
        plot(pv-pvback, zpv);
        linex(0);
        common(runs, tindex);
        xlabel('\Delta PV');

        handles.hax(5) = subplot(326)
        plot(rv/runs.params.phys.f0, zpv);
        linex(0); common(runs, tindex);
        xlabel('rv/f_0');

        reduceSubplotHorizontalSpace(handles.hax(4:5));
    end

    linkaxes(handles.hax(2:end), 'y');
    ylim([-2 0] * runs.eddy.Lgauss(tindex));
    set(groot, 'DefaultAxesColorOrder', backup);
end

function [] = common(obj, tindex)

    ylabel('Z (m)');
    liney(-1 * obj.bathy.hsb, 'H_{sb}');
    liney(-1 * obj.eddy.Lgauss(tindex), 'L_z^{eddy}');
    pbaspect([1.65 1 1]); %axis tight; axis square;
    beautify;
end

function reduceSubplotHorizontalSpace(ax)

    width = ax(1).Position(3);
    dx = ax(2).Position(1) - (ax(1).Position(1) + width);
    ax(1).Position(1) = ax(1).Position(1) + dx*3/4;
    ax(2).Position(1) = ax(2).Position(1) - dx*3/4;

end