% Given a time index, show surface map.
% Then, allow choosing a point and plot up z-profiles
function [handles,xx,yy] = secondary_vortices(runs, tindex, n, opt)

    if ~exist('n', 'var'), n = 1; end
    if ~exist('opt', 'var'), opt = []; end

    figure; insertAnnotation([runs.name '.secondary_vortices']);
    handles.hax(1) = subplot(3,2,[1 2]);
    hfield = runs.animate_field('csdye', handles.hax(1), tindex, 1, opt);
    xlim([-0.5 1]*220 + runs.eddy.mx(tindex)/1000);
    ylim([runs.bathy.xsb/1000-30 120+runs.eddy.my(tindex)/1000]);
    hfield.hcb.LimitsMode = 'auto';
    hfield.hcb.Position(1) = hfield.hcb.Position(1) + 0.04;
    %hfield.hfield.CData = hfield.hfield.CData - runs.bathy.xsb/1000;
    %caxis(hfield.hcb.Limits - runs.bathy.xsb/1000);
    %center_colorbar(hfield.hcb);
    %index = find(hfield.hcb.Ticks == 0);
    %hfield.hcb.TickLabels{index} = '0';
    colormap(handles.hax(1), cbrewer('seq', 'Greys', 32));
    hfield.hcen.Marker = '.';
    hfield.hcen.MarkerSize = 30;
    handles.hfield = hfield;
    handles.hax(1).Title.String = ['(a) ' handles.hax(1).Title.String];
    handles.hfield.hcb.TickDirection = 'out';

    if size(n,2) ~= 2
        [xx,yy] = ginput(n);
    else
        xx = n(:,1);
        yy = n(:,2);

        n = size(n, 1);
    end

    try
        xx = [xx; runs.eddy.mx(tindex)/1000];
        yy = [yy; runs.eddy.my(tindex)/1000];
    catch ME
        xx = [xx runs.eddy.mx(tindex)/1000];
        yy = [yy runs.eddy.my(tindex)/1000];
    end
    n = n + 1;

    %colors = brighten(cbrewer('qual', 'Paired', n), -0.5);
    % colors = [136 204 238; ...
    %           51  34 136; ...
    %           17 199  51; ...
    %           68 170 153; ...
    %           204 102 119; ...
    %           170  68 153]/255;

    colors = flipud([ 146,197,222; ...
                      244,165,130; ...
                      67,147,195; ...
                      214,96,77; ...
                      33,102,172; ...
                      178,24,43;]/255);

    pvback = runs.params.phys.f0 * runs.params.phys.N2/ ...
             runs.params.phys.g;

    tic;
    disp('Reading data');
    for nn = 1:n
        x = xx(nn); y = yy(nn);
        if nn ~= n
            handles.hpt(nn) = plot(x,y,'.','Color',colors(nn,:),'MarkerSize', 30);
        end

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

    [zrho,~] = runs.predict_zpeak(3, 'detect');
    zrho = abs(zrho);
    backup = get(groot, 'DefaultAxesColorOrder');
    set(groot, 'DefaultAxesColorOrder', colors);

    handles.hax(2) = subplot(323);
    handles.hcsd = plot(csdye - runs.bathy.xsb/1000, zr);
    uistack(handles.hcsd(end), 'bottom');
    [handles.hl(1,:), handles.htxt(1,:)] = common(runs, tindex, zrho);
    xlabel('(b) Cross-shelf dye - X_{sb} (km)');
    ylabel('Z (m)');

    handles.hax(3) = subplot(324);
    handles.hrho = plot(rho-rback, zr);
    [handles.hl(2,:), handles.htxt(2,:)] = common(runs, tindex, zrho);
    handles.hax(3).XTickLabelRotation = 0;
    handles.hax(3).XTickMode = 'auto';
    xlabel('(c) \Delta\rho (kg/m^3)')
    reduceSubplotHorizontalSpace(handles.hax(2:3));

    CenterProfileHandles = [handles.hfield.hcen handles.hcsd(end) handles.hrho(end)];
    if ~nopv
        handles.hax(4) = subplot(325);
        handles.hpv = plot(pv-pvback, zpv);
        [handles.hl(3,:), handles.htxt(3,:)] = common(runs, tindex, zrho);
        xlabel('(d) PV - f_0N^2/g');
        ylabel('Z (m)');
        CenterProfileHandles = [CenterProfileHandles handles.hpv(end)];

        handles.hax(5) = subplot(326);
        handles.hrv = plot(rv/runs.params.phys.f0, zpv);
        [handles.hl(4,:), handles.htxt(4,:)] = common(runs, tindex, zrho);
        xlim([-1 1] * max(abs(xlim)));
        xlabel('(e) (Relative vorticity)/f_0');

        reduceSubplotHorizontalSpace(handles.hax(4:5));
        CenterProfileHandles = [CenterProfileHandles handles.hrv(end)];
    end

    linkprop(CenterProfileHandles, 'Color');
    hfield.hcen.Color = [1 1 1] * 0.5;
    linkaxes(handles.hax(2:end), 'y');

    ylim([-1.3 0] * runs.eddy.Lgauss(tindex));
    set(groot, 'DefaultAxesColorOrder', backup);

    for ii=1:length(handles.htxt(:))
        handles.htxt(ii).Units = 'normalized';
        handles.htxt(ii).Position(1) = 0.9;
    end

    % extend horizontal lines
    for ii=1:4
        try
            axes(handles.hax(ii+1));
            for jj=1:size(handles.hl,2)-1
                handles.hl(ii,jj).XData = xlim;
            end
        catch ME
        end
    end
end

function [hl, ht] = common(obj, tindex, zrho)
    [hl(1), ht(1)] = liney(-1 * obj.bathy.hsb, 'H_{sb}');
    [hl(2), ht(2)] = liney(-1 * obj.eddy.Lgauss(tindex), 'L_z');
    [hl(3), ht(3)] = liney(-1 * zrho, 'z_\rho');
    hl(4) = linex(0);
    uistack(hl, 'bottom');
    pbaspect([1.65 1 1]); %axis tight; axis square;
    beautify;
end

function reduceSubplotHorizontalSpace(ax)

    width = ax(1).Position(3);
    dx = ax(2).Position(1) - (ax(1).Position(1) + width);
    ax(1).Position(1) = ax(1).Position(1) + dx*3/4;
    ax(2).Position(1) = ax(2).Position(1) - dx*3/4;

end