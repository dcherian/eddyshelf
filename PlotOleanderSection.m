function [] = PlotOleanderSection(tind, xind, hfig)
    if ~exist('hfig', 'var')
        hfig = figure;
        maximize;
    end

    fname = '../data/All_Oleander_3D.mat';
    ole = load(fname);
    zint = [-500:1:0];

    % do linear vertical interpolation
    Tint = nan([size(ole.Temp_Fix, 2) length(zint)]);
    for aa=1:size(ole.Temp_Fix, 2) % along-shelf direction
        if isempty(cut_nan(squeeze(ole.Depth(tind,aa,:))))
            break;
        end
        [zvec,uind] = unique(cut_nan(squeeze(ole.Depth(tind,aa,:))));
        Tvec = cut_nan(squeeze(ole.Temp_Fix(tind,aa,uind)));
        if length(zvec) > 2 & ~isempty(Tvec)
            Tint(aa,:) = interp1(-1*abs(zvec), Tvec, zint);
        end
    end

    xactual = ole.Dist(tind,1:size(Tint,1))/1000;
    xinterp = [nanmin(xactual):2:nanmax(xactual)];
    Tinterp = nan([length(xinterp) size(Tint,2)]);
    if ole.Year(tind,1) > 2007
        Lgauss = nanmedian(diff(xactual));
    else
        Lgauss = 16;
    end

    for zz=1:size(Tint,2)
        if ~all(isnan(Tint(:,zz)))
            [Tinterp(:,zz),~] = OA1(xinterp, Lgauss, 300, 1, 0.05, xactual, Tint(:,zz), 30);
        end
    end

    insertAnnotation(['PlotOleanderSection.m(' num2str(tind) ',[' num2str(xind) '])']);
    ax(1) = subplot(1,3,[1 2]); cla;
    contourf(xinterp, zint, Tinterp', 20);
    %contourf(xactual, zint, Tint', 20);
    hold on;
    hline = plot(ole.Dist(tind,:)/1000, ole.WaterDepth(tind,:), 'k');
    lon = cut_nan(ole.Lon(tind,:));
    if lon(end) < lon(1)
        set(gca, 'xdir', 'reverse');
    end
    linex(ole.Dist(tind,xind)/1000);
    colorbar; caxis([8 20]);
    ylabel('Z (m)'); xlabel('Along-track Distance (km)');
    beautify; hline.LineWidth = 4;

    ax(2) = axes('position', ax(1).Position);
    ax(2).XAxisLocation = 'top';
    ax(2).XDir = ax(1).XDir;
    ax(2).Box = 'off'; ax(1).Box = 'off';
    ax(2).YTick = [];
    linkaxes(ax(1:2), 'x');
    ax(2).XTick = unique(cut_nan(xactual));
    ax(2).XTickLabel = {'\nabla'};
    ax(2).XAxis.TickLabelGapMultiplier = 0;
    ax(2).XAxis.TickLabelGapOffset = -2;
    title(['XBT Data | Oleander | tind = ' num2str(tind) ' | ' ...
           num2str(nanmin(ole.Day(tind,:))) '-' ...
           num2str(nanmax(ole.Day(tind,:))) ' ' ...
           datestr([num2str(ole.Month(tind,1)) '/' num2str(ole.Month(tind,1))], 'mmm') ...
           ' ' num2str(ole.Year(tind,1))]);
    beautify;
    ax(2).XAxis.TickDirection = 'in';
    axes(ax(1));

    ax(3) = subplot(1,3,3); cla;
    hold on;
    for kk=1:length(xind)
        if ole.Year(tind,1) < 2008
            hh(kk) = plot(squeeze(ole.Temp_Fix(tind,xind(kk),:)), ...
                          -1*squeeze(ole.Depth(tind,xind(kk),:)), 'o', ...
                          'LineWidth', 2);
            plot(Tint(xind(kk),:), zint, 'Color', hh(kk).Color);
        else
            hh(kk) = plot(squeeze(ole.Temp_Fix(tind,xind(kk),:)), ...
                          -1*squeeze(ole.Depth(tind,xind(kk),:)), '-', ...
                          'LineWidth', 2);
        end
    end
    ax(3).XAxisLocation = 'top';
    xlabel('Temp (C)');
    legend(hh, cellstr(num2str(xind')), 'Location', 'SouthEast');
    beautify;

    linkaxes(ax, 'y');
    ylim([-500 0]);
end