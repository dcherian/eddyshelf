function [] = PlotOleanderSection(tind, xind)
    fname = '../data/All_Oleander_3D.mat';
    ole = load(fname);
    zint = [-800:1:0];

    for aa=1:size(ole.Temp_Fix, 2) % along-shelf direction
        if isempty(cut_nan(squeeze(ole.Depth(tind,aa,:))))
            break;
        end
        Tint(aa,:) = interp1(-1*cut_nan(squeeze(ole.Depth(tind,aa,:))), ...
                             cut_nan(squeeze(ole.Temp_Fix(tind,aa,:))), ...
                             zint);
    end

    figure; maximize;
    insertAnnotation(['PlotOleanderSection.m(' num2str(tind) ',[' num2str(xind) '])']);
    ax(1) = subplot(1,3,[1 2]);
    contourf(ole.Dist(tind,1:size(Tint,1))/1000, zint, Tint', 40);
    lon = cut_nan(ole.Lon(tind,:));
    if lon(end) < lon(1)
        set(gca, 'xdir', 'reverse');
    end
    linex(ole.Dist(tind,xind)/1000);
    colorbar;
    ylabel('Z (m)'); xlabel('Along-track Distance (km)');
    title(['XBT Data | Oleander | ' ...
           num2str(nanmin(ole.Day(tind,:))) '-' ...
           num2str(nanmax(ole.Day(tind,:))) ' ' ...
           datestr([num2str(ole.Month(tind,1)) '/' num2str(ole.Month(tind,1))], 'mmm') ...
           ' ' num2str(ole.Year(tind,1))]);
    beautify;

    ax(2) = subplot(1,3,3);
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
    ax(2).XAxisLocation = 'top';
    xlabel('Temp (C)');
    legend(hh, cellstr(num2str(xind')), 'Location', 'SouthEast');
    beautify;

    linkaxes(ax, 'y');
    ylim([-800 0]);
end