function [hplt] = overlay_section(runs, varname1, varname2, tindex, volume, axis, index, ...
                                  commands, hax)

    if ~exist('commands', 'var'), commands = ''; end
    if ~exist('hax', 'var'), figure; hax = gca; end

    tindex = runs.process_time(tindex);

    hplt = runs.plot_section(varname1, tindex, volume, axis, index, commands, hax);
    hplt.h_plot.EdgeColor = 'none';
    insertAnnotation([runs.name '.overlay_section']);
    hold on

    titlestr = [runs.name ' | ' varname1 ' (color) | ' varname2 ...
                ' (contour) | ' axis ' = ' num2str(index)];

    if axis == 's', axis = 'z'; end
    volume{size(volume,1)+1,1} = axis;
    volume{size(volume,1),2} = index;
    volume{size(volume,1),3} = index;

    [v2, xax, yax, zax] = dc_roms_read_data(runs, varname2, tindex, volume);

    switch axis
      case 'x'
        plotx = yax/1000;
        ploty = zax;

      case 'y'
        if size(xax,2) == 1
            plotx = repmat(xax, [1 size(zax,2)])/1000;
        else
            plotx = xax/1000;
        end
        ploty = zax;

      case 'z'
        plotx = xax/1000;
        ploty = yax/1000;
    end

    [~,hplt.h_plot2] = contour(plotx, ploty, v2, 20, 'Color', 'k');

    hplt.h_title.String = titlestr;
end
