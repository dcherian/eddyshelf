function [hplt] = overlay_section(runs, varname1, varname2, tindex, volume, axis, index, ...
                                  commands)

    if ~exist('commands', 'var'), commands = []; end

    hplt = runs.plot_section(varname1, tindex, volume, axis, index, commands);
    hplt.h_plot.EdgeColor = 'none';
    insertAnnotation([runs.name '.overlay_section']);
    hold on

    titlestr = [runs.name ' | ' varname1 ' (color) | ' varname2 ...
                ' (contour) | ' axis ' = ' num2str(index) ' | t = ' ...
                num2str(runs.time(tindex)/86400) ' days'];

    if axis == 's', axis = 'z'; end
    volume{size(volume,1)+1,1} = axis;
    volume{size(volume,1),2} = index;
    volume{size(volume,1),3} = index;

    [v2, xax, yax, zax] = dc_roms_read_data(runs, varname2, tindex, volume);

    switch axis
      case 'x'
        plotx = yax;
        ploty = zax;

      case 'y'
        plotx = xax;
        ploty = zax;

      case 'z'
        plotx = xax;
        ploty = yax;
    end

    [~,hplt.h_plot2] = contour(plotx/1000, ploty/1000, v2, 50, 'Color', 'k');

    hplt.h_title.String = titlestr;
end
