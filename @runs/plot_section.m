% use mod_movie to generate section plots (SINGLE TIME INSTANT)
function [hplt] = plot_section(runs,varname, tindex, volume, axis, ...
                               index, commands, hax)

    if ~exist('commands', 'var'), commands = []; end
    if ~exist('hax', 'var'), figure; hax = gca; end

    hplt = mod_movie(runs, varname, tindex, volume, axis, index, commands, hax);
    insertAnnotation([runs.name '.plot_section']);

    if axis == 'z' || axis == 's'
        hplt.h_bathy = runs.plot_bathy('contour');
        hplt.h_rho = runs.plot_rho_contour('contour', tindex);
        hplt.h_tlabel = runs.add_timelabel(tindex);
        plot(runs.eddy.mx/1000, runs.eddy.my/1000);
        plot(runs.eddy.mx(tindex)/1000, runs.eddy.my(tindex)/1000, 'kx');
    end
end