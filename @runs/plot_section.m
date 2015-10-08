% use mod_movie to generate section plots (SINGLE TIME INSTANT)
function [] = plot_section(runs,varname, tindex, volume, axis, index, commands)

    mod_movie(runs, varname, tindex, volume, axis, index, commands);
    insertAnnotation([runs.name '.plot_section']);

    if axis == 'z' || axis == 's'
        runs.plot_bathy('contour');
        hrho = runs.plot_rho_contour('contour', tindex);
        plot(runs.eddy.mx/1000, runs.eddy.my/1000);
        plot(runs.eddy.mx(tindex)/1000, runs.eddy.my(tindex)/1000, 'kx');
    end
end