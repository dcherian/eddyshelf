function mosaic_zslice(runs, varname, depth, tind)

    n = length(tind);

    figure;
    insertAnnotation([runs.name '.mosaic_zslice(' varname ',' ...
                      num2str(depth) ')']);
    for ii=1:n
        hax(ii) = subplot(2,ceil(n/2),ii);
        runs.animate_zslice(varname, depth, tind(ii), hax(ii));
        title(['z = -' num2str(depth) ' m']);
    end

    linkaxes(hax, 'xy');
end