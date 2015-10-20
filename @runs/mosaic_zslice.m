%        [handles] = mosaic_zslice(runs, varname, depth, tind)
function [handles] = mosaic_zslice(runs, varname, depth, tind)

    n = length(tind);

    figure;
    insertAnnotation([runs.name '.mosaic_zslice(' varname ',' ...
                      num2str(depth) ')']);
    hax = packboth(2,ceil(n/2));
    for ii=1:n
        axes(hax(ii));
        handles(ii) = runs.animate_zslice(varname, depth, tind(ii), hax(ii));
        title('');
        colorbar('off');
    end

    linkaxes(hax, 'xy');

    [handles(1).supax, handles(1).htitle] = ...
        suplabel([runs.name ' | ' varname ' | z = -' num2str(abs(depth)) 'm'], 't');
    handles(1).supax.Position(4) = 0.85;

    set(hax(2), 'YTickLabel', {}, 'XTickLabel', {});
    xlabel(hax(1), ''); xlabel(hax(2), '');
    ylabel(hax(2), ''); ylabel(hax(4), '');
    set(hax(4), 'YTickLabel', {});
    set(hax(1), 'XTickLabel', {});

    handles(1).hcb = colorbar;
    handles(1).hcb.Position(1) = 0.92;
    handles(1).hcb.Position(2) = 0.35;

    axes(hax(1));
    handles(1).htext(1) = text(0.05, 0.83, 'Surface cross-shelf dye', 'FontSize', 13, ...
                               'Units', 'Normalized', 'Color', [1 1 1]*0.65);
    handles(1).htext(2) = text(0.05, 0.79, 'Surface eddy dye', 'FontSize', 13, ...
                               'Units', 'Normalized', 'Color', 'k');
    handles(1).htext(3) = text(0.05, 0.75, 'Eddy core', 'FontSize', 13, ...
                               'Units', 'Normalized', 'Color',[44 162 95]/255);
end