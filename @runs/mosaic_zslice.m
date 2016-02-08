%        [handles] = mosaic_zslice(runs, varname, depth, tind)
function [handles] = mosaic_zslice(runs, varname, depth, tind, opt)

    if ~exist('opt', 'var'), opt = []; end

    n = length(tind);
    letters = 'abcdefghijkl';

    figure;
    insertAnnotation([runs.name '.mosaic_zslice(' varname ',' ...
                      num2str(depth) ')']);
    hax = packfig(2,ceil(n/2));
    for ii=1:n
        axes(hax(ii));
        handles(ii) = runs.animate_zslice(varname, depth, tind(ii), hax(ii), opt);
        title('');
        colorbar('off');
        handles(ii).htime.String = [letters(ii) ') ' ...
                            handles(ii).htime.String];

        uistack([handles(ii).csdsurf handles(ii).rhocont], 'top');
    end

    linkaxes(hax, 'xy');

    [handles(1).supax, handles(1).htitle] = ...
        suplabel([runs.name ' | ' varname ' | z = -' num2str(abs(depth)) 'm'], 't');
    handles(1).supax.Position(4) = 0.85;
    handles(1).htitle.FontWeight = 'normal';

    N = ceil(n/2);

    set(hax(1), 'XTickLabel', {});
    xlabel(hax(1), '');
    for ii=2:N
        set(hax(ii), 'YTickLabel', {}, 'XTickLabel', {});
        set(hax(ii+N), 'YTickLabel', {});

        xlabel(hax(ii), '');
        ylabel(hax(ii), ''); ylabel(hax(ii+N), '');
    end

    handles(1).hcb = colorbar;
    moveColorbarOut2x2(handles(1).hcb);

    pos = [0.83 0.77 0.71];
    if runs.sgntamp < 0
        pos = flip(1-pos);
    end

    axes(hax(1));
    handles(1).htext(1) = text(0.05, pos(1), 'Surface cross-shelf dye', 'FontSize', 16, ...
                               'Units', 'Normalized', 'Color', runs.shelfSlopeColor);
    handles(1).htext(3) = text(0.05, pos(2), 'Eddy core', 'FontSize', 16, ...
                               'Units', 'Normalized', 'Color', runs.rhoContColor);
    handles(1).htext(2) = text(0.05, pos(3), 'Surface eddy dye', 'FontSize', 16, ...
                               'Units', 'Normalized', 'Color', 'k');
end