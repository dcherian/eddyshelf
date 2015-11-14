% mosaic hovmoeller plots
%    [ax] = mosaic_hovmoeller(runArray, varname, axname, loc, iz)
function [ax] = mosaic_hovmoeller(runArray, varname, axname, loc, iz)

    if ~exist('iz', 'var'), iz = []; end

    N = length(runArray.filter);
    figure; maximize; insertAnnotation('runArray.mosaic_hovmoeller');
    ax = packboth(2,N/2);

    clim = [];
    for ii=1:N
        ff = runArray.filter(ii);
        handles(ii) = runArray.array(ff).hovmoeller(varname, axname, loc, iz, ax(ii));
        str = ax(ii).Title.String;
        title('');
        text(0.50, 0.90, str, 'Units', 'normalized');
        clim = [clim caxis];

        if ii == 1 | ii == 3
            ax(ii).XTickLabel = '';
        end
    end
    linkaxes(ax, 'xy');

    clim = [min(clim) max(clim)];
    for ii=1:N
        caxis(clim);
    end
end
