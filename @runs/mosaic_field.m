% mosaic animate_field plots.
%     [ax] = mosaic_field(runs, varname, timesteps, opt)

function [handles] = mosaic_field(runs, varname, timesteps, opt)

    if ~exist('opt', 'var'), opt = []; end

    if ~isfield(opt, 'zetaOnFirstPlot')
        opt.zetaOnFirstPlot = 0;
    end

    N = ceil(length(timesteps)/2);
    letters = 'abcdefghijkl';

    figure; maximize;
    if N~=1
        handles.hax = packfig(2, N);
    else
        handles.hax = packfig(N, 2);
    end

    insertAnnotation([runs.name '.mosaic_field']);
    for ii=1:length(timesteps)
        tstep = runs.process_time(timesteps(ii));
        if opt.zetaOnFirstPlot && ii == 1
            opt.addzeta = 1;
        end

        if opt.zetaOnFirstPlot && ii ~=1
            opt.addzeta = 0;
        end

        handles.hfield{ii} = runs.animate_field(varname, handles.hax(ii), tstep, 1, opt);
        handles.hfield{ii}.htlabel.String = [letters(ii) ') ' ...
                            handles.hfield{ii}.htlabel.String];
        title('');

        if ii == 1 % save first colorbar limit
            clim = caxis;
        end

        if N ~= 1
            if ii <= N % take out x tick labels for top row
                handles.hax(ii).XTickLabel = {};
                xlabel('');
            end
            if  (ii~=1) & (ii~=N+1) % take out y tick labels for first column
                handles.hax(ii).YTickLabel = {};
                ylabel('');
            end
        else
            if  ii == 2 % take out y tick labels for first column
                handles.hax(ii).YTickLabel = {};
                ylabel('');
            end
        end

        if ii ~= N*2
            colorbar('off');
        else
            handles.hcb = colorbar;
        end
    end

    switch N
      case 1
        moveColorbarOut1x2(handles.hcb);
      otherwise
        moveColorbarOut2x2(handles.hcb);
    end

    if runs.bathy.axis == 'y' & N ~= 1
        moveSubplotsCloserInY(2, N, handles.hax);
    end

    linkaxes(handles.hax, 'xy'); axis tight;

    [handles.supax, handles.htitle] = suplabel(['Surface ' varname], 't');
    handles.supax.Position(4) = handles.hax(1).Position(end) ...
        + handles.hax(1).Position(2) - handles.supax.Position(2);
    handles.htitle.FontWeight = 'normal';
end
