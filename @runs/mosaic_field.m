% mosaic animate_field plots.
%     [ax] = mosaic_field(runs, varname, timesteps, opt)

function [handles] = mosaic_field(runs, varname, timesteps, opt, ...
                                  hax)

    if ~exist('opt', 'var'), opt = []; end

    if ~isfield(opt, 'topdown')
        opt.topdown = 0;
    end
    if ~isfield(opt, 'zetaOnFirstPlot')
        opt.zetaOnFirstPlot = 0;
    end

    Nt = length(timesteps);
    N = ceil(Nt/2);
    letters = 'abcdefghijkl';

    if ~exist('hax', 'var')
        figure; maximize;
        if Nt == 2 | Nt == 3
            if opt.topdown
                handles.hax = packfig(Nt, 1);
            else
                handles.hax = packfig(1, Nt);
            end
        else
            handles.hax = packfig(2, N);
        end
    else
        handles.hax = hax;
    end

    if Nt == 2 | Nt == 3
        N = 1;
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

        if ~opt.topdown
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
                if  ii > 1 % take out y tick labels for first column
                    handles.hax(ii).YTickLabel = {};
                    ylabel('');
                end
            end
        else
            % top down single column
            if ii < Nt
                handles.hax(ii).XTickLabels = {};
                handles.hax(ii).XLabel.String = '';
            end
        end
        
        if ii ~= Nt
            colorbar('off');
        else
            handles.hcb = colorbar;
        end
    end

    if N == 1
        moveColorbarOut1x2(handles.hcb);
    else
        moveColorbarOut2x2(handles.hcb);

        if runs.bathy.axis == 'y'
            % moveSubplotsCloserInY(2, N, handles.hax);
        end
    end

    linkaxes(handles.hax, 'xy'); axis tight;

    [handles.supax, handles.hsuptitle] = suplabel(['Surface ' varname], 't');
    handles.supax.Position(4) = handles.hax(1).Position(end) ...
        + handles.hax(1).Position(2) - handles.supax.Position(2);
    handles.hsuptitle.FontWeight = 'normal';

end
