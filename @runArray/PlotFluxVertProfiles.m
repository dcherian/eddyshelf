function [handles] = PlotFluxVertProfiles(runArray, fontSizes, hax)

    if ~exist('fontSizes', 'var'), fontSizes = []; end
    if ~exist('hax', 'var'), hax = []; end

    highlight = ['ew-8392'; 'ew-8381']; hili = [];
    isobath = 1;

    redmap = brighten(cbrewer('seq', 'Greys', runArray.len), -0.5);
    bluemap = brighten(cbrewer('seq', 'Greys', runArray.len), -0.5);
    flatcolor = [0 0 0]; [27,158,119]/255;

    if isempty(hax)
        figure; % maximize;
        handles.hax(1) = subplot(121); hold on;
        handles.hax(2) = subplot(122); hold on;
    else
        handles.hax = hax;
    end

    insertAnnotation('runArray.PlotFluxVertProfiles');

    phio = runArray.print_params('bathy.hsb./(V0./bathy.S_sh/sqrt(phys.N2))');
    runArray.sort(phio);
    phio = runArray.print_params('bathy.hsb./(V0./bathy.S_sh/sqrt(phys.N2))');
    set(handles.hax(1), 'ColorOrder', redmap);

    corder_backup = runArray.sorted_colors;

    kk = 1;
    axes(handles.hax(1));
    colormap(handles.hax(1), redmap(1:end-4,:));
    for ii=1:length(runArray.array)
        run = runArray.array(ii);

        hsb = run.bathy.hsb;

        if run.params.misc.rdrg ~= 0
            continue;
        end

        if run.bathy.sl_shelf == 0 & ~strcmpi(run.name, 'ew-34')
            continue;
        end

        [start, stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
        % offflux = run.csflux.off.slopezt(:,start:stop,1,1); %
        % fullsb
        try
            offflux = squeeze(trapz(run.radius.xvec(:, start), ...
                                    run.radius.off.shelfpos, 1)); % one radius
        catch ME
            continue;
        end
        zivec = run.csflux.vertbins(:,isobath);
        profile = trapz(run.csflux.time(start:stop)*86400, ...
                        offflux, 2);
        profile = profile./max(abs(profile));

        legstr1{kk} = num2str(phio(ii), '%.2f');
        handles.hplt1(kk) = plot(profile, zivec./hsb);
        if run.bathy.sl_shelf == 0 & strcmpi(run.name, 'ew-34')
            handles.hplt1(kk).Color = flatcolor;
            hflat(1) = handles.hplt1(kk);
        end
        if strcmpi(run.name, highlight(1,:)) ...
                | strcmpi(run.name,  highlight(2,:))
            hili = [hili, handles.hplt1(kk)];
        end

        kk = kk+1;
    end

    % phii = runArray.print_params('(bathy.hsb./(V0./bathy.S_sl/sqrt(phys.N2)))');
    % chi = runArray.print_params(['(2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
    %                     'V0/Lz0/(bathy.S_sl*sqrt(phys.N2))']);

    % runArray.sort(phii);
    % phii = runArray.print_params('(bathy.hsb./(V0./bathy.S_sl/sqrt(phys.N2)))');

    % chi = runArray.print_params(['(2/sqrt(pi)*exp(-(bathy.hsb/Lz0)^2)) *' ...
    %                     'V0/Lz0/(bathy.S_sl*sqrt(phys.N2))']);

    ssl = runArray.print_params('bathy.S_sl');
    runArray.sort(ssl);
    ssl = runArray.print_params('bathy.S_sl');
    set(handles.hax(2), 'ColorOrder', bluemap);

    kk = 1;
    axes(handles.hax(2)); hold on;
    colormap(handles.hax(2), bluemap(1:end-4,:));
    for ii=1:length(runArray.array)
        run = runArray.array(ii);

        if run.params.misc.rdrg ~= 0
            continue;
        end

        if run.bathy.sl_shelf == 0 & ~strcmpi(run.name, 'ew-34')
            continue;
        end

        hsb = run.bathy.hsb;

        %[start, stop] = run.flux_tindices(run.csflux.on.slope(:,1,1), 0.2, 0.9);
        %onflux = run.csflux.on.slopeztneg(:,start:stop,1,1); % full sb

        [start, stop] = run.flux_tindices(run.csflux.off.slope(:,1,1));
        try
            onflux = squeeze(trapz(run.radius.xvec(:, start), ...
                                   run.radius.on.nonshelfneg, 1));
            % one radius
        catch ME
            continue;
        end

        zivec = run.csflux.vertbins(:,isobath);
        profile = trapz(run.csflux.time(start:stop)*86400, ...
                        onflux, 2);
        profile = abs(profile./max(abs(profile)));
        %        profile = profile./profile(end);

        % legstr2{kk} = ['(' num2str(phii(ii), '%.1f') ...
        %               ', ' num2str(chi(ii), '%.2f') ')'];

        legstr2{kk} = num2str(ssl(ii), '%.1f');
        handles.hplt2(kk) = plot(profile, zivec./hsb);
        if run.bathy.sl_shelf == 0
            handles.hplt2(kk).Color = flatcolor;
            hflat(2) = handles.hplt2(kk);
        end
        kk = kk+1;
    end

    if isempty(fontSizes)
        legfontsize = 20;
    else
        legfontsize = fontSizes(1);
    end

    axes(handles.hax(1));
    set(gca, 'XAxisLocation', 'Top');
    %liney(-1);%ylim([-1 0]);
    limx = xlim;
    xlim([0 limx(2)]);
    xlabel({'a) Shelf water outflow'; '(normalized transport)'});
    ylabel('Z / H_{sb}');
    axis square;
    handles.hcbar(1) = colorbar('southoutside');
    %    handles.hcbar(1).Position(2) = 0.1;
    handles.hcbar(1).Box = 'off';
    handles.hcbar(1).Label.String = '\phi_o';
    caxis([min(phio) max(phio)]);
    beautify(fontSizes, 'Times');
    % handles.hleg(1) = columnlegend(2, legstr1, 'FontSize', legfontsize-1, 'Location', 'NorthWest');
    % handles.hleg(1).Position(1) = 0.15;
    % handles.hleg(1).Position(2) = 0.23;
    % legpos = handles.hleg(1).Position;
    % handles.htxt(1) = text(legpos(1) + legpos(3), ...
    %                        (legpos(2) + legpos(4))*1.2, ...
    %                        '\phi_o', 'Units', 'normalized', ...
    %                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %                        'FontSize', legfontsize+2);

    axes(handles.hax(2));
    set(gca, 'XAxisLocation', 'Top');
    % liney(-1);
    limx = xlim;
    xlim([0 limx(2)]);
    xlabel({'b) Eddy & slope water inflow';  '(normalized transport)'});
    % handles.hax(2).YTickLabels = {};
    axis square;
    linkaxes(handles.hax, 'y');
    %handles.hax(2).Position(2) = 0.2;

    handles.hcbar(2) = colorbar('southoutside');
    %    handles.hcbar(2).Position(2) = 0.1;
    handles.hcbar(2).Label.String = 'S_{sl}';
    handles.hcbar(2).Box = 'off';
    caxis([min(ssl) 3]);
    beautify(fontSizes, 'Times');

    % handles.hleg(2) = columnlegend(2, legstr2, 'FontSize', legfontsize);
    % handles.hleg(2).Position(1) = 0.75;
    % handles.hleg(2).Position(2) = 0.1;
    % legpos = handles.hleg(2).Position;
    % handles.htxt(2) = text(legpos(1) + legpos(3)*1.25, ...
    %                        legpos(2) + legpos(4)*0.9, ...
    %                        'S_{sl}', 'Units', 'normalized', ...
    %                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %                        'FontSize', legfontsize+2);

    %    linkprop(handles.hleg, 'FontSize');
    % linkprop(handles.htxt, 'FontSize');

    handles.hflat = hflat;
    for hh = 1:length(hili)
        hili(hh).LineWidth = 2;
        hili(hh).Color = [1 0 0];
        uistack(hili(hh), 'top');
    end

    hflat(1).LineWidth = 2;
    hflat(2).LineWidth = 2;
    uistack(hflat(1), 'top');
    uistack(hflat(2), 'top');

end
