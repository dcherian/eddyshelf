function [handles] = PlotSingleXZSection(runs, varname, loc, day, opt, hax)

    if ~exist('opt', 'var'), opt = []; end
    if ~isfield(opt, 'rhocontours')
        opt.rhocontours = 0;
    end
    if ~isfield(opt, 'maskcontour')
        opt.maskcontour = 1;
    end

    varname = runs.process_varname(varname);

    if ~exist('day', 'var'), day = []; end
    if ischar(day)
        day = str2double(day);
        tindex = vecfind(runs.time, day*86400);
    else
        tindex = day;
    end

    if ~exist('loc', 'var')
        loc = runs.eddy.my(tindex);
    end
    if ischar(loc)
        loc = str2double(loc);
    else
        if loc <= length(runs.csflux.x)
            % should be integer isobath
            isobath = loc;
            if runs.params.flags.flat_bottom * isobath > 1
                isobath = isobath - 1;
            end
            loc = runs.csflux.x(isobath);
            if isempty(tindex)
                [~,tindex] = runs.calc_maxflux( ...
                    runs.recalculateFlux(2*runs.bathy.hsb, isobath, isobath));
            end
        else
            if runs.bathy.axis == 'y'
                loc = runs.rgrid.y_rho(loc,1);
            else
                loc = runs.rgrid.x_rho(1,loc);
            end
        end
    end

    if runs.bathy.axis == 'y'
        sx1 = runs.spng.sx1;
        sx2 = runs.spng.sx2;
        ix = find_approx(runs.rgrid.y_rho(:,1), loc);
        bathyax = 2; asax = 1;
        xvec = runs.rgrid.x_rho(1,2:end-1) - runs.eddy.mx(tindex);
        zvec = runs.rgrid.z_r(:, ix, 1);
    else
        sx1 = runs.spng.sy1;
        sx2 = runs.spng.sy2;
        ix = find_approx(runs.rgrid.x_rho(1,:), loc);
        bathyax = 1; asax = 2;
        xvec = (runs.rgrid.y_rho(2:end-1,1) - runs.eddy.my(tindex))';
        zvec = runs.rgrid.z_r(:, 1, ix);
    end

    xvec = xvec(sx1:sx2);
    L = runs.eddy.rhovor.dia(tindex)/2;
    xsb = runs.bathy.xsb;

    if strcmpi(varname, 'rhoanom')
        rhoanomflag = 1;
        rback = dc_roms_read_data(runs.dir, 'rho', ...
                                  1, {runs.bathy.axis ix ix}, ...
                                  [], runs.rgrid, 'his', 'single');
        varname = 'rho';
    else
        rhoanomflag = 0;
    end

    if strcmpi(varname, runs.csvelname)
        var = squeeze(avg1( ...
            dc_roms_read_data(runs.dir, runs.csvelname, tindex, ...
                              {runs.bathy.axis ix-1 ix}, [], runs.rgrid, ...
                              'his', 'single'), bathyax));

    else
        if strcmpi(varname, 'w')
            var = squeeze(avg1( ...
                dc_roms_read_data(runs.dir, 'w', tindex, ...
                                  {runs.bathy.axis ix ix}, [], runs.rgrid, ...
                                  'his', 'single'), 2));
        else
            var = dc_roms_read_data(runs.dir, varname, ...
                                    tindex, {runs.bathy.axis ix ix}, ...
                                    [], runs.rgrid, 'his', 'single');
        end
    end

    if rhoanomflag
        var = var - rback(sx1:sx2,:);
    end

    if opt.rhocontours
        rho = dc_roms_read_data(runs.dir, 'rho', ...
                                    tindex, {runs.bathy.axis ix ix}, ...
                                [], runs.rgrid, 'his', 'single');
        rho = rho(sx1:sx2,:);
    end

    if strcmpi(varname, runs.csdname)
        csdye = var;
    else
        csdye = dc_roms_read_data(runs.dir, runs.csdname, ...
                                  tindex, {runs.bathy.axis ix ix}, ...
                                  [], runs.rgrid, 'his', 'single');
    end

    var = var(sx1:sx2,:);
    csdye = csdye(sx1:sx2,:);

    %if isobath ~= 1
    if runs.sgntamp > 0
            mask = fillnan(bsxfun(@times, csdye < runs.csflux.x(isobath), ...
                                  runs.csflux.offmask(sx1:sx2,tindex,isobath)),0);
        else
            mask = fillnan(bsxfun(@times, csdye > runs.csflux.x(isobath), ...
                                  runs.csflux.offmask(sx1:sx2,tindex,isobath)),0);
        end
        %end

    if ~exist('hax', 'var')
        figure; hax = gca;
    else
        axes(hax);
    end

    insertAnnotation([runs.name '.PlotSingleXZSection']);
    maximize;

    [~,handles.hfield] = contourf(xvec/1000, zvec, var', 20);
    hold on;
    clim = caxis;
    if opt.rhocontours
        [~,handles.hrhocont] = contour(xvec/1000, zvec, rho',10, 'k');
        caxis(clim);
    end
    handles.hfield.EdgeColor = 'none';
    %if isobath ~= 1
    [~,handles.hmask] = ...
        contour(xvec/1000, zvec, repnan(mask',0), [1 1], 'k', 'LineWidth', 2);
    caxis(clim);
        %end
    handles.htitle(1) = title([runs.name ' | ' varname]);
    xlabel('X - X_{eddy} (km)'); ylabel('Z (m)');

    if ~runs.params.flags.flat_bottom
        [handles.hline, handles.htext] = ...
            liney(-1 * [runs.eddy.Lgauss(tindex) runs.bathy.hsb], ...
                              {'vertical scale'; 'h_{sb}'});
    else
        [handles.hline, handles.htext] = ...
            liney(-1 * runs.eddy.Lgauss(tindex), 'vertical scale');
    end

    handles.htime = runs.add_timelabel(tindex);
    handles.htime.Position = [0.68 0.13];
    if exist('isobath', 'var')
        handles.htime.String = [{handles.htime.String;  ...
                            ['y/R = ' ...
                            num2str(runs.csflux.ndloc(isobath), '%.2f')]}];
    end

    beautify;
    if exist('isobath', 'var') & (isobath == 1)
        ylim([min(zvec(:)) 0]);
    end

    if strcmpi(varname, runs.csvelname)
        title('Cross-shelf velocity');
        caxis([-1 1] * max(abs(var(:))));
        handles.hcb = center_colorbar;
    end

    if strcmpi(varname, runs.csdname)
        title('Cross-shelf dye(km)');

        % muck with colorbar
        cmin = round(min(runs.sgntamp*var/1000));
        cmax = round(max(runs.sgntamp*var/1000));
        handles.hcb.TickLabelsMode = 'auto';
        handles.hcb.TickDirection = 'out';
        handles.hcb.Limits = [cmin cmax];
        handles.hcb.Ticks = sort([cmin ...
                            round((runs.bathy.xsl + xsb)/2000 - xsb/1000) ...
                            round(runs.params.eddy.cy/1000 - xsb/1000)]);
        % handles.hcb(4).TickLabels{1} = ['Shelfbreak - ' num2str(-1*handles.hcb(4).Ticks(1)) ' km'];
        handles.hcb.TickLabels{1} = 'Shelf Water';
        % handles.hcb(4).TickLabels{3} = 'Shelfbreak';
        handles.hcb.TickLabels{2} = 'Slope Water';
        handles.hcb.TickLabels{end} = 'Eddy Water';
    end
end