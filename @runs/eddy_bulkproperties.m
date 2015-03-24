% eddy bulk properties - integrated PV, RV, volume, energy,
% bottom pressure torque, volume transport, mass transport
function [] = eddy_bulkproperties(runs, slab)
%%
    if ~exist('slab', 'var') || isempty(slab)
        slab = 10; % read 'n' at a time
    end

    ftype = 'his';
    rho_flag = 1; % read density, not temperature

    nt = length(runs.eddy.t);

    tind = [1 nt];

    sz4dfull = [fliplr(size(runs.rgrid.z_r))-[2 2 0] slab];
    sz4dsp = [prod(sz4dfull(1:3)) slab];
    sz4dspend = [sz4dsp(1) mod(nt,slab)];
    sz3dsp = [sz4dsp(1) 1];

    szpvfull = [fliplr(size(runs.rgrid.z_r))-[2 2 1] slab];
    szpvsp = [prod(szpvfull(1:3)) slab];
    szpvspend = [szpvsp(1) mod(nt,slab)];
    szpv3dsp = [szpvsp(1) 1];

    dVsp = reshape(runs.rgrid.dV(2:end-1,2:end-1,:), sz3dsp);
    dVpvsp = reshape( avg1(runs.rgrid.dV(2:end-1, 2:end-1,:),3), szpv3dsp);

    % store variables to optimize parfor loop
    rgrid = runs.rgrid;
    vormask = runs.eddy.vormask;
    eddname = runs.eddname;
    dirname = runs.dir;
    thresh = runs.eddy_thresh;
    N = runs.rgrid.N;
    g = runs.params.phys.g;
    rho0 = runs.params.phys.rho0;
    zr = permute(runs.rgrid.z_r(:, 2:end-1, 2:end-1), [3 2 1]);
    % bottom slope
    if runs.bathy.axis == 'y'
        slbot = avg1(diff(runs.bathy.h, 1, 2)./diff(runs.rgrid.yr,1,2), 2);
        slbot = slbot(2:end-1,:);
    else
        slbot = avg1(diff(runs.bathy.h, 1, 1)./diff(runs.rgrid.xr,1,1), 1);
        slbot = slbot(:,2:end-1);
    end

    f = runs.rgrid.f(2:end-1,2:end-1)';

    pvname = [runs.dir '/ocean_vor.nc'];
    if exist(pvname,'file')
        dopv = 0;
    else
        dopv = 0;
    end

    % background density field
    if runs.bathy.axis == 'y'
        if ~rho_flag
            tback = permute( dc_roms_read_data(dirname, 'temp', [1 1], ...
                                               {'x' 1 1; 'y' 2 sz4dfull(2)+1}, ...
                                               [], rgrid, ftype, ...
                                               'single'), [3 1 2]);
        else
            rback = permute( dc_roms_read_data(dirname, 'rho', [1 1], ...
                                               {'x' 1 1; 'y' 2 sz4dfull(2)+1}, ...
                                               [], rgrid, ftype, ...
                                               'single'), [3 1 2]);
        end
    else
        error('Not implemented for N-S isobaths');
    end

    if ~isfield(runs.eddy, 'drhothresh')
        % figure out eddy structure in 3D
        % read surface density field
        runs.read_rhosurf;
        rhosurf = runs.rhosurf(2:end-1,2:end-1,:);
        drhosurf = bsxfun(@minus, rhosurf, rback(:,:,end));
        % find what density corresponds to 0 vorticity contour
        runs.eddy.drhothresh = squeeze(nanmax(nanmax(drhosurf.* ...
                                                     fillnan(runs.eddy ...
                                                          .vormask,0), ...
                                                     [], 1), [], ...
                                              2));
        runs.eddy.rhothreshssh = squeeze(nanmax(nanmax(drhosurf.* ...
                                           fillnan(runs.eddy.mask,0), ...
                                           [], 1), [], 2));

        % rhosurf is not needed anymore
        clear rhosurf;

        debug = 0;
        if debug

            % calculate ertel pv
            %[pv,xpv,ypv,zpv] = roms_pv(runs.dir, [], {'z' 70 70},
            %'ocean_pv.nc', 'his')'
            pv = squeeze(ncread([runs.dir '/ocean_pv.nc'], 'pv'));

            runs.eddy.pvthresh = squeeze(nanmean(nanmean(fillnan(pv(:,:,1).* ...
                                                              bwmorph(runs.eddy.vormask(:,:,1), 'remove'), ...
                                                              0), 1), 2));

            % debugging plots
            debug = 0;
            if debug
                tt = 1;
                var = pv; %drhosurf;
                hplt = pcolorcen(runs.rgrid.x_rho(2:end-1, 2:end-1)/1000, ...
                                 runs.rgrid.y_rho(2:end-1, 2:end-1)/1000, ...
                                 var(:,:,tt)');
                clim = caxis;
                hedd = runs.plot_eddy_contour('contour',tt);
                hssh = runs.plot_eddy_sshcontour('contour',tt);
                [~,hpv] = contour(runs.rgrid.x_rho(2:end-1, 2:end-1)/1000, ...
                                  runs.rgrid.y_rho(2:end-1, 2:end-1)/1000, ...
                                  pv(:,:,tt)', [1 1] * runs.eddy.pvthresh(1), ...
                                  'g', 'LineWidth', 2);
                %vend = var(:,:,end);
                %caxis([min(vend(:)) max(vend(:))]);
                caxis(clim);

                for tt=1:2:size(var, 3)
                    set(hplt, 'CData', double(var(:,:,tt)'));
                    runs.update_eddy_contour(hedd, tt);
                    runs.update_eddy_sshcontour(hssh, tt);
                    set(hpv, 'zdata', pv(:,:,tt)');
                    pause(0.1);
                end
            end
        end
    end

    % initial time instant works well - see plot_eddye
    rhothreshes = [runs.eddy.drhothresh(1) ...
                   runs.eddy.drhothreshssh(1)];

    % allocate variables
    %intpe = cell(ceil(nt/slab), length(rhothreshes));
    %intke = intpe;
    volcell = cell(ceil(nt/slab), length(rhothreshes));
    fullmasscell = cell(ceil(nt/slab), length(rhothreshes));
    anommasscell = cell(ceil(nt/slab), length(rhothreshes));
    voltrans = cell(ceil(nt/slab), length(rhothreshes));
    masstrans = cell(ceil(nt/slab), length(rhothreshes));
    btrq = cell(ceil(nt/slab), length(rhothreshes));
    %cor = cell(ceil(nt/slab), length(rhothreshes));

    ticstart = tic;
    for mm=1:ceil(nt/slab)
        tt = (mm-1)*slab + 1;
        disp([' mm= ' num2str(mm) '/' num2str( ...
            ceil(nt/slab))]);
        tend = tt + slab - 1;
        if tend > nt
            tend = nt;
            sz = sz4dspend;
            szpv = szpvspend;
        else
            sz = sz4dsp;
            szpv = szpvsp;
        end
        disp([tt tend]);
        %eddye = dc_roms_read_data(dirname, eddname, ...
        %                          [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
        %                          [],rgrid, ftype, 'single'); %#ok<*PROP>

        %masked  = sparse(reshape(eddye > thresh, sz));
        %maskvor = sparse(reshape( repmat( ...
        %    permute(logical(repnan(vormask(:,:,tt:tend), 0)), [1 2 4 3]), ...
        %    [1 1 N 1]), sz));

        maskvor = 1;

        % integrated energies
        % yes, using sz4dfull(1)+1 IS correct. sz4dfull has interior
        % RHO point counts
        if isnan(runs.params.bg.ubt)
            runs.params.bg.ubt = 0;
        end
        if isnan(runs.params.bg.vbt)
            runs.params.bg.vbt = 0;
        end

        u = avg1(dc_roms_read_data(dirname, 'u', ...
                                   [tt tend],{'y' 2 sz4dfull(2)+1}, ...
                                   [],rgrid, ftype, 'single'),1) - runs.params.bg.ubt; %#ok<*PROP>
        v = avg1(dc_roms_read_data(dirname, 'v', ...
                                   [tt tend],{'x' 2 sz4dfull(1)+1}, ...
                                   [],rgrid, ftype, 'single'),2) - runs.params.bg.vbt; %#ok<*PROP>

        if ~rho_flag
            temp = dc_roms_read_data(dirname, 'temp', ...
                                     [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
                                     [], rgrid, ftype, 'single');

            pe = double(- runs.params.phys.TCOEF* bsxfun(@times, ...
                                                         bsxfun(@minus, temp, tback), zr)  ...
                        .* g .* runs.params.phys.R0);
        else
            rho  = dc_roms_read_data(dirname, 'rho', ...
                                     [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
                                     [], rgrid, ftype, 'single');
        end

        for nnn = 1:length(rhothreshes)
            drhothresh = rhothreshes(nnn);

            disp(['rho threshold ' num2str(nnn) '/' ...
                  num2str(length(rhothreshes))]);

            ranom = bsxfun(@minus, rho, rback);
            masked = sparse(reshape(ranom < drhothresh, sz));

            pe = -1 * double(bsxfun(@times, ranom, zr) .* g);

            intpe{mm, nnn} = full(nansum( bsxfun(@times, ...
                                                 masked.*maskvor.*reshape(pe, sz), dVsp)));

            intke{mm, nnn} = full(nansum( bsxfun(@times, ...
                                                 masked.*maskvor.* ...
                                                 reshape(0.5 * ...
                                                         double((rho+1000) .* ...
                                                              (u.^2 + v.^2)), sz), dVsp)));

            % bottom pressure
            %btrq{mm, nnn} = full(nansum(bsxfun(@times, masked .* maskvor .* ...
            %                            reshape(bsxfun(@times, double(ranom./rho0).*g, ...
            %                                      slbot), sz), ...
            %                                   dVsp)));
            % coriolis term
            %cor{mm, nnn} = full(nansum(bsxfun(@times, masked.*maskvor.* ...
            %                         reshape(bsxfun(@times, double(u), ...
            %                                        f), sz), dVsp)));

            % transport
            voltrans{mm, nnn} = full(nansum(bsxfun(@times, ...
                                                   masked.*maskvor.* ...
                                                   reshape(abs(double(u)), sz), dVsp)));
            masstrans{mm, nnn} = full(nansum(bsxfun(@times, ...
                                                    masked.*maskvor.* ...
                                                    reshape(double((rho+1000).*abs(u)), ...
                                                            sz), dVsp)));

            % calculate total volume
            volcell{mm, nnn} = full(nansum( bsxfun(@times, masked.*maskvor, ...
                                                   dVsp)));
            fullmasscell{mm, nnn} = full(nansum( bsxfun(@times, ...
                                                     masked .* maskvor .* ...
                                                     reshape(double(rho+1000), ...
                                                             sz), ...
                                                        dVsp)));
            % eddy anomaly mass
            anommasscell{mm, nnn} = full(nansum( bsxfun(@times, ...
                                                        masked .* maskvor .* ...
                                                        reshape(double(ranom), ...
                                                              sz), dVsp)));
        end

        % integrated PV, RV
        if dopv
            disp('Reading pv, rv');
            pv = double(ncread(pvname, 'pv',[1 1 1 tt],[Inf Inf Inf tend-tt+1])); %#ok<*PROP>
            rv = double(avg1(avg1(ncread(pvname, 'rv',[1 1 1 tt], ...
                                         [Inf Inf Inf tend-tt+1]),1),2));
            % pv,rv are at N-1 levels in vertical, so we need
            % to calculate masks again
            masked  = sparse(reshape(avg1(eddye,3) > thresh, szpv));
            maskvor = sparse(reshape( repmat( ...
                permute(logical(vormask(:,:,tt:tend)), [1 2 4 3]), ...
                [1 1 N-1 1]), szpv));

            intpv{mm} = full(nansum( bsxfun(@times, ....
                                            masked.*maskvor.*reshape(pv, szpv), dVpvsp))) ...
                ./ volcell{mm};
            intrv{mm} = full(nansum( bsxfun(@times, ....
                                            masked.*maskvor.*reshape(rv, szpv), dVpvsp))) ...
                ./ volcell{mm};
        end
    end
    toc(ticstart);

    % save data to structure
    runs.eddy.fullmass = cell2mat(fullmasscell')';
    runs.eddy.vol = cell2mat(volcell')';
    runs.eddy.voltrans = cell2mat(voltrans')';
    runs.eddy.masstrans = cell2mat(masstrans')';
    runs.eddy.mass = cell2mat(anommasscell')';
    runs.eddy.btrq = cell2mat(btrq')';
    % runs.eddy.cor = cell2mat(cor')';
    if dopv
        runs.eddy.PV = cell2mat(intpv')';
        runs.eddy.RV = cell2mat(intrv')';
    end
    runs.eddy.KE = cell2mat(intke')';
    runs.eddy.PE = cell2mat(intpe')';

    runs.eddy.threshes = rhothreshes;
    runs.eddy.bulkhash = githash([mfilename('fullpath') '.m']);

    runs.eddy.comment = [runs.eddy.comment ' | btrq has 1/ρ0 in it ' ...
                        '| full mass = (1000+ρ) | anommass = mass ' ...
                        'of eddy anomaly | energies are calculated ' ...
                        'used (1000+ρ)'];
    eddy = runs.eddy;
    save([runs.dir '/eddytrack.mat'],'eddy');
end