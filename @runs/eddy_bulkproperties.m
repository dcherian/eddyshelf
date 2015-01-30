% eddy bulk properties - integrated PV, RV, volume, energy
function [] = eddy_bulkproperties(runs)
%%
    slab = 30; % read 10 at a time
    ftype = 'his';

    nt = length(runs.eddy.t);

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
    zr = permute(runs.rgrid.z_r(:, 2:end-1, 2:end-1), [3 2 1]);

    pvname = [runs.dir '/ocean_vor.nc'];
    if exist(pvname,'file')
        dopv = 0;
    else
        dopv = 0;
    end

    % background density field
    if runs.bathy.axis == 'y'
        try
            tback = permute( dc_roms_read_data(dirname, 'temp', [1 1], ...
                                               {'x' 1 1; 'y' 2 sz4dfull(2)+1}, ...
                                               [], rgrid, ftype, ...
                                               'single'), [3 1 2]);
        catch ME
            rback = permute( dc_roms_read_data(dirname, 'rho', [1 1], ...
                                               {'x' 1 1; 'y' 2 sz4dfull(2)+1}, ...
                                               [], rgrid, ftype, ...
                                               'single'), [3 1 2]);
        end
    else
        error('Not implemented for N-S isobaths');
    end

    % figure out eddy structure in 3D
    % read surface density field
    rhosurf = dc_roms_read_data(runs.dir, 'rho', [], {'z' runs.rgrid.N runs.rgrid.N}, ...
                                [], runs.rgrid, ftype, 'single');
    rhosurf = rhosurf(2:end-1, 2:end-1, :);
    % find what density corresponds to 0 vorticity contour
    rhothresh = squeeze(max(max(rhosurf.*runs.eddy.vormask, [], 1), ...
                            [], 2));
    rhothresh = rhothresh(1);

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
        maskvor = sparse(reshape( repmat( ...
            permute(logical(repnan(vormask(:,:,tt:tend), 0)), [1 2 4 3]), ...
            [1 1 N 1]), sz));

        %vol{tt} = runs.domain_integratesp(masked.*maskvor, dVsp);
        % calculate total volume
        volcell{mm} = full(nansum( bsxfun(@times, masked.*maskvor, dVsp)));

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

        try
            temp = dc_roms_read_data(dirname, 'temp', ...
                                     [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
                                     [], rgrid, ftype, 'single');

            pe = double(- runs.params.phys.TCOEF* bsxfun(@times, ...
                                                         bsxfun(@minus, temp, tback), zr)  ...
                        .* runs.params.phys.g .* runs.params.phys.R0);
        catch ME
            rho  = dc_roms_read_data(dirname, 'rho', ...
                                     [tt tend],{'x' 2 sz4dfull(1)+1; 'y' 2 sz4dfull(2)+1}, ...
                                     [], rgrid, ftype, 'single');

            pe = double(bsxfun(@times, rho, zr) .* runs.params.phys.g);
        end

        masked = sparse(reshape(rho > rhothresh, sz));

        intpe{mm} = full(nansum( bsxfun(@times, ...
                                        masked.*maskvor.*reshape(pe, sz), dVsp)));

        intke{mm} = full(nansum( bsxfun(@times, ...
                                        masked.*maskvor.*reshape(0.5 ...
                                        * double((rho+1000) .* u.^2 + v.^2), sz), dVsp)));

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

    runs.eddy.vol = cell2mat(volcell);
    if dopv
        runs.eddy.PV = cell2mat(intpv);
        runs.eddy.RV = cell2mat(intrv);
    end
    runs.eddy.KE = cell2mat(intke);
    runs.eddy.PE = cell2mat(intpe);

    runs.eddy.hash = githash([mfilename('fullpath') '.m']);

    eddy = runs.eddy;
    save([runs.dir '/eddytrack.mat'],'eddy');
end