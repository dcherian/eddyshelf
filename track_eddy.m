function [eddy] = track_eddy(dir1)

    if isobject(dir1)
        runobj = dir1;
        dir1 = runobj.dir;
        fnames = roms_find_file(dir1,'his');
        file = char([dir1 '/' char(fnames(1))]);
        N = runobj.rgrid.N;
        [xr,yr,zr,~,~,~] = dc_roms_var_grid(file,'temp');

        grd = runobj.rgrid;

        if isempty(runobj.zeta)
            runobj.read_zeta;
        end
        if isempty(runobj.usurf)
            runobj.read_velsurf;
        end
        if isempty(runobj.rhosurf)
            runobj.read_rhosurf;
        end

        zeta = runobj.zeta;
        u = runobj.usurf;
        v = runobj.vsurf;
        rho = runobj.rhosurf;

        rthresh.vor = runobj.eddy.drhothresh(1);
        rthresh.ssh = runobj.eddy.drhothreshssh(1);

        f = runobj.rgrid.f(2:end-1,2:end-1)';

        eddy = runobj.eddy;

        eddy.mask = [];
        eddy.vormask = [];

        eddy.h = runobj.bathy.h;
        eddy.t = runobj.time/86400;

        params = runobj.params;
    else
        if isdir(dir1)
            fnames = roms_find_file(dir1,'his');
            file = char([dir1 '/' char(fnames(1))]);
            [xr,yr,zr,~,~,~,grd] = dc_roms_var_grid(file,'temp');
            tic;
            N = size(zr,3);

            trange = [];
            zeta  = dc_roms_read_data(dir1,'zeta',trange,{},[],grd, ...
                                      'his');
            u     = dc_roms_read_data(dir1,'u',trange,{'z' N N},[],grd, ...
                                      'his');
            v     = dc_roms_read_data(dir1,'v',trange,{'z' N N},[],grd, ...
                                      'his');
            rho   = dc_roms_read_data(dir1,'rho',trange,{'z' N N},[],grd, ...
                                      'his');
        else
            fname = dir1;
            index = strfind(dir1,'/');
            dir1 = dir1(1:index(end));
            fnames = [];
            file = fname;
            [xr,yr,zr,~,~,~,grd] = dc_roms_var_grid(file,'temp');
            tic;
            zeta = double(ncread(fname,'zeta'));
            u     = squeeze(double(ncread(fname,'u',[1 1 size(zr,3) 1],[Inf Inf 1 Inf])));
            v     = squeeze(double(ncread(fname,'v',[1 1 size(zr,3) 1],[Inf Inf 1 Inf])));
            toc;
        end
        f = grd.f(2:end-1,2:end-1)';
        params = read_params_from_ini(dir1);
        eddy.h = grd.h';
        eddy.t = dc_roms_read_data(dir1,'ocean_time', [], {}, [], ...
                                   grd, 'his')/86400; % required
                                                      % only for dt
        rthresh.vor = [];
        rthresh.ssh = [];
    end

    if strfind(file, 'his')
        tracer = 'rho';
    else
        tracer = 'temp';
    end

    kk = 2; % if ncread fails, it will use fnames(kk,:)
    tt0 = 0; % offset for new history.average file - 0 initially, updated later

    % search region for tracking eddies (in addition to detected diameter)
    limit_x = 40*1000;
    limit_y = 40*1000;

    % support for both cyclones and anti-cyclones
    sgn = sign(params.eddy.tamp);
    bathyloc = params.bathy.loc;

    %dx = xr(2,1,1) - xr(1,1,1);
    %dy = yr(1,2,1) - yr(1,1,1);

    zeta = sgn * zeta(2:end-1,2:end-1,:);
    vor = sign(params.phys.f0) * sgn * avg1(avg1( ...
                    bsxfun(@rdivide, diff(v,1,1), diff(avg1(xr(:,:,1),2),1,1)) ...
                  - bsxfun(@rdivide, diff(u,1,2), diff(avg1(yr(:,:,1),1),1,2)) ...
                  ,1),2);

    xr   = xr(2:end-1,2:end-1,end);
    yr   = yr(2:end-1,2:end-1,end);

    % make anomaly
    rho  = rho(2:end-1,2:end-1,:) - rho(1,1,1);

    % initial guess for vertical scale fit
    if ~isfield(params.flags,'vprof_gaussian') || params.flags.vprof_gaussian
        initGuess2(2) = params.eddy.depth;
        initGuess3(2) = params.eddy.depth;
        if ~isfield(params,'phys')
            T0 = 20; N2 = 1e-5; g = 9.81; TCOEF = 1.7e-4;
        else
            T0 = params.phys.T0;
            N2 = params.phys.N2;
            g  = params.phys.g;
            TCOEF = params.phys.TCOEF;
        end
        initGuess2(1) = T0;
        initGuess3(1) = T0;
        initGuess3(3) = N2./g./TCOEF;
    else
        % for sine fits
        initGuess(1) = params.eddy.tamp;
        initGuess(2) = 2*pi/params.eddy.depth;
        try
            initGuess(3) = params.eddy.theta0;
        catch ME
            initGuess(3) = pi/2;
        end

        initGuessGauss = [params.phys.T0 params.eddy.depth];
    end

    % detect shelfbreak
    [sbreak,~,~,ax] = find_shelfbreak(file);

    % remove background flow contribution to zeta
    if ax == 'x'
        zeta_bg = zeta(:,end,1);
    else
        zeta_bg = zeta(1,:,1);
    end

    eddy.xr = xr;
    eddy.yr = yr;

    % strategy is to interpolate zeta, vorticity onto uniform grid in the
    % horizontal - then I don't have to recode all the builtin image
    % processing functions that I make use of
    % suffix 'i' indicates interpolated variables
    dxi = min(diff(xr(:,1,1),1,1));
    dyi = min(diff(yr(1,:,1),1,2));
    xrivec = min(xr(:)):dxi:max(xr(:));
    yrivec = min(yr(:)):dyi:max(yr(:));
    [xri,yri] = meshgrid(xrivec,yrivec);
    % for interp2 later
    [xrgrd,yrgrd] = meshgrid(xr(:,1),yr(1,:));
    sz = size(xri');
     %%
    for tt=1:size(zeta,3)
        if tt == 1,
            mask = ones(sz);
            d_sbreak = Inf;
            thresh = nan;
            cx0 = nan;
            cy0 = nan;
        else
            if tt == 394,
                %keyboard % for debugging
            end
            mask = nan(sz);
            lx = eddy.dia(tt-1)/2 + limit_x;
            ly = eddy.dia(tt-1)/2 + limit_y;
            ix1 = find_approx(xri(1,:),eddy.cx(tt-1)-lx,1);
            ix2 = find_approx(xri(1,:),eddy.cx(tt-1)+lx,1);
            iy1 = find_approx(yri(:,1),eddy.cy(tt-1)-ly,1);
            iy2 = find_approx(yri(:,1),eddy.cy(tt-1)+ly,1);

            thresh = eddy.thresh(tt-1);
            cx0 = eddy.vor.cx(tt-1);
            cy0 = eddy.vor.cy(tt-1);
            mask(ix1:ix2,iy1:iy2) = 1;
            % distance to shelfbreak in *m*
%            d_sbreak = eddy.cx(tt-1)-sbreak;
        end

        % interpolate to denser grid
        izeta = interp2(xrgrd,yrgrd,bsxfun(@minus,zeta(:,:,tt), zeta_bg)', ....
                        xrivec', yrivec)';
        ivor = interp2(xrgrd, yrgrd, vor(:,:,tt)', xrivec', ...
                       yrivec)';
        irho = interp2(xrgrd, yrgrd, rho(:,:,tt)', xrivec', yrivec)';

        % find eddy using surface signatures
        fprintf('tt = %3d | ', tt);
        grd.dx = dxi; grd.dy = dyi;
        grd.cxn1 = cx0; grd.cyn1 = cy0;

        temp = eddy_diag(izeta .* mask, ...
                         ivor .* mask, irho, ...
                         grd,sbreak,thresh,[],bathyloc,rthresh); %w(:,:,tt));

        % if eddy detection terminates for whatever reason before
        % the entire simulation has been processed.
        if isempty(temp) || ~isfield(temp.vor, 'mask')
            warning('Eddy not found. Terminating!');
            eddy.tend = tt - 1;
            eddy.t = eddy.t(1:eddy.tend);
            break;
        end

        % let's interpolate the masks back to the coarser grid
        imask = interp2(xri,yri,temp.mask',xrgrd,yrgrd,'nearest')';
        ivormask = interp2(xri,yri,temp.vor.mask',xrgrd,yrgrd, ...
                           'nearest')';
        irhovormask = interp2(xri,yri,temp.rhovor.mask',xrgrd,yrgrd, ...
                           'nearest')';
        irhosshmask = interp2(xri,yri,temp.rhossh.mask',xrgrd,yrgrd, ...
                           'nearest')';

        % do more diagnostics
        % [cx,cy] = location of weighted center (first moment)
        eddy.cx(tt)  = temp.cx;
        eddy.cy(tt)  = temp.cy;
        % [mx,my] = location of maximum closest to shelfbreak
        eddy.mx(tt)  = temp.mx;
        eddy.my(tt)  = temp.my;
        % amplitude defined as max - mean amplitude around perimeter
        eddy.amp(tt) = temp.amp * sgn;
        % diameter of circle with same area
        eddy.dia(tt) = temp.dia;
        % mask for diagnostic purposes
        eddy.mask(:,:,tt) = imask;
        eddy.thresh(tt) = temp.thresh;
        eddy.vormask(:,:,tt) = ivormask;
        % number of pixels in eddy
        eddy.n(tt) = temp.n;
        % north, south, west & east edges
        eddy.we(tt) = temp.we;
        eddy.ee(tt) = temp.ee;
        eddy.ne(tt) = temp.ne;
        eddy.se(tt) = temp.se;

        eddy.Ls(tt) = temp.Ls;

        V = hypot(avg1(u(:,2:end-1,tt), 1) .* ivormask, ...
                  avg1(v(2:end-1,:,tt), 2) .* ivormask);
        eddy.V(tt) = max(V(:));

        thrname = {'vor', 'rhovor', 'rhossh'};
        fldname = {'we', 'ee', 'ne', 'se', 'lmin', 'lmaj', 'angle', ...
                   'dia', 'cx', 'cy', 'amp'};
        for mmm = 1:length(thrname)
            for nnn = 1:length(fldname)
                eval(['eddy.' thrname{mmm} '.' fldname{nnn} '(tt)' ...
                      ' = ' ...
                      'temp.' thrname{mmm} '.' fldname{nnn} ';']);
            end
        end

        eddy.vor.mask(:,:,tt) = ivormask;
        eddy.rhovor.mask(:,:,tt) = irhovormask;
        eddy.rhossh.mask(:,:,tt) = irhosshmask;

        eddy.vor = extradiags(eddy.vor, tt, vor, u, v, f, grd);
        eddy.rhovor = extradiags(eddy.rhovor, tt, vor, u, v, f, grd);
        eddy.rhossh = extradiags(eddy.rhossh, tt, vor, u, v, f, grd);

        if isempty(rthresh.ssh)
            rthresh.ssh = squeeze(nanmax(nanmax( ...
                rho(:,:,1) .* fillnan(eddy.mask,0), [], 1), [], 2));

            rthresh.vor = squeeze(nanmax(nanmax( ...
                rho(:,:,1) .* fillnan(eddy.vor.mask,0), [], 1), [], 2));
        end

        eddy.rthresh = rthresh;

        % diagnose vertical scale (fit Gaussian / sine)
        imx = find_approx(xr(:,1),eddy.mx(tt),1);
        imy = find_approx(yr(1,:),eddy.my(tt),1);
        ze  = squeeze(zr(imx,imy,:)); % z co-ordinate at center of eddy
        try
            eddy.T(tt,:) = double(squeeze(ncread(file,tracer,[imx imy 1 tt-tt0], ...
                            [1 1 Inf 1])));
        catch ME
            disp([' Moving to next file tt = ' num2str(tt) ' - ' char(fnames(kk))]);
            file = [dir1 '/' char(fnames(kk))];
            kk = kk +1;
            tt0 = tt-1;
            eddy.T(tt,:) = double(squeeze(ncread(file,tracer,[imx imy 1 tt-tt0],[ ...
                                                    1 1 Inf 1])));
        end

        try
            eddy.dyecen(tt,:) = double(squeeze(ncread(file, 'dye_03', [imx imy 1 tt-tt0], ...
                                                      [1 1 Inf 1])));
        catch ME
        end

        eddy.zT(tt,:) = ze;
        if params.bathy.axis == 'x'
            Ti = double(squeeze(ncread(file, tracer, ...
                                       [imx  size(xr,2)  1 tt-tt0],[1 1 Inf 1])));
        else
            Ti = double(squeeze(ncread(file, tracer, ...
                                       [size(xr,1)  imy  1 tt-tt0],[1 1 Inf 1])));
        end
        if strcmpi(tracer, 'rho')
            eddy.T(tt,:) = params.phys.T0 - 1./params.phys.TCOEF * ...
                    (1000+eddy.T(tt,:) - params.phys.R0)/params.phys.R0;
            Ti = params.phys.T0 - 1./params.phys.TCOEF * ...
                 (1000 + Ti - params.phys.R0)/params.phys.R0;
        end

        % let's save anomaly instead
        eddy.T(tt,:) = eddy.T(tt,:) - Ti';
        opts = optimset('MaxFunEvals',1e5);
        if ~isfield(params.flags,'vprof_gaussian') || params.flags.vprof_gaussian
            [x2,~,exitflag] = fminsearch(@(x) gaussfit2(x,eddy.T(tt,:)',ze), ...
                initGuess2,opts);
            if ~exitflag
                x2(2) = NaN;
                warning('Vertical scale = nan');
            end
            %[x3,~,exitflag] = fminsearch(@(x) gaussfit3(x,eddy.T(tt,:)',ze),initGuess3,opts);
            %if ~exitflag, x3(2) = NaN; end
            eddy.Lz2(tt)  = abs(x2(2));
            eddy.Lgauss(tt)  = abs(x2(2));
            %eddy.Lz3(tt)  = NaN;%abs(x3(2));
        else
            %fit sinusoid
            [x2,~,exitflag] = fminsearch(@(x) sinefit(x,eddy.T(tt,:)',ze), ...
                                initGuess,opts);
            if ~exitflag, x2(2) = NaN; end
            eddy.Lz2(tt) = abs(2*pi/x2(2));
            % save gaussian fit too
            [x2,~,exitflag] = fminsearch(@(x) gaussfit2(x,eddy.T(tt,:)',ze), ...
                        initGuessGauss,opts);
            if ~exitflag, x2(2) = NaN; end
            eddy.Lgauss(tt) = abs(x2(2));
            %eddy.Lz3(tt) = NaN;
        end

        % calculate center velocity
        if tt == 1
            eddy.mvx(1) = NaN;
            eddy.mvy(1) = NaN;
            eddy.cvx(1) = NaN;
            eddy.cvy(1) = NaN;
        else
            % dt in days; convert dx,dy to km
            dt = eddy.t(tt) - eddy.t(tt-1);
            eddy.mvx(tt) = (eddy.mx(tt) - eddy.mx(tt-1))./dt/1000;
            eddy.mvy(tt) = (eddy.my(tt) - eddy.my(tt-1))./dt/1000;

            % dt in days; convert dx,dy to km
            eddy.cvx(tt) = (eddy.vor.cx(tt) - eddy.vor.cx(tt-1))./dt/1000;
            eddy.cvy(tt) = (eddy.vor.cy(tt) - eddy.vor.cy(tt-1))./dt/1000;

            % repeat for ρ contour
            eddy.rhovor.cvx(tt) = (eddy.rhovor.cx(tt) - eddy.rhovor.cx(tt-1))./dt/1000;
            eddy.rhovor.cvy(tt) = (eddy.rhovor.cy(tt) - eddy.rhovor.cy(tt-1))./dt/1000;

            % repeat for ρ contour
            eddy.rhossh.cvx(tt) = (eddy.rhossh.cx(tt) - eddy.rhossh.cx(tt-1))./dt/1000;
            eddy.rhossh.cvy(tt) = (eddy.rhossh.cy(tt) - eddy.rhossh.cy(tt-1))./dt/1000;

        end
    end
    %%
    eddy.Lz2(abs(eddy.Lz2) > max(abs(zr(:)))) = NaN;
    %eddy.Lz3(abs(eddy.Lz3) > max(abs(zr(:)))) = NaN;
    toc;

    eddy.tmat = repmat(eddy.t', [1 size(eddy.T,2)]);

    eddy.comment = ['(cx,cy) = Location of weighted center (m) | ' ...
                    '(mx,my) = Location of SSH max closest to shelfbreak (m) | ' ...
                    'amp = amplitude (m) | dia = diameter of circle with same area (m) | ' ...
                    'mask = SSH mask to check detection | n = number of pixels | ' ...
                    '(mvx,mvy) = velocity of (mx,my) in km/day | ' ...
                    '(we,ee,ne,se) = West, East, North and South edges of eddy (m) | ' ...
                    'Lgauss = vertical scale (m) when fitting Gaussian - happens with sine fits too | ' ...
                    'Lz2,3 = Vertical scale (m) when fitting without & with linear trend | ' ...
                    'T = temp profile at (mx,my) | L = equiv diameter for vorticity < 0 region '...
                    'Lmin/Lmaj = minor/major axis length | Ls = ' ...
                    'speed based definition in Chelton et al. (2011) ' ...
                    '| tmat = time vector in t x z array form | zT ' ...
                    '= zvector at eddy center.'];

    if exist('runobj', 'var')
        runobj.eddy = eddy;
    end

    % save git hash
    eddy.hash = githash([mfilename('fullpath') '.m']);

    save([dir1 '/eddytrack.mat'],'eddy', '-v7.3');

    disp('Done.');

% Gaussian fit for vertical scale - called by fminsearch
% calculates squared error
function [E] = gaussfit2(x0,T,zr)
    % x = (T0,H,a)
    T0 = x0(1); h = x0(2);

    %E = sum((T - T0 * (1+a*zr) .* (1 + exp(-(zr/h).^2))).^2);
    E = sum((T - T0 .* exp(-(zr/h).^2)).^2);

function [E] = gaussfit3(x0,T,zr)
    % x = (T0,H,a)
    T0 = x0(1); h = x0(2); a = x0(3);

    E = sum((T - (T0+a*zr) .* (exp(-(zr/h).^2))).^2);

function [E] = sinefit(x0,T,zr)
    T0 = x0(1); k = x0(2); theta_0 = x0(3);
    E = sum((T - T0 * (1+sin(-k/4*zr + theta_0))/2).^2);

% Calculates eddy diagnostics as in Chelton et al. (2011)
% doesn't support multiple eddies yet
% only finds anticyclonic eddies
% this routine is called at every timestep
function [eddy] = eddy_diag(zeta, vor, rho, ...
                            grd, sbreak, thresh, w, bathyloc, ...
                            rthresh)

    dx = grd.dx; dy = grd.dy;
    flag_found = 0;

    % algorithm options
    opt.amp_thresh = 0.001; % Amplitude threshold (in m)
    % minimum eddy rad. = 5 km, maximum = 100 km
    opt.low_n  = floor(pi*(5e3)^2/dx/dy);       % minimum number of pixels in eddy
    opt.high_n = floor(pi*(200e3)^2/dx/dy);     % maximum number of pixels in eddy
    opt.connectivity = 8;  % either 4 or 8
    opt.max_dist = 400*1000;
    if isnan(thresh)
        thresh_min = nanmin(zeta(:));
    else
        thresh_min = thresh;
    end
    thresh_loop = linspace(thresh_min,nanmax(zeta(:)),10); % in m
    %thresh_loop = linspace(nanmin(zeta(:)),nanmax(zeta(:)),10);
    %%
    for ii=1:length(thresh_loop)
        threshold = thresh_loop(ii);

        % Criterion 1 - first get everything above threshold
        mask = zeta > threshold;

        % first find simply connected regions
        regions = bwconncomp(mask,opt.connectivity);

        % sort regions by n
        nn = [];
        for ll=1:regions.NumObjects
            nn(ll) = length(regions.PixelIdxList{ll});
        end
        [~,ind] = sort(nn,'descend');

        % loop over these regions and apply the criteria
        for kk=1:regions.NumObjects
            jj = ind(kk);
            % Criterion 2 - either too big or too small
            n = length(regions.PixelIdxList{jj});
            if n < opt.low_n || n > opt.high_n, continue; end

            % reconstruct zeta for identified region
            maskreg = zeros(regions.ImageSize);
            maskreg(regions.PixelIdxList{jj}) = 1;
            zreg   = repnan(zeta .* maskreg,0); % zeta in region
            zperim = zeta .* bwmorph(maskreg,'remove'); % zeta in perimeter

            % Criterion 3 - need local maximum in the
            % see http://stackoverflow.com/questions/1856197/how-can-i-find-local-maxima-in-an-image-in-matlab
            grd.local_max = imregionalmax(zreg);
            if grd.local_max == zeros(size(zreg)), continue; end
            %if sum(grd.local_max(:)) > 1, continue; end % skip if more than one maximum

            props = regionprops(maskreg,zeta,'WeightedCentroid','Solidity', ...
                        'EquivDiameter','Area','BoundingBox', ...
                        'MinorAxisLength','MajorAxisLength');

            % Criterion 4 - amplitude is at least > opt.amp_thresh
            indices = find(grd.local_max == 1);
            if length(indices) > 20; continue; end
            % make sure all local maxima found satisfy this criterion
            for mm=1:length(indices)
                if (zreg(indices(mm)) - nanmean(zperim(:))) < opt.amp_thresh
                   grd.local_max(indices(mm)) = 0;
                end
            end
            amp = max(zreg(grd.local_max)) - nanmean(zperim(:));
            if isempty(amp) || (amp < opt.amp_thresh) ...
                    || isequal(grd.local_max,zeros(size(grd.local_max)))
                continue;
            end

            % Criterion 5 - distance between any pair of points must be
            % less than a specified maximum
            % need to trace out outline, first find one point to start from
            % (convex hull doesn't do the job)
            [x0,y0] = ind2sub(size(maskreg),find(maskreg == 1,1,'first'));
            % Trace boundary points - sometimes the initial
            % starting direction doesn't work, so I try all 4 in
            % succession to find one that does
            try
                points = bwtraceboundary(maskreg,[x0,y0],'E');
            catch ME
                try
                    points = bwtraceboundary(maskreg,[x0,y0],'W');
                catch ME
                    try
                        points = bwtraceboundary(maskreg,[x0,y0], ...
                                                 'N');
                    catch ME
                        try
                            points = bwtraceboundary(maskreg,[x0,y0], ...
                                                     'S');
                        catch ME
                            disp(['All 4 bwtraceboundary directions ' ...
                                  'failed. This should not happen.']);
                            rethrow(ME);
                        end
                    end
                end
            end
            points = bsxfun(@times,points,[dx dy]);
            if maximumCaliperDiameter(points) > opt.max_dist, continue; end
           % if minimumCaliperDiameter(points) < min_dist, continue; end

            % Criterion 6 - Obuko-Weiss parameter must make sense
            if exist('w','var') && ~isempty(w)
                thresh = 0.2 * std(w(:)); % see Isern-Fontanet et al. (2006)
                wreg = w .* maskreg;
                if mean(wreg(:)) > thresh
                    break
                end
            end

            % Criterion 7 - solidity must be good - helps get rid of some
            % thin 'isthumuses'
            %if props.Solidity  < 0.75, continue; end

            % Criterion 8 - low 'rectangularity' - gets rid of rectangles
            % that are sometime picked up
            rectarea = props.BoundingBox(3) * props.BoundingBox(4);
            if props.Area/rectarea > 0.85, continue; end

            % Calculate properties of region
            cx = dx/2 + props.WeightedCentroid(2) * dx;
            cy = dy/2 + props.WeightedCentroid(1) * dy;

            % make grid index matrices
            grd.ix = bsxfun(@times,ones(size(zeta)),[1:size(zeta,1)]');
            grd.iy = bsxfun(@times,ones(size(zeta)),[1:size(zeta,2)]);

            [~,imax] = max(zeta(indices));
            imax = indices(imax);
            grd.indx = grd.ix(imax);
            grd.indy = grd.iy(imax);

            % I have an eddy!!!
            imagesc(zreg');

            xmax = fillnan(maskreg(:).*grd.ix(:),0);
            ymax = fillnan(maskreg(:).*grd.iy(:),0);

            % store eddy properties for return
            eddy.cx   = cx; % weighted center
            eddy.cy   = cy; %     "
            eddy.mx   = grd.ix(grd.indx,grd.indy) * dx; % maximum
            eddy.my   = grd.iy(grd.indx,grd.indy) * dy; %    "

            if eddy.my < sbreak && bathyloc ~= 'h';
                warning('eddy moving below shelfbreak!');
            end
            eddy.we   = dx/2 + nanmin(xmax) * dx; % west edge
            eddy.ee   = dx/2 + nanmax(xmax) * dx; % east edge
            eddy.ne   = dy/2 + nanmax(ymax) * dy; % south edge
            eddy.se   = dy/2 + nanmin(ymax) * dy; % north edge
            eddy.amp  = amp;
            eddy.dia  = props.EquivDiameter * sqrt(dx*dy);
            eddy.mask = maskreg;
            eddy.n    = n; % number of points
            eddy.jj   = jj;
            eddy.thresh = threshold;

            % find 0 rel. vor (max. speed) contour within SSH mask
            try
                eddy.vor = detect_eddy(vor.*eddy.mask < 0, zeta, opt, ...
                                       grd);

                % drhothresh based on ssh mask if it doesn't exist
                if isempty(rthresh.ssh)
                    rthresh.ssh = squeeze(nanmax(nanmax( ...
                        rho .* fillnan(eddy.mask,0), [], 1), [], 2));

                    rthresh.vor = squeeze(nanmax(nanmax( ...
                        rho .* fillnan(eddy.vor.mask,0), [], 1), [], 2));
                end
                eddy.rhovor = detect_eddy(rho < rthresh.vor, zeta, opt, grd);
                eddy.rhossh = detect_eddy(rho < rthresh.ssh, zeta, ...
                                          opt, grd);
                flag_found = 1;
                fprintf('Eddy found with threshold %.3f \n', threshold);
            catch ME
                flag_found = 0;
            end

            % stop when eddy is found
            if flag_found, break; end
        end

        if exist('flag_found','var') && flag_found == 1
            break;
        end
    end

    % eddy was found but i'm testing what Chelton's Ls would look like
    % calculate geostrophic velocity magnitude
    % then do thresholding and find SSH contour with max. avg speed.
    % then figure out equivDiameter for that threshold
    if flag_found
        ugeo = -1 * 9.81 .* 1/1e-4 * diff(zeta,1,2) / dy;
        vgeo =      9.81 .* 1/1e-4 * diff(zeta,1,1) / dx;

        V = hypot( avg1(ugeo,1), avg1(vgeo,2) );
        Vmax = 0;
        threshmax = 0;
        thresh_loop = linspace(threshold, nanmax(zeta(:)), 10);
        for iii=1:length(thresh_loop)-1
            zmask = (zeta .* eddy.mask) > thresh_loop(iii);
            try
                [x0,y0] = ind2sub(size(zmask),find(zmask == 1,1,'first'));
                points = bwtraceboundary(zmask,[x0,y0],'E');

                Vboundary = V(sub2ind(size(V), points(:,1), points(:,2)));
                if mean(Vboundary) > Vmax
                    Vmax = mean(Vboundary);
                    threshmax = thresh_loop(iii);
                end
            catch ME
                %disp(ME)
            end
        end

        % reuse zmask
        zmask = (zeta .* eddy.mask) > threshmax;
        props = regionprops(zmask, 'EquivDiameter', 'Centroid');
        min_index = 1;
        if length(props) ~= 1
            warning(['multiple regions found while calculating Ls. ' ...
                     'choosing southernmost one']);

            % assume first region is southernmost and compare that
            % to rest in a loop
            % min_index is assigned earlier to allow for case when
            % length(props) == 1
            for zzz=2:length(props)
                if props(zzz).Centroid < props(min_index).Centroid
                    min_index = zzz;
                end
            end
        end

        eddy.Ls = props(min_index).EquivDiameter/2 * sqrt(dx*dy);
    else
        % eddy was not found
        eddy = [];
    end

    if ~exist('eddy','var')
        warning('Eddy not found! Going to smaller threshold');
        eddy.cx   = NaN;
        eddy.cy   = NaN;
        eddy.mx   = NaN;
        eddy.my   = NaN;
        eddy.we   = NaN;
        eddy.ee   = NaN;
        eddy.ne   = NaN;
        eddy.se   = NaN;
        eddy.amp  = NaN;
        eddy.dia  = NaN;
        eddy.mask = NaN;
        eddy.n    = NaN;
        eddy.jj   = NaN;
        eddy.L    = NaN;
        eddy.lmaj = NaN;
        eddy.lmin = NaN;
        eddy.vormask = NaN;
        eddy.thresh =  NaN;

        % let's try again with lower threshold = 0
        % unless we have already tried that
        if thresh_min > 0
            eddy = eddy_diag(zeta, vor, dx, dy, sbreak, 0, w, cxn1, cyn1);
        end
    end

    try
        rmfield(eddy,'jj');
    catch ME
    end

% extra diagnostics
function [in] = extradiags(in, tt, vor, u, v, f, grd)
    Roeddy = in.mask(:,:,tt) .* vor(:,:,tt) ./ f;
    dA = 1./grd.pm(2:end-1,2:end-1)' .* 1./grd.pn(2:end-1,2:end-1)' ...
         .* in.mask(:,:,tt);
    in.area(tt) = nansum(dA(:));
    in.Ro(tt) = abs(nansum(nansum(Roeddy .* dA, 1),2) / ...
                          in.area(tt));

    KE = (avg1(u(:,2:end-1,tt),1).^2 + avg1(v(2:end-1,:, tt),2).^2);
    in.Vke(tt) = sqrt(nansum(nansum(KE .* dA, 1), 2) / in.area(tt));

function [out] = detect_eddy(maskin, zeta, opt, grd)
% This routine finds a contour for a non-SSH field that contains
% the SSH maximum and satisfies Chelton's criteria
% Used for 0vor and ρ fields.
    regions = bwconncomp(maskin,opt.connectivity);
    nn = [];
    for zz=1:regions.NumObjects
        nvor(zz) = length(regions.PixelIdxList{zz});
    end
    [~,ind] = sort(nvor,'descend');

    dx = grd.dx; dy = grd.dy;

    % ugh. now I have to pick the right region, like earlier
    for mm=1:regions.NumObjects
        ll = ind(mm);
        maskreg = zeros(regions.ImageSize);
        maskreg(regions.PixelIdxList{ll}) = 1;
        % reconstruct zeta for identified region
        zreg   = repnan(zeta .* maskreg,0); % zeta in region
        zperim = zeta .* bwmorph(maskreg,'remove'); % zeta in perimeter

        % only 1 criterion for now
        n = length(regions.PixelIdxList{ll});
        if n < opt.low_n || n > opt.high_n, continue; end

        % make sure that region contains eddy.mx,eddy.my
        if maskreg(grd.ix(grd.indx,grd.indy), ...
                   grd.iy(grd.indx,grd.indy)) == 0
            continue;
        end

        % extract information
        props  = regionprops(maskreg,zeta,'EquivDiameter', ...
                             'MinorAxisLength','MajorAxisLength', ...
                             'BoundingBox', 'WeightedCentroid', ...
                             'Area', 'Solidity', 'Orientation');

        % store eddy properties for return
        out.cx = dx/2 + props.WeightedCentroid(2) * dx;
        out.cy = dy/2 + props.WeightedCentroid(1) * dy;

        % check displacement of center
        % should be less than 10 grid cells
        if ~isnan(out.cx) && ~isnan(out.cy)
            answer = 1;
            if hypot(out.cx-grd.cxn1, out.cy-grd.cyn1) > ...
                    10*hypot(dx,dy)
                disp(['eddy center > 10 dx from last time ' ...
                      'instant.']);

                clf;
                pcolorcen(zeta');
                clim = caxis;
                hold on;
                contour(maskreg', 1, 'k', 'LineWidth', 2);
                caxis(clim);
                plot(grd.cxn1/dx, grd.cyn1/dy, '*');
                plot(out.cx/dx, out.cy/dy, 'k*');
                legend('zeta','current contour', 'earlier center', ['present ' ...
                                    'center']);
                answer = input([' enter 0 to skip this region, 1 to ' ...
                                'accept (' ...
                                num2str(regions.NumObjects ...
                                        - mm) ' regions left): ']);
                if ~answer
                    disp('region skipped');
                    flag_found == 0;
                    continue;
                end
            end
        end

        xmax = fillnan(maskreg(:).*grd.ix(:),0);
        ymax = fillnan(maskreg(:).*grd.iy(:),0);

        out.we   = dx/2 + nanmin(xmax) * dx; % west edge
        out.ee   = dx/2 + nanmax(xmax) * dx; % east edge
        out.ne   = dy/2 + nanmax(ymax) * dy; % south edge
        out.se   = dy/2 + nanmin(ymax) * dy; % north edge
        out.amp  = max(zreg(grd.local_max) - nanmean(zperim(:)));
        out.dia  = props.EquivDiameter * sqrt(dx*dy);
        out.mask = maskreg;
        out.lmaj = props.MajorAxisLength * sqrt(dx*dy);
        out.lmin = props.MinorAxisLength * sqrt(dx*dy);
        out.angle = props.Orientation;
    end