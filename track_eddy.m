function [eddy] = track_eddy(dir1)

    if isdir(dir1)
        fnames = roms_find_file(dir1,'his');
        file = char([dir1 '/' char(fnames(1))]);
        file
        [xr,yr,zr,~,~,~] = dc_roms_var_grid(file,'temp');
        tic;
        disp('Reading data');
        zeta  = roms_read_data(dir1,'zeta');
        u     = roms_read_data(dir1,'u',[1 1 size(zr,3) 1],[Inf Inf 1 Inf]);
        v     = roms_read_data(dir1,'v',[1 1 size(zr,3) 1],[Inf Inf 1 Inf]);
        toc;
    else
        fname = dir1;
        index = strfind(dir1,'/');
        dir1 = dir1(1:index(end));
        fnames = [];
        file = fname;
        [xr,yr,zr,~,~,~] = dc_roms_var_grid(file,'temp');
        tic;
        disp('Reading data');
        zeta = double(ncread(fname,'zeta'));
        u     = squeeze(double(ncread(fname,'u',[1 1 size(zr,3) 1],[Inf Inf 1 Inf])));
        v     = squeeze(double(ncread(fname,'v',[1 1 size(zr,3) 1],[Inf Inf 1 Inf])));
        toc;
    end
    kk = 2; % if ncread fails, it will use fnames(kk,:)
    tt0 = 0; % offset for new history.average file - 0 initially, updated later
    
    % search region for tracking eddies (in addition to detected diameter)
    limit_x = 40*1000;
    limit_y = 40*1000;

    eddy.h = ncread(file,'h');
    eddy.t = roms_read_data(dir1,'ocean_time')/86400; % required only for dt
    dt = eddy.t(2)-eddy.t(1);
    dx = xr(2,1,1) - xr(1,1,1);
    dy = yr(1,2,1) - yr(1,1,1);

    zeta = zeta(2:end-1,2:end-1,:);
    vor = avg1(avg1(diff(v,1,1)./dx - diff(u,1,2)./dy,1),2);
    clear u v
    xr   = xr(2:end-1,2:end-1,end);
    yr   = yr(2:end-1,2:end-1,end);
    
    sz = size(zeta(:,:,1));
    
    % initial guess for vertical scale fit
    params = read_params_from_ini(dir1);
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
            initGuess(3) = 0;
        end
    end

    % detect shelfbreak
    [sbreak,~,~,ax] = find_shelfbreak(file);

    % remove background flow contribution to zeta
    if ax == 'x'
        zeta_bg = zeta(:,end,1);
    else
        zeta_bg = zeta(1,:,1);
    end
    
    for tt=1:size(zeta,3)
        if tt == 1,
            mask = ones(sz);
            d_sbreak = Inf;
        else 
            if tt ==  73,
                disp('debug time!');
            end
            mask = nan(sz);
            lx = eddy.dia(tt-1)/2 + limit_x;
            ly = eddy.dia(tt-1)/2 + limit_y;
            ix1 = find_approx(xr(:,1),eddy.cx(tt-1)-lx,1);
            ix2 = find_approx(xr(:,1),eddy.cx(tt-1)+lx,1);
            iy1 = find_approx(yr(1,:),eddy.cy(tt-1)-ly,1);
            iy2 = find_approx(yr(1,:),eddy.cy(tt-1)+ly,1);

            mask(ix1:ix2,iy1:iy2) = 1;
            % distance to shelfbreak in *m*
%            d_sbreak = eddy.cx(tt-1)-sbreak;
        end
        fprintf('tt = %3d | ', tt);
        temp = eddy_diag(bsxfun(@minus,zeta(:,:,tt),zeta_bg) .* mask, ...
                            vor(:,:,tt).*mask,dx,dy,sbreak); %w(:,:,tt));

        % [cx,cy] = location of weighted center (first moment)
        eddy.cx(tt)  = temp.cx;
        eddy.cy(tt)  = temp.cy;
        % [mx,my] = location of maximum closest to shelfbreak
        eddy.mx(tt)  = temp.mx;
        eddy.my(tt)  = temp.my;
        % amplitude defined as max - mean amplitude around perimeter
        eddy.amp(tt) = temp.amp;
        % diameter of circle with same area
        eddy.dia(tt) = temp.dia;
        % mask for diagnostic purposes
        eddy.mask(:,:,tt) = temp.mask;
        eddy.vormask(:,:,tt) = temp.vormask;
        % number of pixels in eddy
        eddy.n(tt) = temp.n;
        % north, south, west & east edges
        eddy.we(tt) = temp.we;
        eddy.ee(tt) = temp.ee;
        eddy.ne(tt) = temp.ne;
        eddy.se(tt) = temp.se;
        eddy.L(tt)  = temp.L;
        eddy.lmin(tt) = temp.lmin;
        eddy.lmaj(tt) = temp.lmaj;
        
        % diagnose vertical scale (fit Gaussian)
        imx = find_approx(xr(:,1),eddy.mx(tt),1);
        imy = find_approx(yr(1,:),eddy.my(tt),1);
        ze  = squeeze(zr(imx,imy,:)); % z co-ordinate at center of eddy
        try
            eddy.T(tt,:)   = double(squeeze(ncread(file,'temp',[imx imy 1 tt-tt0],[1 1 Inf 1])));
        catch ME
            disp([' Moving to next file tt = ' num2str(tt) ' - ' char(fnames(kk))]);
            file = [dir1 '/' char(fnames(kk))];
            kk = kk +1;
            tt0 = tt-1;
            eddy.T(tt,:)   = double(squeeze(ncread(file,'temp',[imx imy 1 tt-tt0],[ ...
                                                    1 1 Inf 1])));
        end
        
        if params.bathy.axis == 'x'
            Ti = double(squeeze(ncread(file,'temp',[imx  size(xr,2)  1 tt-tt0],[1 1 Inf 1])));
        else
            Ti = double(squeeze(ncread(file,'temp',[size(xr,1)  imy  1 tt-tt0],[1 1 Inf 1])));
        end
        opts = optimset('MaxFunEvals',1e3);
        if ~isfield(params.flags,'vprof_gaussian') || params.flags.vprof_gaussian
            [x2,~,exitflag] = fminsearch(@(x) gaussfit2(x,eddy.T(tt,:)'-Ti,ze),initGuess2,opts);
            if ~exitflag, x2(2) = NaN; end
            %[x3,~,exitflag] = fminsearch(@(x) gaussfit3(x,eddy.T(tt,:)'-Ti,ze),initGuess3,opts);
            %if ~exitflag, x3(2) = NaN; end
            eddy.Lz2(tt)  = abs(x2(2));
            %eddy.Lz3(tt)  = NaN;%abs(x3(2));
        else
            %fit sinusoid
            [x2,~,exitflag] = fminsearch(@(x) sinefit(x,eddy.T(tt,:)'-Ti,ze),initGuess,opts);
            eddy.Lz2(tt) = abs(2*pi/x2(2));
            %eddy.Lz3(tt) = NaN;
        end
        
        % pcolor(xr,yr,eddy.mask(:,:,tt).*zeta(:,:,tt)); linex(eddy.mx(tt));title(num2str(tt))
        % calculate center velocity
        if tt == 1
            eddy.mvx(1) = NaN;
            eddy.mvy(1) = NaN;
            eddy.cvx(1) = NaN;
            eddy.cvy(1) = NaN;
        else
            % dt in days; convert dx,dy to km
            eddy.mvx(tt) = (eddy.mx(tt) - eddy.mx(tt-1))./dt/1000;
            eddy.mvy(tt) = (eddy.my(tt) - eddy.my(tt-1))./dt/1000;
            
            % dt in days; convert dx,dy to km
            eddy.cvx(tt) = (eddy.cx(tt) - eddy.cx(tt-1))./dt/1000;
            eddy.cvy(tt) = (eddy.cy(tt) - eddy.cy(tt-1))./dt/1000;
        end
    end
    eddy.xr = xr;
    eddy.yr = yr;
    
    eddy.Lz2(abs(eddy.Lz2) > max(abs(zr(:)))) = NaN;
    %eddy.Lz3(abs(eddy.Lz3) > max(abs(zr(:)))) = NaN;
    toc;
    
    eddy.comment = ['(cx,cy) = Location of weighted center (m) | ' ...
                    '(mx,my) = Location of SSH max closest to shelfbreak (m) | ' ...
                    'amp = amplitude (m) | dia = diameter of circle with same area (m) | ' ...
                    'mask = SSH mask to check detection | n = number of pixels | ' ...
                    '(mvx,mvy) = velocity of (mx,my) in km/day | ' ...
                    '(we,ee,ne,se) = West, East, North and South edges of eddy (m) | ' ...
                    'Lz2,3 = Vertical scale (m) when fitting without & with linear trend | ' ...
                    'T = temp profile at (mx,my) | '...
                    'Lmin/Lmaj = minor/major axis length'];
    
    save([dir1 '/eddytrack.mat'],'eddy');
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
function [eddy] = eddy_diag(zeta,vor,dx,dy,sbreak,w)

    % algorithm options
    amp_thresh = 0.001; % Amplitude threshold (in m)
    
    % minimum eddy rad. = 5 km, maximum = 100 km
    low_n  = floor(pi*(5e3)^2/dx/dy);       % minimum number of pixels in eddy
    high_n = floor(pi*(100e3)^2/dx/dy);     % maximum number of pixels in eddy
    connectivity = 8;  % either 4 or 8
    max_dist = 400*1000;
    thresh_loop = nanmin(zeta(:)):0.005:nanmax(zeta(:)); % in m

    for ii=1:length(thresh_loop)
        threshold = thresh_loop(ii);
        
        % Criterion 1 - first get everything above threshold
        mask = zeta > threshold;
        %(mask');
        
        % first find simply connected regions
        regions = bwconncomp(mask,connectivity);
        
        % sort regions by n
        nn = [];
        for ii=1:regions.NumObjects
            nn(ii) = length(regions.PixelIdxList{ii});
        end
        [~,ind] = sort(nn,'descend');
        
        % loop over these regions and apply the criteria
        for kk=1:regions.NumObjects
            jj = ind(kk);
            % Criterion 2 - either too big or too small
            n = length(regions.PixelIdxList{jj});
            if n < low_n || n > high_n, continue; end

            % reconstruct zeta for identified region
            maskreg = zeros(regions.ImageSize);
            maskreg(regions.PixelIdxList{jj}) = 1;
            zreg   = repnan(zeta .* maskreg,0); % zeta in region
            zperim = zeta .* bwmorph(maskreg,'remove'); % zeta in perimeter
            
            % Criterion 3 - need local maximum in the 
            % see http://stackoverflow.com/questions/1856197/how-can-i-find-local-maxima-in-an-image-in-matlab
            local_max = imregionalmax(zreg);
            if local_max == zeros(size(zreg)), continue; end
            %if sum(local_max(:)) > 1, continue; end % skip if more than one maximum
            
            props = regionprops(maskreg,zeta,'WeightedCentroid','Solidity', ...
                        'EquivDiameter','Area','BoundingBox', ...
                        'MinorAxisLength','MajorAxisLength');
            
            % Criterion 4 - amplitude is at least > amp_thresh
            indices = find(local_max == 1);
            if length(indices) > 20; continue; end
            % make sure all local maxima found satisfy this criterion
            for mm=1:length(indices)
                if (zreg(indices(mm)) - nanmean(zperim(:))) < amp_thresh
                   local_max(indices(mm)) = 0;
                end
            end
            amp = max(zreg(local_max)) - nanmean(zperim(:));
            if isempty(amp) || (amp < amp_thresh) || isequal(local_max,zeros(size(local_max)))
                continue; 
            end
            
            % Criterion 5 - distance between any pair of points must be
            % less than a specified maximum            
            % need to trace out outline, first find one point to start from
            % (convex hull doesn't do the job)
            [x0,y0] = ind2sub(size(maskreg),find(maskreg == 1,1,'first'));
            points = bwtraceboundary(maskreg,[x0,y0],'E');
            points = bsxfun(@times,points,[dx dy]);
            if maximumCaliperDiameter(points) > max_dist, continue; end
           % if minimumCaliperDiameter(points) < min_dist, continue; end
            
            % Criterion 6 - Obuko-Weiss parameter must make sense
            if exist('w','var')
                thresh = 0.2 * std(w(:)); % see Isern-Fontanet et al. (2006)
                wreg = w .* maskreg;
                if mean(wreg(:)) > thresh
                    break
                end
            end
            
            % Criterion 7 - solidity must be good - helps get rid of some
            % thin 'isthumuses'
            if props.Solidity  < 0.90, continue; end
            
            % Criterion 8 - low 'rectangularity' - gets rid of rectangles
            % that are sometime picked up
            rectarea = props.BoundingBox(3) * props.BoundingBox(4);
            if props.Area/rectarea > 0.85, continue; end
            
            % Calculate properties of region
            %c = regionprops(maskreg,zeta,'WeightedCentroid');
            
            % Criterion 7 - if multiple regions (eddies), only store the one
            % closest to shelfbreak
            cx = props.WeightedCentroid(2) * dx;
            cy = props.WeightedCentroid(1) * dy;
            
            if exist('eddy','var') 
                if (cx-sbreak) > (eddy.cx-sbreak) && (eddy.cx-sbreak > 0)
                    continue
                else
                    % current eddy is closer but could be bigger than we want
                    % mask out the eddy that is farther away and
                    % recursively call this function
%                     mask = ones(regions.ImageSize);
%                     mask(regions.PixelIdxList{eddy.jj}) = 0; 
%                     eddy = eddy_diag(zeta .* mask,dx,dy,sbreak,w);
%                     return
                    
                end
            end

            % I have an eddy!!!
            flag_found = 1;
            fprintf('Eddy found with threshold %.3f \n', threshold);
            %imagesc(zreg'); pause
            
            % find location of maximum that is closest to shelfbreak
            ix = bsxfun(@times,ones(size(local_max)),[1:size(local_max,1)]');
            iy = bsxfun(@times,ones(size(local_max)),[1:size(local_max,2)]);
            indx = nanmin(nanmin(fillnan(local_max.*ix,0)));
            [~,indy] = max(local_max(indx,:));
            
            xmax = fillnan(maskreg(:).*ix(:),0);
            ymax = fillnan(maskreg(:).*iy(:),0);
            
            % store eddy properties for return
            eddy.cx   = cx; % weighted center
            eddy.cy   = cy; %     "
            eddy.mx   = ix(indx,indy) * dx; % maximum
            eddy.my   = iy(indx,indy) * dy; %    "
            eddy.we   = nanmin(xmax) * dx; % west edge
            eddy.ee   = nanmax(xmax) * dx; % east edge
            eddy.ne   = nanmax(ymax) * dy; % south edge
            eddy.se   = nanmin(ymax) * dy; % north edge
            eddy.amp  = amp;
            eddy.dia  = props.EquivDiameter * sqrt(dx*dy);
            eddy.mask = maskreg;
            eddy.n    = n; % number of points
            eddy.jj   = jj;
            
            % find 0 rel. vor (max. speed) contour
            vormask   = vor.*eddy.mask < 0;
            vorregions = bwconncomp(vormask,connectivity);
            for zz=1:vorregions.NumObjects
                nvor(zz) = length(vorregions.PixelIdxList{zz});
            end
            [~,ind] = sort(nvor,'descend');
            
            % ugh. now I have to pick the right region, like earlier
            for mm=1:vorregions.NumObjects
                ll = ind(mm);
                vormaskreg = zeros(vorregions.ImageSize);
                vormaskreg(vorregions.PixelIdxList{ll}) = 1;
                
                % only 1 criterion for now
                n = length(vorregions.PixelIdxList{ll});
                if n < low_n || n > high_n, continue; end
                
                vorprops  = regionprops(vormaskreg,zeta,'EquivDiameter', ...
                        'MinorAxisLength','MajorAxisLength');
                eddy.L    = vorprops.EquivDiameter * sqrt(dx*dy);
                eddy.lmaj = vorprops.MajorAxisLength * sqrt(dx*dy);
                eddy.lmin = vorprops.MinorAxisLength * sqrt(dx*dy);
                eddy.vormask = vormaskreg;
                % if i get here, i'm done
                break;
            end
            
%             imagesc(zreg');
%             linex(cx./dx); liney(cy./dy);
%             linex(indx,[],'r'); liney(indy,[],'r');
%             linex(nanmin(xmax),[],'b'); linex(nanmax(xmax),[],'b');
%             liney(nanmin(ymax),[],'b'); liney(nanmax(ymax),[],'b');
            
            % stop when eddy is found
            break;
        end
        if exist('flag_found','var') && flag_found == 1, break; end
    end 
    
    if ~exist('eddy','var')
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
        disp('Eddy not found!');
    end
    
    try
        rmfield(eddy,'jj');
    catch ME
    end
    
    
%             % first find convex hull vertices
%             chull = regionprops(maskreg,'ConvexHull');
%             chull = chull.ConvexHull;
%             chull = bsxfun(@times,chull,[dx dy]);
