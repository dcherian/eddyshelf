function [eddy] = track_eddy(dir)

    fnames = ls([dir '/*his*.nc']);
    if isempty(ls([dir '/*his*.nc']))
        fnames = ls([dir '/*avg*.nc']);
    end

    fname = [dir '/' fnames(1,:)];
    % search region for tracking eddies (in addition to detected diameter)
    limit_x = 40*1000;
    limit_y = 40*1000;

    zeta = roms_read_data(dir,'zeta');
    [xr,yr,~,~,~,~] = roms_var_grid(fname,'zeta');
    eddy.h = ncread(fname,'h');
    eddy.t = roms_read_data(dir,'ocean_time'); % required only for dt
    dt = eddy.t(2)-eddy.t(1);

    zeta = zeta(2:end-1,2:end-1,:);
    xr = xr(2:end-1,2:end-1);
    yr = yr(2:end-1,2:end-1);

    dx = xr(2,1)-xr(1,1);
    dy = yr(1,2)-yr(1,1);

    sz = size(zeta(:,:,1));

    % detect shelfbreak
    sbreak = find_shelfbreak(fname);

    % remove background flow contribution to zeta
    zeta_bg = zeta(:,end,1);

    tic;
    for tt=1:size(zeta,3)
%         if tt == 126
%             disp('time to debug');
%         end
        if tt == 1,
            mask = ones(sz);
            d_sbreak = Inf;
        else 
            mask = nan*ones(sz);
            lx = eddy.dia(tt-1)/2 + limit_x;
            ly = eddy.dia(tt-1)/2 + limit_y;
            ix1 = find_approx(xr(:,1),eddy.cx(tt-1)-lx,1);
            ix2 = find_approx(xr(:,1),eddy.cx(tt-1)+lx,1);
            iy1 = find_approx(yr(1,:),eddy.cy(tt-1)-ly,1);
            iy2 = find_approx(yr(1,:),eddy.cy(tt-1)+ly,1);

            mask(ix1:ix2,iy1:iy2) = 1;

            d_sbreak = eddy.cx-sbreak;
        end
        temp = eddy_diag(bsxfun(@minus,zeta(:,:,tt),zeta_bg) .* mask, ...
                            dx,dy,sbreak); %w(:,:,tt));

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
        % number of pixels in eddy
        eddy.n(tt) = temp.n;
        % north, south, west & east edges
        eddy.we(tt) = temp.we;
        eddy.ee(tt) = temp.ee;
        eddy.ne(tt) = temp.ne;
        eddy.se(tt) = temp.se;
        
        % pcolor(xr,yr,eddy.mask(:,:,tt).*zeta(:,:,tt)); linex(eddy.mx(tt));title(num2str(tt))
        % calculate center velocity
        if tt == 1
            eddy.mvx(1) = NaN;
            eddy.mvy(1) = NaN;
        else
            eddy.mvx(tt) = (eddy.mx(tt) - eddy.mx(tt-1))./dt;
            eddy.mvy(tt) = (eddy.my(tt) - eddy.my(tt-1))./dt;
        end
    end
    eddy.xr = xr;
    eddy.yr = yr;
    toc;
    
    save([dir '/eddytrack.mat'],'eddy');
    disp('Done.');

% Calculates eddy diagnostics as in Chelton et al. (2011)
% doesn't support multiple eddies yet
% only finds anticyclonic eddies

% this routine is called at every timestep
function [eddy] = eddy_diag(zeta,dx,dy,sbreak,w)

    % algorithm options
    amp_thresh = 0.01; % Amplitude threshold (in m)
    thresh_loop = nanmin(zeta(:)):0.01:nanmax(zeta(:)); % in m
    low_n  = 400;       % minimum number of pixels in eddy
    high_n = 2500;     % maximum number of pixels in eddy
    connectivity = 8;  % either 4 or 8
    max_dist = 400*1000;

    for ii=1:length(thresh_loop)
        threshold = thresh_loop(ii);
        
        % Criterion 1 - first get everything above threshold
        mask = zeta > threshold;
        %(mask');
        
        % first find simply connected regions
        regions = bwconncomp(mask,connectivity);
        
        % loop over these regions and apply the criteria
        for jj=1:regions.NumObjects
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
            
            % Criterion 4 - amplitude is at least > amp_thresh
            indices = find(local_max == 1);
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
            % first find convex hull vertices
            chull = regionprops(maskreg,'ConvexHull');
            chull = chull.ConvexHull;
            chull = bsxfun(@times,chull,[dx dy]);
            if maximumCaliperDiameter(chull) > max_dist, continue; end
            
            % Criterion 6 - Obuko-Weiss parameter must make sense
            if exist('w','var')
                thresh = 0.2 * std(w(:)); % see Isern-Fontanet et al. (2006)
                wreg = w .* maskreg;
                if mean(wreg(:)) > thresh
                    break
                end
            end
            
            % Calculate properties of region
            c = regionprops(maskreg,zeta,'WeightedCentroid');
            
            % Criterion 7 - if multiple regions (eddies), only store the one
            % closest to shelfbreak
            cx = c.WeightedCentroid(2) * dx;
            cy = c.WeightedCentroid(1) * dy;

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
            % ASSUMES ISOBATHS ARE NORTH-SOUTH
            ix = bsxfun(@times,ones(size(local_max)),[1:size(local_max,1)]');
            iy = bsxfun(@times,ones(size(local_max)),[1:size(local_max,2)]);
            indx = nanmin(nanmin(fillnan(local_max.*ix,0)));
            [~,indy] = max(local_max(indx,:));
            
            xmax = ix(logical(maskreg(:,indy)),1);
            ymax = iy(1,logical(maskreg(indx,:)));
            
            % store eddy properties for return
            eddy.cx   = cx; % weighted center
            eddy.cy   = cy; %     "
            eddy.mx   = ix(indx,indy) * dx; % maximum
            eddy.my   = iy(indx,indy) * dy; %    "
            eddy.we   = min(xmax) * dx; % west edge
            eddy.ee   = max(xmax) * dx; % east edge
            eddy.ne   = min(ymax) * dy; % south edge
            eddy.se   = max(ymax) * dy; % north edge
            eddy.amp  = amp;
            eddy.dia  = regionprops(maskreg,'EquivDiameter');
            eddy.dia  = eddy.dia.EquivDiameter * sqrt(dx*dy);
            eddy.mask = maskreg;
            eddy.n    = n; % number of points
            eddy.jj   = jj;
            
%             imagesc(zreg');
%             linex(cx./dx); liney(cy./dy);
%             linex(indx,[],'r'); liney(indy,[],'r');
            
            % stop when eddy is found
            %break;
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
        disp('Eddy not found!');
    end
    
    try
        rmfield(eddy,'jj');
    catch ME
    end
