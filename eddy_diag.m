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
    connectivity = 4;  % either 4 or 8
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
            amp = max(zreg(:)) - mean(zperim(:));
            if amp < amp_thresh, continue; end
            
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
                if (cx-sbreak) > (eddy.cx-sbreak)
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
            
            % store eddy properties for return
            eddy.cx   = cx; % weighted center
            eddy.cy   = cy; %     "
            eddy.mx   = ix(indx,indy) * dx; % maximum
            eddy.my   = iy(indx,indy) * dy; % "
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
        eddy.amp  = NaN;
        eddy.dia  = NaN;
        eddy.mask = NaN;
        eddy.n    = NaN;
    end
    
    try
        rmfield(eddy,'jj');
    catch ME
    end
