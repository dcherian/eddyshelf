% Calculates eddy diagnostics as in Chelton et al. (2011)
% doesn't support multiple eddies yet
% only finds anticyclonic eddies

function [eddy] = eddy_diag(zeta,dx,dy)

    zeta = zeta - min(zeta(:));
    
    % algorithm options
    amp_thresh = 0.01; % Amplitude threshold (in m)
    thresh_loop = linspace(amp_thresh,max(zeta(:)),100); % in m
    low_n  = 15;       % minimum number of pixels in eddy
    high_n = 1300;     % maximum number of pixels in eddy
    connectivity = 8;  % either 4 or 8
    max_dist = 400*1000;

    for ii=1:length(thresh_loop)
        threshold = thresh_loop(ii);
        
        % Criterion 1 - first get everything above threshold
        mask = zeta > threshold;
        %imagesc(mask');
        
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
            zreg   = zeta .* maskreg; % zeta in region
            zperim = zeta .* bwmorph(maskreg,'remove'); % zeta in perimeter
            
            % Criterion 3 - need local maximum in the 
            % see http://stackoverflow.com/questions/1856197/how-can-i-find-local-maxima-in-an-image-in-matlab
            if imregionalmax(zreg) == zeros(size(zreg)), continue; end

            % Criterion 4 - amplitude is at least > amp_thresh
            if min(fillnan(zreg(:),0)) < amp_thresh, continue; end
            
            % Criterion 5 - distance between any pair of points must be
            % less than a specified maximum
            % first find convex hull vertices
            chull = regionprops(maskreg,'ConvexHull');
            chull = chull.ConvexHull;
            chull = bsxfun(@times,chull,[dx dy]);
            if maximumCaliperDiameter(chull) > max_dist, continue; end
            
            % I have an eddy!!!
            fprintf('Eddy found with threshold %.2f \n', threshold);
            
            flag_found = 1;
            
            % Calculate properties
            c = regionprops(maskreg,zeta,'WeightedCentroid');
            eddy.cx   = c.WeightedCentroid(2) * dx;
            eddy.cy   = c.WeightedCentroid(1) * dy;
            eddy.amp  = max(zreg(:)) - mean(zperim(:));
            eddy.dia  = regionprops(maskreg,'EquivDiameter');
            eddy.dia  = eddy.dia.EquivDiameter * sqrt(dx*dy);
            eddy.mask = maskreg;
            
            % stop when eddy is found
            break;
        end
        if exist('flag_found','var') && flag_found == 1, break; end
    end 
    
    if ~exist('eddy','var')
        eddy.cx   = NaN;
        eddy.cy   = NaN;
        eddy.amp  = NaN;
        eddy.dia  = NaN;
        eddy.mask = NaN;
    end
