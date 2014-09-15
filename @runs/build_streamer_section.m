% extract points for streamer section
function [] = build_streamer_section(runs)

% make plots to check?
    debug_plot = 1;

    runs.streamer.west.fit_circle = 1;

    if ~isfield(runs.streamer,'yend')
        runs.detect_streamer_mask();
    end
    yend = runs.streamer.yend;
    xr = runs.rgrid.xr(:,1:yend)/1000;
    yr = runs.rgrid.yr(:,1:yend)/1000;

    cx = runs.eddy.mx/1000;
    cy = runs.eddy.my/1000;
    cy(cy > max(yr(:))) = max(yr(:));

    cxind = vecfind(xr(:,1),runs.eddy.mx/1000);
    cyind = vecfind(yr(1,:),cy)';

    for tind=1:size(runs.streamer.west.mask,2)
        % now pick ONLY SURFACE
        stream = reshape(full(runs.streamer.west.mask(:,tind)), ...
                         runs.streamer.sz3dfull);
        stream = stream(:,:,end); % SURFACE ONLY

        % if no streamer or too small, skip
        if isequal(stream,zeros(size(stream))) ...
                || numel(find(stream(:) ~= 0)) < 150
            continue;
        end

        % code from
        % http://blogs.mathworks.com/steve/2014/01/07/automating-data-extraction-2/x
        skeleton = bwmorph(stream,'skel','inf');
        skel = breakapart(skeleton);
        skelcomps = bwconncomp(skel);
        % find distance from eddy center?
        distcen = sqrt( (skel.*xr - runs.eddy.cx(tind)).^2 +  ...
                        (skel.*yr - runs.eddy.cy(tind)).^2 );
        distcen = distcen .* fillnan(skel,0);
        meandist = nan([skelcomps.NumObjects 1]);

        icen = nan(skelcomps.NumObjects,2);

        % process the branches for mean distance, centroid, and sort
        % clockwise
        for mm = 1:skelcomps.NumObjects
            meandist(mm) = nanmean(distcen(skelcomps.PixelIdxList{mm}));

            [ixtemp,iytemp] = ind2sub(size(skel), ...
                                      skelcomps.PixelIdxList{mm});
            % don't remap to preserve order of points crossing the
            % horizontal axis
            %[~,sorttang] = angleSort([ixtemp iytemp], ...
            %                [cxind(tind) cyind(tind)],-pi/2);
            %sorttang = flipdim(sorttang,1);
            % works with 0 crossing
            tempang = atan2(iytemp-cyind(tind),ixtemp-cxind(tind));

            if max(diff(tempang) > 5.9)
                tempang = mod(tempang + 2*pi,2*pi);
            end
            %tempang(tempang < 0) = tempang(tempang < 0) + 360;
            %[~,sorttang] = sort(tempang,'descend');
            %skelcomps.PixelIdxList{mm} = skelcomps.PixelIdxList{mm}(sorttang);
            %if ~isclockwise(ixtemp,iytemp)

            % sort points in each branch clockwise. This is imposed by
            % setting the reference angle (w.r.t eddy center) to be the
            % minimum of all point angles in the branch
            refAngle = min(tempang(:));
            [out,~] = angleSort([ixtemp iytemp], ...
                                [cxind(tind) cyind(tind)],refAngle);
            skelcomps.PixelIdxList{mm} = sub2ind(size(skel), ...
                                                 flipud(out(:,1)),flipud(out(:,2)));
            %testBranch(skelcomps.PixelIdxList{mm},size(skel));

            % store centroid and find it's angle w.r.t eddy center
            icen(mm,:) = centroid([ixtemp(:) iytemp(:)]);
        end

        % sort by distance
        [~,sortdist] = sort(meandist);
        %, then chuck top 20%
        %indices = cat(1, ...
        %    skelcomps.PixelIdxList{ sortdist(1: floor(0.8*length(sortdist)) ) });
        %indices = skelcomps.PixelIdxList{sortdist(1)};

        % measure number of pixels in each branch and
        % throw out small branches
        numPixels = cellfun(@numel,skelcomps.PixelIdxList);
        numPixels(numPixels < 5) = NaN;
        [~,sortnum] = sort(numPixels);
        nanindices = cut_nan(fillnan(isnan(numPixels) ...
                                     .* (1:skelcomps.NumObjects),0));
        for mm=1:length(nanindices)
            sortnum(sortnum == nanindices(mm)) = NaN;
        end

        % remove farthest away segment for sure
        if skelcomps.NumObjects > 1
            sortnum( sortnum == sortdist(end) ) = NaN;
        end
        % chuck out indices I'm not interested in
        sortnum = cut_nan(sortnum);

        % if region is too small, exit
        if isempty(sortnum)
            warning(['skipping @ tt=' num2str(tind)]);
            continue;
        end

        if runs.streamer.west.fit_circle
            % first get discrete points
            % old version without joining
            indices = cat(1,skelcomps.PixelIdxList{sortnum});
            [ixstr,iystr] = ind2sub(size(skel),indices);
            xstr = xr(ixstr,1);
            ystr = yr(1,iystr)';

            % fit circle
            circ = CircleFitByPratt([xstr ystr]);
            Cx = circ(1); Cy = circ(2); R = circ(3);
            theta0 = unwrap(atan2(ystr-Cy,xstr-Cx));
            % i want 2 km resolution i.e., R * dtheta = 2 km
            dtheta = 2/R;
            theta = min(theta0(:)):dtheta:max(theta0(:));
            xstr = Cx + R .* cos(theta);
            ystr = Cy + R .* sin(theta);

            strmask = round(interp2(xr',yr',stream',xstr,ystr));
            xstr(strmask == 0) = [];
            ystr(strmask == 0) = [];

            if ~isclockwise(xstr,ystr)
                xstr = fliplr(xstr)';
                ystr = fliplr(ystr)';
            end

            ixstr = []; iystr = [];
        else
            % use angleSort on branch centroids to order regions appropriately
            [~,sortcen] = angleSort(icen(sortnum,:), ...
                                    [cxind(tind) cyind(tind)],-pi/2);
            sortcen = flipdim(sortcen,1);
            % alternative to above - sortcen code
            %[~,sortang] = sort(meanangle(sortnum),'descend');

            % now actually select the remaining regions and figure out
            % (x,y) co-ordinates
            sortnum = sortnum(sortcen);

            % sortnum should be final sorted order here
            [ixstr, iystr] = ind2sub(size(skel), skelcomps.PixelIdxList{sortnum(1)});
            for mm = 1:length(sortnum)-1
                ix1 = ixstr(end);
                iy1 = iystr(end);

                [ix2,iy2] = ind2sub(size(skel), ...
                                    skelcomps.PixelIdxList{sortnum(mm+1)});

                % use Bresenham's algorithm to join
                [jx,jy] = bresenham(ix1,iy1,ix2(1),iy2(1));

                ixstr = [ixstr; jx; ix2];
                iystr = [iystr; jy; iy2];
            end

            xstr = xr(ixstr,1);
            ystr = yr(1,iystr)';
        end

        % fix the starting!!!
        dstr = [0; cumsum(hypot(diff(xstr),diff(ystr)))];

        % distance from perimeter - NOT QUITE AS GOOD
        %{
        distper = bwdist(~stream);
        [~,index1] = max(distper(1:cxind(tt),:));
        [~,index2] = max(distper(cxind(tt):end,:));
        index1(index1 == 1) = NaN;
        index2(index2 == 1) = NaN;
        index2 = index2+cxind(tt);
        idxx = [index1(:); fliplr(index2(:))]';
        idxy = [1:size(stream,2) fliplr(1:size(stream,2))];
        %}
        %polyline = [cut_nan(idxx)' (cut_nan(idxy .* idxx)./cut_nan(idxx))'];

        % testing streamer cross-section detection
        if debug_plot
            clf
            subplot(211)
            pcolorcen(xr,yr,double(stream));
            hold on;
            plot(cx(tind),cy(tind),'ko','MarkerSize',16);
            %plot(idxx,idxy,'bx','markersize',8);
            plot(xstr,ystr,'k*');

            drawCircle(circ(1),circ(2),circ(3));
            %ell = EllipseDirectFit([ixstr iystr]);
            %a = ell(1); b = ell(2); c = ell(3);
            %d = ell(4); e = ell(5); f = ell(6);
            %x0 = (c*d - b*f)/(b^2-a*c);
            %y0 = (a*f - b*d)/(b^2-a*c);

            title(num2str(tind));
            subplot(212)
            for mm=1:length(sortnum)
                testBranch(skelcomps.PixelIdxList{sortnum(mm)}, ...
                           size(skel),sortnum(mm));
            end
            set(gca,'ydir','normal');
            hold on;
            plot(cxind(tind),cyind(tind),'ko','MarkerSize',16);
            plot(icen(:,1),icen(:,2),'k*','MarkerSize',8);
            pause();
        end

        % save locations in runs object
        runs.streamer.west.xstr{tind}  = xstr;
        runs.streamer.west.ystr{tind}  = ystr;
        runs.streamer.west.ixstr{tind} = ixstr;
        runs.streamer.west.iystr{tind} = iystr;
        runs.streamer.west.dstr{tind}  = dstr;
        runs.streamer.comment   = ...
            [' contour = 1 in streamer, 0 outside |\n ' ...
             ' (xstr,ystr) = cross-section through streamer (cell array) |\n ' ...
             ' (ixstr,iystr) = indices corresponding to (xstr,ystr) ' ...
             ' - (cell array) |\n dstr = along-streamer distance (cell array)'];
    end

    % save to file
    disp('Writing to file');tic;
    streamer = runs.streamer;
    streamer.hash = githash([mfilename('fullpath') '.m']);
    save([runs.dir '/streamer.mat'],'streamer');
    toc;
end
