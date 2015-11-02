% detect streamer contours
function [] = detect_streamer_mask(runs)

% upper y-limit to save memory
    yend = find_approx(runs.rgrid.y_rho(:,1),130*1000);
    t0 = runs.eddy.trevind;
    %read_start = [1 1 1 t0-20];
    slab = 40;

    runs.streamer.yend = yend;

    szeta = size(runs.zeta);

    % allocate memory
    nanvec = nan(size(runs.time));
    runs.streamer.time = nanvec;
    runs.streamer.off.vol = nanvec;
    runs.streamer.off.zcen = nanvec;
    runs.streamer.off.zdcen = nanvec;

    % grid matrices required for plotting
    xsb = runs.bathy.xsb/1000;
    xr = runs.rgrid.xr(:,1:yend)/1000;
    yr = runs.rgrid.yr(:,1:yend)/1000;
    zr = permute(runs.rgrid.z_r(:,1:yend,:),[3 2 1]);
    %ix = repmat([1:size(xr,1)]',[1 yend]);
    %iy = repmat([1:yend],[size(xr,1) 1]);

    runs.streamer.xr = xr;
    runs.streamer.yr = yr;
    runs.streamer.zr = zr;

    % size matrices to make processing easier
    runs.streamer.sz4dfull = [size(zr) szeta(3)];
    runs.streamer.sz4dsp = [numel(zr) szeta(3)];
    runs.streamer.sz3dsp = [numel(zr) 1];
    runs.streamer.sz3dfull = size(zr);

    % allocate streamer mask variable
    runs.streamer.off.mask = sparse(runs.streamer.sz4dsp(1),szeta(3));

    % grid cell volume
    dVs = reshape(runs.rgrid.dV(:,1:runs.streamer.yend,:), ...
                  runs.streamer.sz3dsp);

    for ii=1:floor(szeta(3)/slab)

        tstart = t0+slab*(ii-1);
        tend = tstart+slab-1;
        if tend > szeta(3), tend = szeta(3); end

        sz4dfull = [runs.streamer.sz4dfull(1:3) tend-tstart+1];
        sz4dsp = [runs.streamer.sz4dsp(1) tend-tstart+1];
        sz3dsp = runs.streamer.sz3dsp;
        sz3dfull = runs.streamer.sz3dfull;


        runs.streamer.time(tstart:tend) = runs.time(tstart:tend);
        tindices = [tstart tend];

        csdye = dc_roms_read_data(runs.dir, runs.csdname, tindices, ...
                                  {'y' 1 yend},[],runs.rgrid)/1000;
        zdye  = dc_roms_read_data(runs.dir, runs.zdname, tindices, ...
                                  {'y' 1 yend},[],runs.rgrid);
        eddye = dc_roms_read_data(runs.dir,runs.eddname, tindices, ...
                                  {'y' 1 yend},[],runs.rgrid);
        %asdye = dc_roms_read_data(runs.dir, runs.asdname, tindices, ...
        %            {'y' 1 yend});

        % identify streamer with 4D data
        % preliminary detection
        % I use cross-shore label to roughly filter first

        % eliminate this step somehow and remove the temporary array?
        streamer1 = (csdye > xsb-10) & (csdye < xsb+30) & (eddye < 0.2);

        % (xs,ys,zs) are the Eulerian x,y,z values
        %xs = bsxfun(@times, streamer, grd.xax)/1000;
        %ys = bsxfun(@times, streamer1, yr);
        %zs = bsxfun(@times, streamer, grd.zax);

        % (as,cs,z) dyes contain the Lagrangian labels
        % some distance metric between the two will give me an idea of
        % what's happening
        %if runs.bathy.axis == 'y'
        %    das = asdye - xs;
        %    dcs = csdye - ys;
        %else
        %    das = asdye - ys;
        %    dcs = csdye - xs;
        %end
        %dz = zdye - zs;

        warning('DO I NEED TO ACCOUNT FOR TILTING IN VERTICAL?');
        %cx = runs.eddy.cx(tstart:tend)/1000;
        %cy = runs.eddy.cy(tstart:tend)/1000;
        ee = runs.eddy.ee(tstart:tend)/1000;
        % hack if eddy center is outside extracted domain
        %cy(cy > max(yr(:))) = max(yr(:));
        %cxind = vecfind(xr(:,1),cx);
        %cyind = vecfind(yr(1,:),cy)';

        % pick out western streamer by chucking points that are >
        % eastern edge + initial radius. This allows streamer to wrap
        % around eddy and not be chopped off.
        west_mask = bsxfun(@lt, repmat(xr,[1 1 runs.rgrid.N]), ...
                           permute(ee + runs.params.eddy.dia/2000,[3 4 1 2]));

        % mask of points west of eddy center - OLD ATTEMPT
        %west_mask = bsxfun(@lt,repmat(runs.rgrid.x_rho',[1 1 runs.rgrid.N]), ...
        %               permute(runs.eddy.cx(runs.eddy.trevind:runs.eddy.trevind+19), [1 3 4 2]));


        %r = sqrt(bsxfun(@minus,xr,permute(cx,[3 1 2])).^2 ...
        %       + bsxfun(@minus,yr,permute(cy,[3 1 2])).^2);

        % picking only western streamer
        streamer1 = squeeze(streamer1  ... % original streamer
                            ... % parcels have moved more than 5 km
        ... %   in the cross-shelf dirn.
        .* (abs(csdye - bsxfun(@times,streamer1,yr))>5)) ...
            ... % remove eastern half
        .* (west_mask);
        %...     % streamer depth is not total depth
        %.* squeeze(bsxfun(@lt,max(abs(zs),[],3), runs.rgrid.h(1:yend,:)'));

        % pick out biggest surface piece
        % it looks like the surface is the biggest so we look for
        % only look for everything under it - i.e., hopefully no
        % tilting

        stream = streamer1(:,:,runs.rgrid.N,:);
        for tt=1:size(stream,4)
            % get biggest part - assume it's what i'm interested in
            strtemp = stream(:,:,1,tt);
            strcomps = bwconncomp(strtemp);
            numPixels = cellfun(@numel,strcomps.PixelIdxList);
            [~,bigidx] = max(numPixels);
            strtemp(strcomps.PixelIdxList{bigidx}) = 2;
            strtemp(strtemp < 2) = 0;
            strtemp(strtemp == 2) = 1;
            stream(:,:,1,tt) = strtemp;
        end

        % filter and save
        runs.streamer.off.mask(:,tstart:tend) = sparse(reshape( ...
            bsxfun(@times,streamer1,stream), ...
            sz4dsp));
        %clear west_mask streamer1 stream strtemp;

        % compress somehow
        %streamnan = fillnan(runs.streamer.off.mask,0);
        % calculate statistics
        %xs = bsxfun(@times, streamnan, xr);
        %ys = bsxfun(@times, streamnan, yr);
        zs = bsxfun(@times, runs.streamer.off.mask(:,tstart:tend), ...
                    reshape(zr,sz3dsp));

        zdyestr = runs.streamer.off.mask(:,tstart:tend) .* ...
                  reshape(zdye,sz4dsp);
        %csdyestr = bsxfun(@times, streamnan, csdye);

        %dcs  = abs(csdyestr - ys);
        %dzd  = abs(zdyestr - zs);

        % calculate volume
        runs.streamer.off.vol(tstart:tend) = runs.domain_integratesp( ...
            runs.streamer.off.mask(:,tstart:tend), dVs);

        % Haven't used temperature yet

        % suffix cen = just centroids
        % suffix dcen = centroid weighted by dye value
        runs.streamer.off.zcen(tstart:tend) = bsxfun(@rdivide, ...
                                                      runs.domain_integratesp(zs,dVs), ...
                                                      runs.streamer.off.vol(tstart:tend));
        runs.streamer.off.zdcen(tstart:tend) = bsxfun(@rdivide,...
                                                       runs.domain_integratesp(zdyestr,dVs), ...
                                                       runs.streamer.off.vol(tstart:tend));

        % volume v/s depth plot for streamer
        % VECTORIZE SOMEHOW
        disp('Binning streamer volume...');
        tic;
        dbin = 20;
        bins = -1*[0:dbin:1000];
        % required so that 0 bin doesn't get a ton of points
        %zsf = fillnan(full(zs),0);
        %sz = size(runs.streamer.off.mask(:,tstart:tend));
        parfor kk=1:length(bins)-1
            temparray(kk,:) = sum(bsxfun(@times, ...
                                         (zs < bins(kk) & zs >= bins(kk+1)), ...
                                         dVs),1);
        end
        runs.streamer.off.Vbin(:,tstart:tend) = temparray;
        runs.streamer.bins = bins;
        toc;
    end
end
