% dir = ['runs/topoeddy/runteb-04-hires-7/'];
% run = runs(dir);
% 
% %% masks 
% 
% try
%     run.dye = squeeze(double(ncread(run.out_file,'dye_02',[1 1 run.rgrid.N 1],[Inf Inf 1 Inf])));
% catch ME
% end

% first bin floats?
function [] = transport_floats(run)
    xvec = run.rgrid.xr(2:end-1,1);
    yvec = run.rgrid.yr(1,2:end-1)';
    sz = size(run.rgrid.xr(2:end-1,2:end-1));
    floats = run.roms;
    eddy = run.eddy;
    bathy = run.bathy;
    
    % rossby radius
    rr = sqrt(run.params.phys.N2)*bathy.hsb/run.rgrid.f(bathy.isb,1);
    distance = 5*rr; % 5 times rossby radius
    
    clim = [bathy.xsb/1000 bathy.xsb/1000+distance/1000]; % colorbar
    disp(['limit = ' num2str(clim)]);
    
    % mask based on initial position
    initmask = ones([1 size(floats.x,2)]);
    % mask if floats have been in eddy before or at present instant
    fineddymask = initmask;
    
    if bathy.axis == 'y'
        initmask(floats.init(:,2) < clim(1)*1000) = NaN;
        initmask(floats.init(:,2) > clim(2)*1000) = NaN;
    else
        initmask(floats.init(:,1) < bathy.xsb) = NaN;
    end
    initmask(floats.init(:,3) < -100) = NaN;
    

    % transform to (r,theta) about eddy center
    % interpolate center location to float times
    cxi = interp1(eddy.t,eddy.cx,floats.time/86400);
    cyi = interp1(eddy.t,eddy.cy,floats.time/86400);
    % is r used elsewhere?
    [th,rr] = cart2pol(bsxfun(@minus,floats.x,cxi), ...
                      bsxfun(@minus,floats.y,cyi));
    th = th .*180/pi;
    dth = diff(th,1,1);
    % find Delta-theta > 350 as location where float crosses center line of
    % eddy
    indices = bsxfun(@times,double(ceil(dth) > 350),[1:size(dth,1)]');
    indices = nanmin(fillnan(indices,0),[],1) + 1; % add one to account for diff earlier
    
    % find change in angle by dtheta after the float has crossed the
    % center line
    dtheta = 90;
    tic;
    for a = 1:size(th,2)
        if ~isnan(indices(a))
            dtt = th(indices(a):end,a)-th(indices(a),a);
            ind = find(abs(dtt) >= dtheta,1,'first') + indices(a);
            floats.x(ind:end,a) = NaN;
            floats.y(ind:end,a) = NaN;
        end
    end
    toc;
    disp([num2str(length(cut_nan(indices))) ' float tracks have been truncated.']);
    
    %% test plot
%     figure
%     a = 324; 
%     dtt = th(indices(a):end,a)-th(indices(a),a);
%     %dtheta = 90;
%     ind = find(abs(dtt) >= dtheta,1,'first') + indices(a);
%     subplot(211)
%     scatter(floats.x(:,a)/1000,floats.y(:,a)/1000,40,th(:,324),'filled'); 
%     hold on; plot(cxi/1000,cyi/1000);
%     colormap(flipud(copper)); colorbar
%     liney(run.bathy.xsb/1000);
%     plot(floats.x(ind,a)/1000,floats.y(ind,a)/1000,'kx','MarkerSize',20);
%     title(['float track ' num2str(a) ' with theta=color']);
%     legend('theta (deg)','eddy center','shelfbreak','Location','Northwest');
%     xlabel('X (km)'); ylabel('Y (km)');
%     beautify([14 14 16]);
%     subplot(212)
%     plot(floats.time/86400,th(:,324));
%     liney([-180 180]);
%     xlabel('time (days)');
%     ylabel('theta (deg)');
%     beautify([14 14 16]);
    
    %%
    if run.makeVideo
%                 runs.mm_instance = mm_setup;
%                 runs.mm_instance.pixelSize = [1600 900];
%                 runs.mm_instance.outputFile = 'mm_output.avi';
%                 runs.mm_instance.ffmpegArgs = '-q:v 1 -g 1';
%                 runs.mm_instance.InputFrameRate = 3;
%                 runs.mm_instance.frameRate = 3;
        aviobj = VideoWriter('ptfloats','MPEG-4');
        set(aviobj,'Quality',100,'FrameRate',4);
        open(aviobj);
    end
    figure
    for i=1:size(eddy.mask,3)
        ff = ceil(i*floats.fac);
        if ff > size(floats.x,1), break; end
        % determine grid indices & make float locations a matrix
        % probably useful later
        xf = (floats.x(ff,:) .* initmask);
        yf = (floats.y(ff,:) .* initmask);
        zf = (floats.z(ff,:) .* initmask);
        
        % remove points just outside grid
        outsideGrid = fillnan((xf<max(xvec(:))),0);
        xf = xf .* outsideGrid; yf = yf .* outsideGrid; zf = zf .* outsideGrid;
        
        % remove floats downstream of eddy
        downstreamEddy = fillnan(xf < eddy.ee(i),0);
        xf = xf .* downstreamEddy; yf = yf .* downstreamEddy; zf = zf .*downstreamEddy;
        
        % keep only 'shallow' floats
        shallow = fillnan(abs(zf) < 10,0);
        xf = xf .* shallow; yf = yf .* shallow; zf = zf .* shallow;
        
        % remove floats inside eddy contour
        netmask = eddy.mask(:,:,i);
        % find float indices
        xg = vecfind(xvec,cut_nan(xf) );
        yg = vecfind(yvec,cut_nan(yf) );
        if ~isempty(xg)
            % nan if float is currently in eddy
            vec_eddy_mask = fillnan(~netmask(sub2ind(sz,xg,yg)),0);  
        else
            vec_eddy_mask = nan(size(xf));
        end
        
        if sum(~isnan(xf)) > 1
            xf = cut_nan(xf) .* vec_eddy_mask; yf = cut_nan(yf) .* vec_eddy_mask;
        end
        % make float matrix
        %fltloc = make_matrix_mask(sz,xvec,yvec,xf,yf);
        
        % after all filtering is done, get area over which I want to
        % integrate
        if ~isempty(cut_nan(xf))
            fltloc = make_matrix_mask(sz,xvec,yvec,cut_nan(xf),cut_nan(yf));
        end
                
        if ~isempty(run.dye)
            if i == 1
                clf;
                hf = pcolor(run.rgrid.xr/1000,run.rgrid.yr/1000,run.dye(:,:,i)/1000);         
                caxis(clim);colormap(flipud(cbrewer('div', 'RdYlGn', 32)));
                shading flat;
                colorbar; hold on; 
                [~,hc] = contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,i),1,'w');
                set(hc,'LineWidth',2);
                hp = plot(xf/1000,yf/1000,'k.','MarkerSize',14);        
                if run.bathy.axis == 'y'
                   hsb =  liney(bathy.xsb/1000,'shelfbreak','w');
                else
                   hsb =  linex(bathy.xsb/1000,'shelfbreak','w');
                end
                ht = title(['Passive tracer + Floats | t = ' num2str(run.rgrid.ocean_time(1)/86400) ' days']);
                %axis image;
                xlabel('X (km)'); ylabel('Y (km)');
                beautify([16 16 18]);
                if run.makeVideo
                    %shading(gca,'interp');
                    disp('maximize!');
                    pause; 
    %               mm_addFrame(run.mm_instance,gcf);
                    F = getframe(gcf);
                    writeVideo(aviobj,F);
                end
            else
                set(hf,'CData',run.dye(:,:,i)/1000);
                set(hc,'ZData',eddy.mask(:,:,i));
                if isempty(hp) % happens if initial xf,yf are NaN
                    hp = plot(xf/1000,yf/1000,'k.','MarkerSize',12);  
                end
                set(hp,'XData',xf/1000);
                set(hp,'YData',yf/1000);
                set(ht,'String',['Passive tracer + Floats | t = ' num2str(run.rgrid.ocean_time(i)/86400) ' days']);
                if run.makeVideo
                    F = getframe(gcf);
                    writeVideo(aviobj,F);
                end
            end
        end
        pause(0.03);
    end
    if run.makeVideo
       % mm_render(run.mm_instance);
       close(aviobj);
    end
end

function [mask] = make_matrix_mask(sz,xvec,yvec,xin,yin)
        xg = vecfind( xvec, xin );
        yg = vecfind( yvec, yin );
        mask = zeros(sz);
        mask(sub2ind(sz,xg,yg)) = 1;
end

function make_vector_mask()

end
