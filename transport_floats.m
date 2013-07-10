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
    

    for i=1:size(eddy.mask,3)
        fltloc = zeros(sz);
        ff = ceil(i*floats.fac);
        if ff > size(floats.x,1), break; end
        % determine grid indices & make float locations a matrix
        % probably useful later
        xf = (floats.x(ff,:) .* initmask);
        yf = (floats.y(ff,:) .* initmask);
        
        % remove points just outside grid
        outsideGrid = fillnan((xf<max(xvec(:))),0);
        xf = xf .* outsideGrid; yf = yf .* outsideGrid;
        
        % remove floats downstream of eddy
        downstreamEddy = fillnan(xf < eddy.ee(i),0);
        xf = xf .* downstreamEddy; yf = yf .* downstreamEddy;
        
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
        
        if ~isempty(run.dye)
            if i == 1
                clf;
                hf = pcolor(run.rgrid.xr/1000,run.rgrid.yr/1000,run.dye(:,:,i)/1000);         
                caxis(clim);
                shading flat;
                colorbar; hold on
                [~,hc] = contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,i),1,'k');
                set(hc,'LineWidth',2);
                hp = plot(xf/1000,yf/1000,'k.','MarkerSize',12);        
                if run.bathy.axis == 'y'
                    liney(bathy.xsb/1000,'shelfbreak','w');
                else
                    linex(bathy.xsb/1000,'shelfbreak','w');
                end
                ht = title(num2str(i));
                %axis image;
                beautify;
            else
                set(hf,'CData',run.dye(:,:,i)/1000);
                set(hc,'ZData',eddy.mask(:,:,i));
                if isempty(hp) % happens if initial xf,yf are NaN
                    hp = plot(xf/1000,yf/1000,'k.','MarkerSize',12);  
                end
                set(hp,'XData',xf/1000);
                set(hp,'YData',yf/1000);
                set(ht,'String',num2str(i));
            end
        end
        pause(0.04);
    end
end

function [mask] = make_matrix_mask(sz,xvec,yvec,xin,yin)
        xg = vecfind( xvec, xin );
        yg = vecfind( yvec, yin );
        mask(sub2ind(sz,xg,yg)) = 1;
end

function make_vector_mask()

end
