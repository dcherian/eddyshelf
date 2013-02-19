% create horizontal gradient of tracer about shelfbreak
% ax_as = axis along slope
% ax_cs = axis cross slope

function [Tx] = hor_grad_tracer(axmat_as,ax_as,ax_cs,zrmat,i_cs,i_as,front,B)

    Tx  = nan(size(axmat_as)); %nan([size(axmat_as,1) size(axmat_as,2)]);
    Tx0 = front.Tx0;
    
    scale = (front.LTleft + front.LTright );
    sign  = 1;
        
    n_points = 4; % maybe account for smoothing on shelfbreak 
    
    % code below works if slope is at low end of ax
    % need to flip if high end
    if B.loc == 'h'
        B.h = flipdim(B.h,i_cs);
        %sign = -1; NOT SURE IF THIS IS NEEDED
    end
    
    % first locate shelfbreak as a cross-shelf axis index 
    % at every along-shelf line for reference
    sbreak = zeros([size(B.h,i_as) 1]);
    zsb    = zeros(size(sbreak));
    
    for yy = 1:size(B.h,i_as)
        if i_cs == 1 % vertical isobaths
            % location of shelfbreak (index)
            index = find_approx(B.h(:,yy),B.H_shelf + B.L_shelf * B.sl_shelf,1) ...
                            + floor(n_points/2);
            % location of shelfbreak (distance)
            sbreak(yy,:) = ax_cs(index);
            % water depth at shelfbreak (distance)
            zsb(yy) = squeeze(zrmat(index,1,1));
        else % horizontal isobaths
            index = find_approx(B.h(yy,:),B.H_shelf + B.L_shelf * B.sl_shelf,1) ...
                            + floor(n_points/2);
            sbreak(yy,:) = ax_cs(index);
            % water depth at shelfbreak (distance)
            zsb(yy) = squeeze(zrmat(1,index,1));
        end
    end
    
    if i_cs == 2
        Tx = permute(Tx,[2 1 3]);
        zrmat = permute(zrmat,[2 1 3]);
    end
    
    % move front away from boundary if flat bottom
    if (B.h / max(B.h(:)) == ones(size(B.h)))
        npoints = 60;
        sbreak = ax_cs(npoints * ones(size(sbreak)));
    end
    
    % build gradient
    for jj = length(ax_as):-1:1
        for ii = 1:length(ax_cs)
            for kk=size(axmat_as,3):-1:1
                % introduce frontal slope as shift of "shelfbreak"
                % with depth (shift = 0 at shelfbreak itself)
                shift  = (zrmat(ii,jj,kk) - zsb(jj)) ./ front.slope;
                % for front.slope = 0
                if isinf(shift) || isnan(shift), shift = 0; end
                % change sbreak so that cos-curve center moves also
                center = sbreak(jj) + sign*shift;
                Lleft  = center - front.LTleft;
                Lright = center + front.LTright;
                if ax_cs(ii) <= Lleft || ax_cs(ii) >= Lright
                    Tx(ii,jj,kk) = 0;
                else if ax_cs(ii) > Lleft && ax_cs(ii) < Lright
                        %Tx(ii,jj,kk) = Tx0/2 + Tx0/2*cos(pi*(ax_cs(ii)-center)/scale);
                        Tx(ii,jj,kk) = Tx0*sech((ax_cs(ii)-center)/scale*5).^2;
                    else
                        Tx(ii,jj,kk) = nan;
                    end
                end
                % for debugging
                %sback(ii,jj,kk) = shift;
            end
        end
    end

    if exist('sback','var')
        figure
        mask = squeeze(sback(:,40,:)) > 0;
        contourf(repmat(ax_cs,[1 size(zrmat,3)]),squeeze(zrmat(:,40,:)),mask.*squeeze(sback(:,40,:))/1000,40)
        title('shift in km');colorbar;
    end
    if size(Tx,1) ~= size(axmat_as,1), Tx = permute(Tx,[2 1 3]); end
    
    % flip back if flipped earlier
    if B.loc == 'h'
        Tx = flipdim(Tx,i_cs);
    end