% create horizontal gradient of tracer about shelfbreak
% ax_as = axis along slope
% ax_cs = axis cross slope

function [Tx] = hor_grad_tracer(axmat_as,ax_as,ax_cs,i_cs,i_as,front,B)

    Tx = nan([size(axmat_as,1) size(axmat_as,2)]); nan(size(axmat_as));
    Tx0 = front.Tx0;
    
    scale = (front.LTleft + front.LTright );
    sign  = 1;
    
    % cross-shelf 'axis'
    % ax_cs = [1:size(axmat_as,ax2)]';
    
    n_points = 4; % maybe account for smoothing on shelfbreak 
    
    % code below works if slope is at low end of ax
    % need to flip if high end
    if B.loc == 'h'
        B.h = flipdim(B.h,i_cs);
        sign = -1;
    end
    
    % first locate shelfbreak as a cross-shelf axis index 
    % at every along-shelf line for reference
    sbreak = zeros([size(B.h,i_as) 1]);
    for yy = 1:size(B.h,i_as)
        if i_cs == 1 % vertical isobaths
            sbreak(yy,:) = ax_cs(find_approx(B.h(:,yy),B.H_shelf + B.L_shelf * B.sl_shelf,1) ...
                            + floor(n_points/2));
        else % horizontal isobaths
            sbreak(yy,:) = ax_cs(find_approx(B.h(yy,:),B.H_shelf + B.L_shelf * B.sl_shelf,1) ...
                            + floor(n_points/2));
        end
    end
        
    if i_cs == 2
        Tx = permute(Tx,[2 1 3]);
    end
    
    % build gradient
    for jj = 1:length(ax_as)
        for ii = 1:length(ax_cs)
            %for kk=1:size(sbreak,2)
                % introduce frontal slope as varying reference 'shelfbreak' in depth
                % also convert sbreak to distance
%                sbreak(jj,kk) = 
                Lleft  = sbreak(jj) - front.LTleft;
                Lright = sbreak(jj) + front.LTright;
                if ax_cs(ii) <= Lleft || ax_cs(ii) >= Lright
                    Tx(ii,jj) = 0;
                else if ax_cs(ii) > Lleft && ax_cs(ii) < Lright
                        Tx(ii,jj) = Tx0/2 + Tx0/2*cos(pi*(ax_cs(ii)-sbreak(jj))/scale);
                    else
                        Tx(ii,jj) = nan;
                    end
                end
            %end
        end
    end

    if size(Tx,1) ~= size(axmat_as,1), Tx = permute(Tx,[2 1 3]); end
    
    % flip back if flipped earlier
    if B.loc == 'h'
        Tx = flipdim(Tx,i_cs);
    end