% create horizontal gradient of tracer about shelfbreak
function [Tx] = hor_grad_tracer(axmat,ax,ax1,ax2,front,B)

    Tx = nan([size(axmat,1) size(axmat,2)]);
    Tx0 = front.Tx0;
    
    scale = (front.LTleft + front.LTright );
    
    % code below works if slope is at low end of ax
    % need to flip if high end
    if B.loc == 'h'
        B.h = flipdim(B.h,ax1);
    end
    
    n_points = 3;
    
    % first locate shelfbreak
    sbreak = zeros([size(B.h,ax2) 1]);
    for yy = 1:size(B.h,ax2)
        if ax1 == 1 % vertical isobaths
            sbreak(yy) = find_approx(B.h(:,yy),B.H_shelf + B.L_shelf * B.sl_shelf,1) ...
                            + floor(n_points/2);
        else % horizontal isobaths
            sbreak(yy) = find_approx(B.h(yy,:),B.H_shelf + B.L_shelf * B.sl_shelf,1) ...
                            + floor(n_points/2);
        end
    end
    
    if ax1 == 2
        Tx = Tx';
    end
    
    % build gradient
    for jj = 1:length(sbreak)
        for ii = 1:size(axmat,ax1)        
            Lleft = ax(sbreak(jj)) - front.LTleft;
            Lright = ax(sbreak(jj)) + front.LTright;
            if ax(ii) <= Lleft || ax(ii) >= Lright
                Tx(ii,jj) = 0;
            else if ax(ii) > Lleft && ax(ii) < Lright
                    Tx(ii,jj) = Tx0/2 + Tx0/2*cos(pi*(ax(ii)-ax(sbreak(jj)))/scale);
                else
                    Tx(ii,jj) = nan;
                end
            end
        end
    end

    if size(Tx,1) ~= size(axmat,1), Tx = Tx'; end
    
    % flip back if flipped earlier
    if B.loc == 'h'
        Tx = flipdim(Tx,ax1);
    end