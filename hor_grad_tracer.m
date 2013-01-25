% create horizontal gradient of tracer about shelfbreak
function [Tx] = hor_grad_tracer(axmat,ax,ax1,ax2,sbreak,front)

    Tx = nan([size(axmat,1) size(axmat,2)]);
    Tx0 = front.Tx0;
    
    scale = (front.LTleft + front.LTright );
    
    for jj = 1:size(axmat,ax2)
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

    if size(Tx,1) ~= size(ax,1), Tx = Tx'; end