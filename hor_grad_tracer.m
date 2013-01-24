% create horizontal gradient of tracer about shelfbreak
function [Tx] = hor_grad_tracer(h,ax,ax1,ax2,sbreak,front)

    Tx = nan(size(h));
    Tx0 = front.Tx0;
    
    for jj = 1:size(h,ax2)
        for ii = 1:size(h,ax1)        
            Lleft = ax(sbreak(jj)) - front.LT;
            Lright = ax(sbreak(jj)) + front.LT;
            if ax(ii) <= Lleft || ax(ii) >= Lright
                Tx(ii,jj) = 0;
            else if ax(ii) > Lleft && ax(ii) < Lright
                    Tx(ii,jj) = Tx0/2 + Tx0/2*cos(pi*(ax(ii)-ax(sbreak(jj)))/front.LT);
                else
                    Tx(ii,jj) = nan;
                end
            end
        end
    end

    if size(Tx,1) ~= size(h,1), Tx = Tx'; end