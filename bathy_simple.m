% straight isobaths
% axis = x - vertical isobath
%        y - horizontal isobath

function [S] = bathy_simple(S,B,X,Y,axis)

switch axis
    case 'x'
        L = X;
        ax1 = S.x_rho(:,1);
        ax  = S.x_rho;
        
    case 'y'
        L = Y;
        ax1 = S.y_rho(1,:);
        ax  = S.y_rho;
        
    otherwise
        error('invalid axis label provided to bathy_simple_y.m');
end

S.h = zeros(size(ax));

S.h = (H_shelf + B.sl_shelf * (L - ax)) .* (ax > (L-L_shelf));
hmax = max(S.h(:)); %subplot(311); plot(S.h(1,:));
S.h = S.h + (B.sl_slope * ((L-L_shelf) - ax) + hmax) .* ~(ax > (L-L_shelf));  
hdeep = hmax + B.sl_slope * L_slope; %subplot(312); plot(S.h(1,:));
S.h = S.h .* ~(ax < (L-L_shelf-L_slope)) + hdeep * (ax < (L-L_shelf-L_slope));
% subplot(313);plot(S.h(1,:));



% repeat
    