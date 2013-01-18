% straight isobaths
%       [S] = bathy_simple(S,B,X,Y,axis,loc)
%               axis = x - vertical isobath
%                      y - horizontal isobath
%               loc  = h - coast where axis = high
%                      l - coast where axis = low

function [S] = bathy_simple(S,B,X,Y,axis,loc)

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

    switch loc
        case 'h' % coast at high end of axis
            mask_shelf = (ax > (L-B.L_shelf));
            mask_deep  = (ax < (L-B.L_shelf-B.L_slope));
            ax_shelf = L-ax;

        case 'l'
            mask_shelf = (ax < (B.L_shelf));
            mask_deep  = (ax > (B.L_shelf+B.L_slope));
            ax_shelf = ax;
    end

    S.h = zeros(size(ax));

    S.h = (B.H_shelf + B.sl_shelf * ax_shelf) .* mask_shelf;
    hmax = max(S.h(:)); 
    S.h = S.h + (B.sl_slope *(ax_shelf-B.L_shelf) + hmax) .* ~mask_shelf; 
    hdeep = hmax + B.sl_slope * B.L_slope;
    S.h = S.h .* ~mask_deep + hdeep * mask_deep;
    