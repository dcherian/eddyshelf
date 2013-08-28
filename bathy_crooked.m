function [hx,hy,hdeep] = bathy_crooked(x_rho,y_rho,B,X,Y)

    hx = nan(size(y_rho));
    hy = hx;
    
    hx = (B.H_shelf + B.sl_shelf * (Y - y_rho)) .* (y_rho > (Y-B.L_shelf));
    hmax = max(hx(:)); %subplot(311); plot(hx(1,:));
    hx = hx + (B.sl_slope * ((Y-B.L_shelf) - y_rho) + hmax) .* ~(y_rho > (Y-B.L_shelf));  
    hdeep = hmax + B.sl_slope * B.L_slope; %subplot(312); plot(hx(1,:));
    hx = hx .* ~(y_rho < (Y-B.L_shelf-B.L_slope)) + hdeep * (y_rho < (Y-B.L_shelf-B.L_slope));
   % subplot(313);plot(hx(1,:));
   
  	%B.L_slope = B.L_slope - 2500;
    hy = 0 * (x_rho > (X-B.L_entry));
    hy = hy + (-B.sl_slope2 * (x_rho- (X-B.L_entry))./B.L_slope) .* (x_rho < (X-B.L_entry) & x_rho > (X-B.L_entry-B.L_slope));  
    hmax = max(hy(:));
    hy = hy + (hmax) .* (x_rho < (X-B.L_entry-B.L_slope));
    hy = hy./max(hy(:));
    hy(y_rho < Y-B.L_slope-B.L_shelf) = 0;
    
    
%     % old version
%     hx = zeros(size(y_rho));
%     hy = hx;
%     
%     hx = (B.H_shelf + B.sl_shelf * (Y - y_rho)) .* (y_rho > (Y-B.L_shelf));
%     hmax = max(hx(:)); %subplot(311); plot(hx(1,:));
%     hx = hx + (B.sl_slope * ((Y-B.L_shelf) - y_rho) + hmax) .* ~(y_rho > (Y-B.L_shelf));  
%     hdeep = hmax + B.sl_slope * B.L_slope; %subplot(312); plot(hx(1,:));
%     hx = hx .* ~(y_rho < (Y-B.L_shelf-B.L_slope)) + hdeep * (y_rho < (Y-B.L_shelf-B.L_slope));
%    % subplot(313);plot(hx(1,:));
%    
%   	B.L_slope = B.L_slope - 2500;
%     hy = 0 * (x_rho > (X-B.L_entry));
%     hy = hy + (-B.sl_slope * (x_rho- (X-B.L_entry))./B.L_slope) .* (x_rho < (X-B.L_entry) & x_rho > (X-B.L_entry-B.L_slope));  
%     hmax = max(hy(:));
%     hy = hy + (hmax) .* (x_rho < (X-B.L_entry-B.L_slope));
%     hy = hy./max(hy(:));
%     hy(y_rho < Y-B.L_slope-B.L_shelf) = 0;