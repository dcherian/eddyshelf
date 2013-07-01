% finds shelfbreak
%       [xsb,isb,hsb] = find_shelfbreak(fname)

function [xsb,isb,hsb] = find_shelfbreak(fname)

    h = ncread(fname,'h');
    if h(2,1)-h(1,1) == 0
        yr = ncread(fname,'y_rho');
        dx = yr(1,2)-yr(1,1);
        hvec = h(1,:)';
    else
        xr = ncread(fname,'x_rho'); 
        dx = xr(2,1)-xr(1,1);
        hvec = h(:,1);
    end
   
    dh2dx2 = diff(hvec,2,1)./dx^2;
    [~,isb] = max(dh2dx2);
    isb = isb-1;
    xsb = xr(isb,1);
    hsb = hvec(isb);