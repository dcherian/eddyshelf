function [sbreak,ind] = find_shelfbreak(fname)

    h = ncread(fname,'h');
    xr = ncread(fname,'x_rho');    
    dx = xr(2,1)-xr(1,1);    
    dhdx = diff(h(:,1),1,1)./dx;
    [~,ind] = max(dhdx);
    sbreak = xr(ind,1);