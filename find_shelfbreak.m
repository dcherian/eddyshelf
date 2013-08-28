% finds shelfbreak or end of continental slope (type = 'slope')
%       [xsb,isb,hsb,ax] = find_shelfbreak(fname,type)

function [xsb,isb,hsb,ax] = find_shelfbreak(fname,type)

    h = ncread(fname,'h');
    if h(2,1)-h(1,1) < 1e-3
        xr = ncread(fname,'y_rho')';
        hvec = h(1,:)';
        ax = 'y';
    else
        xr = ncread(fname,'x_rho'); 
        hvec = h(:,1);
        ax = 'h';
    end
    dx = xr(2,1)-xr(1,1);
    dh2dx2 = diff(hvec,2,1)./dx^2;
    if exist('type','var') && strcmpi(type,'slope')
        [~,isb] = min(dh2dx2);
    else
        [~,isb] = max(dh2dx2);
    end
    
    isb = isb-1;
    xsb = xr(isb,1);
    hsb = hvec(isb);