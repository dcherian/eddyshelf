function [out] = interpolate(var, z_in, z_out)
    
    sz = size(var);
    if length(sz) == 3
        sz(4) = 1;
    end
    
    var = reshape(var,[sz(1)*sz(2) sz(3) sz(4)]);
    z_in = reshape(z_in, [sz(1)*sz(2) sz(3)]);
    
    len_new = length(z_out);
    
    out = nan([size(var,1) len_new]);
    
    tic;
    parfor ii=1:sz(1)*sz(2)
        out(ii,:) = interp1(z_in(ii,:), var(ii,:), z_out, 'linear');
    end
    toc;
    
    out = reshape(out, [sz(1) sz(2) len_new 1]);
end


% function [out] = interpolate(xmat,ymat,zmat,v,xi,yi,zi)
% 
%     [sz(1),sz(2),sz(3),sz(4)] = size(v);
%     vnew = reshape(v,[sz(1)*sz(2) sz(3)]);
%     
%     % interpolate to same z-levels first
%     for ii=1:length(xi)
%         
%     end
%     
%     % then interpolate in the horizontal on same vertical grid levels
%     for ii=1:sz(3)
%         F = scatteredInterpolant(reshape(xmat(:,:,ii),[sz(1)*sz(2) 1]) ...
%                     ,reshape(ymat(:,:,ii), [sz(1)*sz(2) 1]),vnew(:,ii));
%         out(:,ii) = F(xi,yi); 
%     end
%     
%     if ~isempty(zmat)
%        
%     end
% 
% end