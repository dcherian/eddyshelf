function [out] = interpolate(var, z_in, z_out)
    
    sz = size(var);
    if length(sz) == 3
        sz(4) = 1;
    end

    if sz(4) == 1
        var = reshape(var,[sz(1)*sz(2) sz(3)]);
        z_in = reshape(z_in, [sz(1)*sz(2) sz(3)]);
    else
        var = permute(var, [1 2 4 3]);
        var = reshape(var, [sz(1)*sz(2)*sz(4) sz(3)]);
        if size(z_in, 4) == 1
            z_in = repmat(z_in, [1 1 1 sz(4)]);
            z_in = reshape(z_in, [sz(1)*sz(2)*sz(4) sz(3)]);
        end
    end

    %var = single(var);
    %z_in = single(z_in);
    %z_out = single(z_out);
    
    len_new = length(z_out);
    
    out = nan([len_new size(var,1)])';
    
    tic;
    parfor ii=1:sz(1)*sz(2)*sz(4)
        out(ii,:) = interp1q(z_in(ii,:)', var(ii,:)', z_out);
    end
    toc;

    if sz(4) == 1
        out = reshape(out, [sz(1) sz(2) len_new 1]);
    else
        out = reshape(out, [sz(1) sz(2) sz(4) len_new]);
        out = permute(out, [1 2 4 3]);
    end
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