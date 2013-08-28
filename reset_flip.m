function [S,axmat] = reset_flip(S,axmat)

s = size(S.x_rho);

% flip h, zeta, temp, Tx

if size(S.h) ~= s, S.h = S.h'; end
if size(S.zeta) ~= s, S.zeta = S.zeta'; end
if size(S.Tra,1) ~= s(1) && size(S.Tra,2) ~= s(2)
    S.Tra = permute(S.Tra,[2 1 3]); 
end
if size(S.Trax,1) ~=s(1) ||  size(S.Trax,2) ~=s(2), S.Trax = permute(S.Trax,[2 1 3]); end

if exist('axmat','var') && size(axmat,1) ~=s(1) && size(axmat,2) ~= s(2)
    axmat = permute(axmat,[2 1 3]);
end

disp('Reset flip');
