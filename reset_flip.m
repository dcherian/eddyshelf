function [S,axmat] = reset_flip(S,axmat)

s = size(S.x_rho);

% flip h, zeta, temp, Tx

if size(S.h) ~= s, S.h = S.h'; end
if size(S.zeta) ~= s, S.zeta = S.zeta'; end
if size(S.temp,1) ~= s(1) && size(S.temp,2) ~= s(2)
    S.temp = permute(S.temp,[2 1 3]); 
end
if [size(S.Tx,1) size(S.Tx,2)] ~=s, S.Tx = permute(Tx,[2 1 3]); end

if exist('axmat','var') && size(axmat,1) ~=s(1) && size(axmat,2) ~= s(2)
    axmat = permute(axmat,[2 1 3]);
end

disp('Reset flip');
