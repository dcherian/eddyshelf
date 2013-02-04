% flip variables
function [S,axmat] = flip_vars(flip_flag,S,axmat)

    if flip_flag
         S.Tx = permute(S.Tx,[2 1 3]);
         S.h = S.h';
         S.zeta = S.zeta';
         S.temp = permute(S.temp,[2 1 3]);
         if exist('axmat','var')
             axmat = permute(axmat,[2 1 3]);
         end
         
         disp('Flipped vars');
    end