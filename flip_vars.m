% flip variables
function [S] = flip_vars(flip_flag,S,axmat)

    if flip_flag
        S.Tx = S.Tx';
         S.h = S.h';
         S.zeta = S.zeta';
         S.temp = permute(S.temp,[2 1 3]);
         if exist('var','axmat')
             axmat = permute(axmat,[2 1 3]);
         end
    end