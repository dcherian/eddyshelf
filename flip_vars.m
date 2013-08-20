% flip variables
function [S,axmat] = flip_vars(flip_flag,S,axmat)

    if flip_flag
         S.Trax = permute(S.Trax,[2 1 3]);
         try
            S.Trax_sig = permute(S.Trax_sig,[2 1 3]);
         catch ME
         end
         S.Tra = permute(S.Tra,[2 1 3]);
         
         S.h = S.h';
         S.zeta = S.zeta';
         if exist('axmat','var')
             axmat = permute(axmat,[2 1 3]);
         end
         
         disp('Flipped vars');
    end