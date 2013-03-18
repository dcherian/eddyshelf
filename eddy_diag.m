% Calculates eddy diagnostics as in Chelton et al. (2011)

function [] = eddy_diag()

    load zeta0.mat % in m
    zeta = zeta - mean(zeta(:));
    [L,M,N] = size(zeta);

    thresh_loop = [-100:0]; % in m
    low_n = 10;        % minimum number of pixels in eddy
    high_n = 1000;     % maximum number of pixels in eddy
    amp_thresh = 0.01; % Amplitude threshold (in m)

    for ii=1:length(thresh_loop)
        threshold = thresh_loop(ii);

        % Criterion 1 - first get everything above threshold - is this needed?
        mask = (fillnan(double(zeta > threshold),0));
        
        % Criterion 2 - either too big or too small
        n = sum(isnan(mask(:)));
        if n < low_n || n > high_n, continue; end
        
        % Criterion 3
        if ~local_maximum(z1), continue; end

        % Criterion 4 - amplitude is at least > amp_thresh
        
    
    end
    
    
function [] = local_maximum(zeta)
    
    return 1
