% Integrates to get 2 layer system
%   [dlayer1,dlayer2] = integ2layer(run,varname,tindices)
%           run -> class
%           varname -> variable name : can be 'vor'
%           tindice -> time indices

function [dlayer1,dlayer2] = integ2layer(run,varname,tindices,volume)

    if ~exist('volume','var'), volume = {}; end
    if ~exist('tindices','var'), tindices = []; end
           
    data = dc_roms_read_data(run.dir,varname,tindices,volume);
    stride = [1 1 1 1];
    h = run.bathy.h(1:stride(1):end,1:stride(2):end);
    zeta = run.zeta(1:stride(1):end,1:stride(2):end,:);

    for i=1:size(run.bathy.h,1)-size(data,1)
        h = avg1(h,1);
        zeta = avg1(zeta,1);
    end
    for i=1:size(run.bathy.h,2)-size(data,2)
        h = avg1(h,2);
        zeta = avg1(zeta,2);
    end

    Z = max(h(:));
    disp('Integrating...');
    tic;
    
    dlayer1 = nan([size(data,1) size(data,2) size(data,4)]);
    dlayer2 = dlayer1;
    
    for ii=1:size(data,4)
        [~,dlayer2(:,:,ii)] = roms_depthIntegrate(data(:,:,:,ii),run.rgrid.Cs_r, ...
                run.rgrid.Cs_w, h,zeta(:,:,ii),[-Z -Z/2]);
        [~,dlayer1(:,:,ii)] = roms_depthIntegrate(data(:,:,:,ii),run.rgrid.Cs_r, ...
            run.rgrid.Cs_w, h,zeta(:,:,ii),[-Z/2 0]);
    end
    toc;
    disp('Done');    
    end