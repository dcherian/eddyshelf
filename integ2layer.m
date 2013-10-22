function [dlayer1,dlayer2] = integ2layer(run,varname,tindices)
           
           read_start = [ 1 1 1 tindices(1)];
           read_count = [Inf Inf Inf tindices(2)-tindices(1)+1];
           stride = [1 1 1 1];
           
           data = roms_read_data(run.dir,varname,read_start,read_count,stride);
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