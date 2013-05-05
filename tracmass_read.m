function [floats] = tracmass_read(fname)

% from the manual
% The trajectories are stored in <outDataDir>/<outDataFile>_run.asc, which has to be
% specified in <project>_run.in
% The trajectories are stored in colums of
%
%              ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens
%
% where
% ntrac is the trajectory number
% niter is the TRACMASS code iteration (only important for TRACMASS modellers)
% x1 is the zoonal position of the trajectory particle
% y1 is the meridional position of the trajectory particle
% z1 is the vertical position of the trajectory particle
% tt the time of the trajectory particle (in days)
% t0 the initial time of the trajectory particle
% subvol is the the "volume" transport in m3/s of the trajectory
% temp is the temperature of the trajectory particle
% salt is the salinity/specific humidity of the trajectory particle
% dens is the density of the trajectory particle
    
   [ntrac,~,x,y,z,tt,t0,~,temp,salt,~] = ...
                textread(fname,'%d%d%f%f%f%f%f%f%f%f%s');
            
   floats.time = unique(tt);
   floats.t0   = unique(t0); % needs to be fixed
   
   for i = 1:length(unique(ntrac)) % ith drifter
      k=1;
      for j=1:length(ntrac)    
        if ntrac(j) == i
           floats.x(k,i) = x(j);
           floats.y(k,i) = y(j);
           floats.z(k,i) = z(j);
           floats.z(k,i) = temp(j);
           floats.z(k,i) = salt(j);
           floats.t(k,i) = tt(j);
           k=k+1;
        end
      end
   end