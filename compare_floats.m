function [] = compare_floats()

%dir = ['runs/topoeddy/runteb-04-hires-6/'];
dir = ['runs/topoeddy/runbathysouth-02/'];
track_file = [dir '/ltrans-compare.nc'];
out_file = [dir 'ocean_avg.nc'];
flt_file = [dir 'ocean_flt.nc'];
fname = out_file;

h = ncread(out_file,'h');
rgrid = roms_get_grid(out_file,out_file,1,1);
zrmat = permute(rgrid.z_r,[3 2 1]);
time = ncread(fname,'ocean_time')/86400;
xr = rgrid.x_rho';
yr = rgrid.y_rho';

% locate shelfbreak
isb = find_approx(h(:,1),120,1);
xsb = xr(isb,1);

zeta = double(ncread(out_file,'zeta'));
time = double(ncread(out_file,'ocean_time'))/86400;

%% read all float data
clear ltrans roms
ltrans = floats('ltrans',track_file,rgrid);
roms   = floats('roms',flt_file,rgrid);
ltrans.plot_stats; roms.plot_stats
ltrans.plot_displacements; roms.plot_displacements;
% fig;
% plot(ltrans.time/86400,ltrans.N); hold on
% plot(roms.time/86400,roms.N,'r'); 
% xlabel('time (days)'); ylabel('N = number of particles');
% legend('LTRANS','ROMS');


%% find tracmass floats that correspond to ROMS deployment
% dx = xr(2,1)-xr(1,1);
% dy = yr(1,2)-yr(1,1);
% 
% if roms.init(1,4) > 1000
%     roms.init(:,4) = floor(roms.init(:,4)/86400);
% end
% mask = zeros([1 size(trac.x,2)]);
% match = 0;
% 
% for mm = 1:size(roms.x,2)
%     for nn = 1:size(trac.x,2)
%         if hypot(floor(roms.init(mm,1)) - floor(trac.init(nn,1)), ...
%                  floor(roms.init(mm,2)) - floor(trac.init(nn,2))) < 20*sqrt(dx*dy) && ...
%            floor(roms.init(mm,4)) == floor(trac.init(nn,4))
%             match = match + 1;
%             mask(1,nn) = 1;
%         end
%     end
% end
% mask = repmat(mask,[size(trac.x,1) 1]);
% %mask = ones(size(mask));
% %
% i = 1;
% figure
% hold on
% plot3(roms.x/1000,roms.y/1000,roms.z,'r')
% plot3(trac.x/1000 .* mask,trac.y/1000 .*mask ,trac.z .*mask,'k');
% plot3(trac.init(:,1)/1000,trac.init(:,2)/1000,trac.init(:,3),'kx','MarkerSize',12);
% linex(xsb/1000,'Shelfbreak');
% title('red = ROMS, black = TRACMASS');
