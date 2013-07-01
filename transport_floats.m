dir = ['runs/topoeddy/runteb-04-hires-6/'];
track_file = [dir '/ltrans.nc'];
out_file = [dir 'ocean_avg.nc'];
flt_file = [dir 'ocean_flt.nc'];
fname = out_file;

h = ncread(out_file,'h');
rgrid = roms_get_grid(out_file,out_file,1,1);
zrmat = permute(rgrid.z_r,[3 2 1]);
time = ncread(fname,'ocean_time')/86400;
xr = rgrid.x_rho';
yr = rgrid.y_rho';
zeta = double(ncread(out_file,'zeta'));
time = double(ncread(out_file,'ocean_time'))/86400;

[xsb,isb,hsb] = find_shelfbreak(out_file);

roms = floats('roms',flt_file,rgrid);
ltrans = floats('ltrans',track_file,rgrid);
eddy = load([dir '/eddytrack.mat'],'eddy');
eddy = eddy.eddy;

%% masks - doesn't work properly, 

% first bin floats?
xvec = xr(:,1);
yvec = yr(1,:)';
S = size(xr);
floats = ltrans;

initmask = zeros(size(xr));
xi = vecfind(xvec, cut_nan(floats.init(:,1)) );
yi = vecfind(yvec, cut_nan(floats.init(:,2)) );
initmask(sub2ind(S,xi,yi)) = 1;


shelfmask = initmask .* (xr >= xsb); 

for i=1:size(eddy.mask,3)
    clf;
    fltloc = zeros(size(xr));
    ff = i*floats.fac+1;
    % determine grid indices & make float locations a matrix
    xg = vecfind( xvec, cut_nan(floats.x(ff,:)) );
    yg = vecfind( yvec, cut_nan(floats.y(ff,:)) );
    fltloc(sub2ind(S,xg,yg)) = 1;
    
    % remove floats inside eddy contour and those that started offshore of
    % shelbreak
    netmask = ~eddy.mask(:,:,i) .* shelfmask(2:end-1,2:end-1);
    
    % truncate because of eddy.mask
    fltloc = fltloc(2:end-1,2:end-1);    
    
    pcolorcen(eddy.xr/1000,eddy.yr/1000,fltloc ); hold on
    contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,i),1);
    %plot(floats.x(i,:)/1000,floats.y(i,:)/1000,'k.');
    title(num2str(i));
    pause(0.01);
end
