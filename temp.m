dir = 'runs/runte-11/';
varname = 'temp';

for i=1:11
    fname = [dir sprintf('ocean_his_00%02d',i) '.nc'];
    if i == 1
        v = ncread(fname,varname);
        rgrid = roms_get_grid(fname,fname,0);
    else
        temp = ncread(fname,varname);
        v(:,:,:,end+1:end+size(temp,4)) = temp;
    end
end
%%
%for ii=1:size(v,4)
%    temp = v(:,:,:,ii);
    metric= domain_integrate(v,rgrid.x_rho',rgrid.y_rho',permute(rgrid.z_r,[3 2 1]));
%end

plot(metric);