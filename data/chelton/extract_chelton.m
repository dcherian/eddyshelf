%% read data

fname = 'E:\Work\eddyshelf\data\chelton\AVISO_v3.04wks.nc';

% read and convert to SI units
data.track  = double(ncread(fname,'track'));
data.lon    = double(ncread(fname,'lon')) - 360;
data.lat    = double(ncread(fname,'lat'));
data.A      = double(ncread(fname,'A')) ./ 100;
data.U      = double(ncread(fname,'U')) ./ 100;
data.L      = double(ncread(fname,'L')) .* 1000;
data.cyc    = double(ncread(fname,'cyc'));
data.n      = double(ncread(fname,'n'));
data.j1     = double(ncread(fname,'j1'));

[yr,mo,da,hh,mm,ss] = jd2date(data.j1);
data.date = datenum(yr,mo,da,hh,mm,ss);

% extract region

lonmin = -76;
lonmax = -60;
latmin = 32;
latmax = 46;

mask = (data.lon >= lonmin & data.lon <= lonmax & data.lat >= latmin & data.lat <= latmax);
mask = double(mask);

data = apply_mask(data,mask);

comment = {'extract region from Chelton et al. (2011) dataset + convert to SI units + convert to matlab datenum'};

save chelton_extract.mat data lonmin lonmax latmin latmax comment