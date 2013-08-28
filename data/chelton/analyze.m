load E:\Work\eddyshelf\data\chelton\chelton_extract.mat

cyc_mask = fillnan(double(data.cyc == 1),0);
acyc_mask = fillnan(double(data.cyc == -1),0);

[y,m,d,h,mi,s] = datevec(data.date);

%%
date_mask = fillnan(double(y == 2003),0);
r_mask = ones(size(data.L)); fillnan(double(data.L >= 40 *1000),0);
n_mask = fillnan(double(data.n >= 3),0);

% shelf mask
lon_min = -75; lon_max = -60;
lat_min =  32; lat_max =  42;
slope = (lat_max - lat_min)/(lon_max-lon_min);
cc = (lon_max * lat_min - lon_min * lat_max)/(lon_max-lon_min);

shelf_mask = (data.lat - slope*data.lon - cc) > 0;

% total mask
total_mask = acyc_mask .* date_mask .* n_mask.* shelf_mask; %.*r_mask;

% make plot
plot_chelton(data,total_mask);
%maximize;
%pause;export_fig -zbuffer 'chelton_test.png';