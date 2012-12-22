load chelton_extract.mat

cyc_mask = fillnan(double(data.cyc == 1),0);
acyc_mask = fillnan(double(data.cyc == -1),0);

[y,m,d,h,mi,s] = datevec(data.date);

%%
date_mask = fillnan(double(y == 2002),0);
r_mask = fillnan(double(data.L >= 70 *1000),0);
n_mask = fillnan(double(data.n >= 3),0);

total_mask = acyc_mask.*date_mask.*n_mask; %.*r_mask;

% should make a shelf_mask

plot_chelton(data,total_mask);