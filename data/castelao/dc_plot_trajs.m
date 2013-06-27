% plots trajectories from renato's output

%% load data
load topo_useast
load aviso_data_1992_2011_eddies_anticyclones_7days_40km_gap_Leff_mod_MAB

%% do traj plots

% plot coast & topo
[cc,hh] = contour(Xuseast,Yuseast,-Zuseast,-[100 200 500 1000 2000 3000 5000],'k');
clabel(cc,hh);
Z_dar
hold on
plot(coast(:,1),coast(:,2),'k')
ylim([34.5 43]);

