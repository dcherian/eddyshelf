% plots trajectories from renato's output

%% load data
load topo_useast
load ('03eddy_trajectories_7days_gap40km_1992_Leff_mod_MAB_v5','ea*');
%% do traj plots

% plot coast & topo
plot_coast();

% plot trajs
hold on
colors = distinguishable_colors(length(eaX));
set(gca,'ColorOrder',colors);
for ii=1:length(eaX)/2
    plot(eaX{ii},eaY{ii},'color',colors(ii,:));
    plot(eaX{ii}(1),eaY{ii}(1),'r.','MarkerSize',12);
end


