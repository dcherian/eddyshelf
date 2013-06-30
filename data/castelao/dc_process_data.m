%% process anticyclonic tracks
load topo_useast
load ('03eddy_trajectories_7days_gap40km_1992_Leff_mod_MAB_v5','ea*');

%% year by year

nfigs = 2;
years = 1992:2011;

aa = 2; bb = 5;

for kk = 1:nfigs
    figure;
    for ii=1:length(years)/nfigs;
        nn = (kk-1)* length(years)/nfigs + ii;
        subplot(aa,bb,ind2sub([aa bb],ii));
        plot_coast(20,0);
        
        % now plot trajs
        
        % clean up fig
        Z_dar;
        title(num2str(years(nn)));
    end
    spaceplots
end