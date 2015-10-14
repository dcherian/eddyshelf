% plots trajectories from renato's output

%% load data
load topo_useast
load ('03eddy_trajectories_7days_gap40km_1992_Leff_mod_MAB_v5','ea*');
%% do traj plots
figure
% plot coast & topo
plot_coast();

% plot trajs
colors = distinguishable_colors(length(eaX));
set(gca,'ColorOrder',colors);
for ii=1:length(eaX)/2
    plot(eaX{ii},eaY{ii},'b');%'color',colors(ii,:));
    plot(eaX{ii}(1),eaY{ii}(1),'r.','MarkerSize',12);
end
xlabel('lon'); ylabel('lat'); beautify([14 14 16]);

%% build mask
load topo_useast
N = ceil(length(eaX));
mask = ones([N 1]);
cx = mask; cy = mask;
% get centers
for ii=1:N
    cx(ii) = eaX{ii}(1); cy(ii) = eaY{ii}(1);
    
    % remove shelf
    ix = find_approx(Xuseast(:,1),cx(ii));
    iy = find_approx(Yuseast(1,:),cy(ii));
    
    if Zuseast(ix,iy) < 110,
        mask(ii) = 0;
    end
    
    % remove south of 36N
    if cy(ii) < 36
        mask(ii) = 0;
    end
end
%%
figure
% plot coast & topo
plot_coast();

% plot trajs
colors = distinguishable_colors(length(eaX));
set(gca,'ColorOrder',colors);
for ii=1:length(eaX)/2
    if mask(ii)
        plot(eaX{ii}*mask(ii),eaY{ii}*mask(ii),'b');%'color',colors(ii,:));
        plot(eaX{ii}(1)*mask(ii),eaY{ii}(1)*mask(ii),'r.','MarkerSize',12);
    end
end
xlabel('lon'); ylabel('lat'); beautify([14 14 16]);