%% read runs
rootdir = 'runs/topoeddy/runew-';
dirs = {... %'02-bg-dt2', '03'
    '02','03','04','05','06-vis3','08','09'};

leg = dirs;
% run(1) = runs('runs/topoeddy/runew-02-bg-dt2/');
% run(2) = runs('runs/topoeddy/runew-03/');
% run(3) = runs('runs/topoeddy/runew-06/');
% run(3) = runs('runs/topoeddy/runew-05/');
% run(4) = runs('runs/topoeddy/runew-06/');
% %run(5) = runs('runs/topoeddy/runew-07/');

for ii=1:length(dirs)
    run(ii) = runs([rootdir dirs{ii} '/']);
end
colors = distinguishable_colors(length(run));

%% plot prox

clf;
hold on
for ii=1:length(run)
    plot(run(ii).eddy.prox/1000,'Color',colors(ii,:));
end
legend(leg);
ylabel('Proximity (km)');
xlabel(' Time (days)');

%% plot initial vel. 

tend = 30; nsmooth = 3;

fig;
subplot(121)
hold on
for ii=1:length(run)
    plot(run(ii).eddy.t(1:tend),smooth(run(ii).eddy.cvy(1:tend),nsmooth),'Color',colors(ii,:));
end
legend(cellstr(dirs),'Location','SouthWest')
ylabel('weighted center meridional vel (km/day)');
xlabel('time (days)');


vs = nan([1 length(run)]); v = vs; x = vs; y = vs;
% plot parameterization
subplot(122);
hold on
for ii=1:length(run)
    % deformation radius
    Lr = run(ii).params.eddy.dia/2;
    % rossby wave speed = beta LR^2
    vr = run(ii).params.phys.beta * (Lr)^2;
    % gravity wave speed = NH
    vgw = sqrt(run(ii).params.phys.N2) * max(run(ii).bathy.h(:));
    % eddy amplitude
    amp = run(ii).eddy.amp(1);
    
    H = vgw^2 / run(ii).params.phys.g;
    % v scale based on early et al. (2011)
    vs(ii) = 14.5*86.4*vr * (vr/vgw) * run(ii).params.nondim.eddy.Rh;%
    Rh(ii) = run(ii).params.nondim.eddy.Rh;
    % obs vel
    v(ii) = nanmax(abs(smooth(run(ii).eddy.cvy(1:tend),nsmooth)));
    %nanmean(run(ii).eddy.cvy(1:15))
    
    x(ii) = abs(v(ii));
    y(ii) = abs(vs(ii));
    text(x(ii),y(ii),dirs{ii});
    %set(gca,'XTickLabel', cellstr(dirs));
end
ratio = v./vs;
plot(x,y,'.','MarkerSize',15);
axis square
%plot(x,y/x*x,'r*-')
plot(xlim,xlim,'k-'); 
title(['slope = ' num2str(x/y)]);
xlabel('max. v (from model, km/day)');
ylabel('v_s = 14.5 * v_r^2/v_{gw} *Rh + 0.5 (km/day)');


%% plot tracks
% 
% run(1) = run1;
% run(2) = run2;
% leg = {'07'; '10'}

h1 = figure;
hold on
h2 = figure;
hold on
h3 = figure;
hold on
for ii=1:length(run)
    figure(h1)
    plot(run(ii).eddy.cx/1000 - (run(ii).params.bg.ubt ...
        )...%- run(ii).params.phys.beta * (run(ii).params.eddy.dia/2).^2) ...
        *[1:length(run(ii).eddy.cx)]*86.4, ...
        run(ii).eddy.cy/1000-run(ii).bathy.xsb/1000,'*','Color',colors(ii,:));
    
    figure(h3)
    %ndtime = (run(ii).eddy.cx - run(ii).eddy.cx(1))./ ...
    %                (run(ii).params.bg.ubt -  ...
    %                run(ii).params.phys.beta/2*(run(ii).eddy.dia/2).^2);
    plot(run(ii).eddy.prox/1000,run(ii).eddy.Lz2./run(ii).eddy.Lz2(1),'*','Color',colors(ii,:));
    
    figure(h2)
    plot(run(ii).params.nondim.eddy.Rh,run(ii).eddy.trev/86400,'*','Color',colors(ii,:));
end
loc = 'SouthWest';
figure(h1)
title('Weighted center (eddy''s frame of reference)');
legend(leg,'Location',loc);
xlabel('X - u_{bt} x time (km)');
ylabel('Y - Y_{shelfbreak} (km)');
figure(h2)
ylabel('Day of reversal');
xlabel('Rh');
legend(leg,'Location',loc);
figure(h3)
ylabel('Vertical scale fraction');
xlabel('Time (days)');
legend(leg,'Location',loc);

