%% read runs
rootdir = 'runs/topoeddy/runew-';
dirs = {... %'02-bg-dt2', '03'
    '06','07','08','09','10'};

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

%% plot tracks
% 
% run(1) = run1;
% run(2) = run2;
% leg = {'07'; '10'}

figure;
subplot(211)
hold on
subplot(212)
hold on
for ii=1:length(run)
    subplot(211)
    plot(run(ii).eddy.cx/1000 - run(ii).params.bg.ubt*[1:length(run(ii).eddy.cx)]*86.4 ...
        - run(ii).params.phys.beta * run(ii).eddy.dia/2, ...
        run(ii).eddy.cy/1000-run(ii).bathy.xsb/1000,'*','Color',colors(ii,:));
    
    subplot(212)
    plot(run(ii).params.nondim.eddy.Rh,run(ii).eddy.trev/86400,'*','Color',colors(ii,:));
end
subplot(211);
title('Weighted center (eddy''s frame of reference)');
legend(leg,'Location','NorthWest');
xlabel('X - u_{bt} x time (km)');
ylabel('Y - Y_{shelfbreak} (km)');
subplot(212)
ylabel('Day of reversal');
xlabel('Rh');

