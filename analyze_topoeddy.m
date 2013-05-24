% analyze topoeddy runs
% makes subplots of eddy tracks

function [] = analyze_topoeddy()
%% read data
dir1 = 'runs/topoeddy/';
dir2 = {'runteb-02';'runteb-03';'runteb-04';'runteb-05';'runteb-06';'runteb-07'};
redo = 0;
for ii=1:length(dir2)
    dir3 = [dir1 dir2{ii} '/'];
    % first look for eddy tracks
    fname = [dir3 '/eddytrack.mat'];
    if ~exist(fname,'file') || redo == 1
        disp(['Calculating track for ' dir3]);
        track{ii}.eddy = track_eddy(dir3);
    else
        track{ii} = load(fname,'eddy');
    end
    % get topo parameters
    params = read_params_from_ini(dir3);
    track{ii}.bathy = params.bathy;
    
    % get non-dim parameters
    fname = [dir3 '/nondim.mat'];
    if ~exist(fname,'file') || redo == 1
       disp(['Calculating non-dim params for ' dir3]);
       track{ii}.nondim = calc_nondim(dir3);
    else
       track{ii}.nondim = load(fname,'nondim');
       track{ii}.nondim = track{ii}.nondim.nondim;
    end
end

%% make plot
nn = length(track);
if mod(nn,2) == 0
    a = factor(nn);
else
    a = factor(nn+1);
end

figure
for ii=1:nn
    % get data
    eddy = track{ii}.eddy;
    nondim = track{ii}.nondim;
    bathy  = track{ii}.bathy;    
    vscale = eddy.Lz2;
    
    % make plot
    subplot(a(1),a(2),ii);
    [c,h] = contour(eddy.xr/1000,eddy.yr/1000,eddy.h(2:end-1,2:end-1),10,'k');
    hold on; caxis([min(vscale(:)) max(vscale(:))]); colorbar;
    scatter(eddy.mx/1000,eddy.my/1000,10,vscale,'filled'); 
    plot(eddy.ee/1000,eddy.my/1000,'k');
    plot(eddy.we/1000,eddy.my/1000,'k');
    plot(eddy.mx(1)/1000,eddy.my(1)/1000,'ko','MarkerSize',7);    
    axis image; xlim([0 220]);
    title(['runteb-0' num2str(ii+1) '| ' num2str(length(eddy.t)) ' days ' ...
            '| sl = ' num2str(bathy.sl_slope) ' | L = ' num2str(bathy.L_slope/1000)]);
        
    if strcmpi(char(dir2{ii}),'runteb-02')
        print_nondim(120,180+100,nondim);
    else
        print_nondim(120,180,nondim);
    end
    xlabel('X (km)');
    ylabel('Y (km)');
end
% correct for bigger domain in runteb-02
if strcmpi(char(dir2{1}),'runteb-02')
    limy = ylim;
    subplot(a(1),a(2),1)
    ylim(limy+100);
end
colormap(pmkmp);

function [] = print_nondim(ix,iy,nondim)
    a = fieldnames(nondim);
    for kk = 1:length(a)
        if ~strcmpi(a{kk},'comment')
            str = sprintf('%s = %.2e',char(a{kk}),nondim.(char(a{kk})));
            if strfind(str,'beta'), str = ['\' str]; end
            ht = text(ix,iy+20*(kk-1),str);
            set(ht,'FontSize',10);
        end
    end