%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% xymap
name = 'ew-64361';
if ~exist('run', 'var') || ~strcmpi(run.name, name)
    run = runs(['../topoeddy/run' name '/']);
end
fontSize = [20 22 24];
ms = 12; % marker size
trackcolor = [1 1 1]*0.65;
sbslcolor = trackcolor;
tt = [1 250];
seqcolor = flipud(cbrewer('div','RdYlBu',32));
figure; maximize(); pause(0.5);

ax3 = subplot(121);
run.animate_field('eddye', ax3, tt(1), 1);
limx = xlim;
plot(run.eddy.mx/1000, run.eddy.my/1000, 'Color', trackcolor);
plot(run.eddy.mx(tt(1))/1000, run.eddy.my(tt(1))/1000, 'x', ...
     'MarkerSize', ms, 'Color', trackcolor);
text(0.15*diff(limx), run.bathy.xsl/1000, 'slopebreak', ...
     'VerticalAlignment', 'Bottom', 'FontSize', fontSize(1)-4, ...
     'Color', sbslcolor);
text(0.15*diff(limx), run.bathy.xsb/1000, 'shelfbreak', ...
     'VerticalAlignment', 'Top', 'FontSize', fontSize(1)-4, ...
     'Color', sbslcolor);
title('Dyes and SSH'); caxis([-1 1]); beautify(fontSize);
colormap(ax3,seqcolor);
correct_ticks('y',[],[3 6]);

ax4 = subplot(122);
run.animate_field('eddye', ax4, tt(2), 1);
plot(run.eddy.mx/1000, run.eddy.my/1000, 'Color', trackcolor);
plot(run.eddy.mx(tt(2))/1000, run.eddy.my(tt(2))/1000, 'x', ...
     'MarkerSize', ms, 'Color', trackcolor);
caxis([-1 1]);
colormap(ax4,seqcolor); beautify(fontSize);
title('Dyes and SSH'); ylabel([]);
ax4.YTickLabel = [];
correct_ticks('y',[],[3 6]);

export_fig('-r450','images/xymap-poster.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x-y plots of shelfbreak runs
% add colorbar
% figure out best time indices
tinds = [90 90 100 100]; % figure out what these are
figure; maximize();
sb.filter = [2 4 3 5];
fontSize = [16 16 18];
len = length(sb.filter);
for ii=1:len
    run = sb.array(sb.filter(ii));
    run.read_zeta;

    ix = run.spng.sx1:run.spng.sx2;
    iy = run.spng.sy1:run.spng.sy2;

    xr = run.rgrid.x_rho(iy,ix)'/1000;
    yr = run.rgrid.y_rho(iy,ix)'/1000;

    subplot(floor(len/2),2,ii);
    run.animate_field('eddye', gca, tinds(ii), 1);
    %pcolorcen(xr, yr, run.zeta(ix,iy,tinds(ii)));
    hold on;
    plot(run.eddy.mx/1000, run.eddy.my/1000, 'k');
    plot(run.eddy.mx(tinds(ii))/1000, run.eddy.my(tinds(ii))/1000, ...
         'kx');
    if ii == 1, clim = caxis; end
    hold all
    caxis(clim);

    title(['\lambda = ' num2str(run.bathy.hsb/run.eddy.Lgauss(1), ...
                                2) ...
          ' (H_{sb} = ' num2str(run.bathy.hsb, '%.0f') ' m)']);
    ylabel('Y (km)');
    xlabel('X (km)');
    if ii == 3
        correct_ticks('y',[],[4 6]);
    else
        correct_ticks('y',[],[5]);
    end
    beautify(fontSize);
end

tic; export_fig('-r450','images/paper1/sb-maps.png'); toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image effect?
co = image.sorted_colors;
figure; maximize();
ax(2) = subplot(223); hold on;
ax(3) = subplot(224); hold on;
image.flip_colors = 1;

% track
ax(1) = subplot(2,2,[1 2]);
image.plot_penetration(gca);
title([]);
% add timestamps
text(-1.85, 2.75, '50d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-1.85*[1 1], [2.75 2.3], 'k-');
text(-3.45, 2.75, '100d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-3.45*[1 1], [2.75 1.9], 'k-');
text(-3.45, 2.75, '100d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-3.45*[1 1], [2.75 1.9], 'k-');
text(-5.45, 2.75, '200d', 'HorizontalAlignment', 'Center', ...
     'VerticalAlignment', 'bottom');
plot(-5.45*[1 1], [2.75 1.5], 'k-');

% volume, speed diagnostics
co = image.sorted_colors;
for ii=1:image.len
    run = image.array(ii);
    tvec = run.time/run.eddy.turnover;
    nsmooth = 3; run.eddy.turnover./diff(run.time(1:2));

    subplot(223);
    cvx = run.eddy.mvx;
    %cvx(cvx < -0.06) = NaN;
    plot(tvec, smooth(cvx, 8));
    ylabel({'Eddy center along-isobath', 'velocity (km/day)'});

    subplot(224);
    plot(tvec, run.eddy.vol(:,1)./run.eddy.vol(1,1));
    ylabel('Volume / Volume(t=0)');
    xlabel('Time / Turnover Time');
end
subplot(223);
hl = liney(0); ylim([-0.08 0.05]);
uistack(hl, 'bottom'); axis tight;
ax(2).YTick = sort(unique([ax(2).YTick -2.3 -1.5 -0.8]));
correct_ticks('y', '%.1f', {'-2','-1'});
beautify; pbaspect([1.618 1 1]);

subplot(224);
ylim([0 1]);
%liney([0.75 0.4 0.3 0.25]);
ax(3).YTick = sort(unique([ax(3).YTick 0.3  0.75]));
ax(3).XTick = sort(unique([ax(3).XTick 70]));
correct_ticks('y', '%.2f', '0.8');
beautify; pbaspect([1.618 1 1]);
image.reset_colors(co);

linkaxes(ax(2:3), 'x');
export_fig('images/paper1/image-effect.pdf');

%% image effect velocity - sb
image = runArray({ 'runew-2360-fb', 'runew-2360-20km', ...
                   'runew-2360', 'runew-2360_wider'});
for ii=1:image.len
    sh(ii) = image.array(ii).params.bathy.L_shelf;
    if image.array(ii).params.flags.flat_bottom
        sh(ii) = 0;
    end
    image.name{ii} = [num2str(sh(ii)/1000, '%.0f') ' km shelf'];
end
image.sort(sh);
