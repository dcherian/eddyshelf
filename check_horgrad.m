Tgrad = diff_cgrid(tgrid,S.temp,1);
Tx0 = avg1(S.Tx,1); %diff(S.temp,1,1)./diff(xrmat,1,1);

Tgrad = squeeze(Tgrad(:,1,:));
Tx0   = squeeze(Tx0(:,1,:));

dT = abs(Tgrad-Tx0);

h(1) = subplot(132);
contourf(avg1(squeeze(xrmat(:,1,:)),1)/1000,avg1(squeeze(zrmat(:,1,:)),1),Tgrad);
title('diff\_cgrid'); colorbar; cx = caxis;

h(2) = subplot(131);
contourf(avg1(squeeze(xrmat(:,1,:)),1)/1000,avg1(squeeze(zrmat(:,1,:)),1),Tx0);
title('diff/imposed'); colorbar; caxis(cx);

h(3) = subplot(133);
contourf(avg1(squeeze(xrmat(:,1,:)),1)/1000,avg1(squeeze(zrmat(:,1,:)),1),dT./max(dT(:))*100);
title('% error'); colorbar;
linkaxes(h,'xy');

%% check_Hz

% check Hz calculation with vertical derivatives
Tgradz = avg1(diff_cgrid(tgrid,S.temp,3),3);
Tz0 = diff(S.temp,1,3)./diff(zrmat,1,3);

Tgradz = squeeze(Tgradz(:,1,:));
Tz0 = squeeze(Tz0(:,1,:));
dTz = abs(Tgradz - Tz0);

figure
h(1) = subplot(131);
contourf(avg1(squeeze(xrmat(:,1,:)),2)/1000,avg1(squeeze(zrmat(:,1,:)),2),Tgradz,20);
title('diff\_cgrid'); colorbar; cx = caxis;

h(2) = subplot(132);
contourf(avg1(squeeze(xrmat(:,1,:)),2)/1000,avg1(squeeze(zrmat(:,1,:)),2),Tz0,20);
title('diff'); colorbar; caxis(cx);

h(3) = subplot(133);
contourf(avg1(squeeze(xrmat(:,1,:)),2)/1000,avg1(squeeze(zrmat(:,1,:)),2),dTz./max(Tz0(:))*100);
title('% error'); colorbar;
linkaxes(h,'xy');