% compare eddy structures for Katsman et al. (2003) & Zhang et al. (2013)

load compare_eddy.mat

figure;
plot(xy_katsman(:,iy));
hold on
plot(xy_zhang(:,iy),'r');
legend('Katsman et al. (2003)','Zhang et al. (2013)');
beautify;
%export_fig('images/compare-eddystruct.png');

%% plot zhang structure

rn = 0:0.05:10;
R = (1-rn.^2/2) .* exp(-rn.^2/2);
dR = -rn .* exp(-rn.^2/2) .* (2 - rn.^2/2);
ar = 3.337*rn*sqrt(2);
aviso = (1+ar + 1/6 *ar.^2 - 1/6*ar.^3) .* exp(-ar);

figure
hold on
plot(rn,R,'b')
plot(avg1(rn),diff(R)./diff(rn),'r');
plot(rn,aviso,'k');
legend('T = R(r_n)','V = d/dr_n R(r_n)','AVISO correlation fn');
liney(0);
title('Zhang et al. (2013)');
linex( sqrt(7/2 - sqrt(33)/2));
linex( sqrt(7/2 + sqrt(33)/2));
linex(4.4)
xlabel('r_n = r/r_0');
beautify;
export_fig('images/zhangeddy.png');

%% plot Katsman structure

rn = 0:0.05:10;
eddy.a = 3;
rmnorm = nthroot( (eddy.a-2)/(eddy.a-1) , eddy.a);
gamma = -2 * (eddy.a-2)./eddy.a ./ (rmnorm)^2;
exponent = (eddy.a - 1)/eddy.a .* (rn.^(eddy.a) - rmnorm.^(eddy.a)); 
R = (gamma/2 * rn.^2 + 1) .* (rn <= rmnorm) ...
                       + (gamma/2 *rmnorm^2 + 1) .* exp( -1 * exponent ) .* ...
                                                           (rn > rmnorm);
figure
hold on
plot(rn,R,'b')
plot(avg1(rn),diff(R)./diff(rn),'r');
linex(2.4);
title('Katsman et al. (2013)');
legend('T = R(r_n)','V = d/dr_n R(r_n)');
xlabel('r_n = r/r_0');
beautify;
export_fig('images/katsmaneddy.png');