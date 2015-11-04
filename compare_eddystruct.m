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
ar = 3.337*rn;
ar2 = ar/sqrt(2);
aviso = (1+ar + 1/6 *ar.^2 - 1/6*ar.^3) .* exp(-ar);
aviso2 = (1+ar2 + 1/6 *ar2.^2 - 1/6*ar2.^3) .* exp(-ar2);

figure
hold on
plot(rn,R,'b')
plot(avg1(rn),diff(R)./diff(rn),'r');
plot(rn,aviso,'k');
plot(rn,aviso2,'k--');
legend('T = R(r_n)','V = d/dr_n R(r_n)','AVISO correlation fn','AVISO | r=r/1.414');
liney(0);
title('Zhang et al. (2013)');
%linex( sqrt(7/2 - sqrt(33)/2));
%linex( sqrt(7/2 + sqrt(33)/2));
%linex(4.4)
xlabel('r_n = r/r_0');
beautify([18 18 20]);
maximize
pause(0.1);
export_fig('images/zhangeddy.png');

%% plot Katsman structure

rn = 0:0.05:10;
eddy.a = 2.1;
rmnorm = nthroot( (eddy.a-2)/(eddy.a-1) , eddy.a);
gamma = -2 * (eddy.a-2)./eddy.a ./ (rmnorm)^2;
exponent = (eddy.a - 1)/eddy.a .* (rn.^(eddy.a) - rmnorm.^(eddy.a)); 
R = (gamma/2 * rn.^2 + 1) .* (rn <= rmnorm) ...
                       + (gamma/2 *rmnorm^2 + 1) .* exp( -1 * exponent ) .* ...
                                                           (rn > rmnorm);
                                                       
gaussian = exp(-1/2 * rn.^2);
V = diff(R)./diff(rn);

[~,ind] = min(V);

figure
hold on
plot(rn,R,'b')
plot(avg1(rn),V,'r');
plot(rn,gaussian, 'k');
plot(avg1(rn), diff(gaussian)./diff(rn), 'k--');
linex([2 2.2 2.8 4.3]);
liney(0);
linex(rn(ind));
title('Katsman et al. (2013)');
legend('T = R(r_n)','V = d/dr_n R(r_n)', 'T (gaussian)', 'V (gaussian)');
xlabel('r_n = r/r_0');
%beautify;
%export_fig('images/katsmaneddy.png');

%% Hassanzadeh et al. (2011) paper

z = 0:-0.05:-4;
p = exp(-z.^2);
dpdz = diff(p)./diff(z)

hold all
plot(p, z)
plot(dpdz, avg1(z));

%% gaussian eddy
syms r R
v(r) = diff(exp(-(r/R)^2),r);
vor(r) = diff(v,r);
a = solve(diff(vor,r) == 0, r);
vor(a)