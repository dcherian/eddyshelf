% Figure out relation b/w horizontal scales and vertical scales.
Z = 2000; % total depth
zvec = 1:10:Z;

f = 1e-4;


N2 = 1e-5 * avg1(exp(-zvec/Z*20));
[Vmode,Hmode,c] = vertmode(N2,zvec,1,1);

% rossby radius from eigenvalue (km)
lc = c/f/1000;
% rossby radius from profile
vscale = zvec(find_approx(Hmode,0,1));
%vscale = zvec(find_approx(N2,max(N2)/2.718,1));
liney(vscale);
lp = mean(sqrt(N2)*vscale/f)/1000;
fprintf(' Theoretical (cn/f) = %.2f km | Inferred = NZ/f = %.2f km | error= %.2f %%\n', ...
    lc,lp,abs(lp-lc)/lc*100);
