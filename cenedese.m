cenedese_data;

Q = gc .* hc.^2./(2*f);
Rv = sqrt(gv.*hv)./f;
Rc = sqrt(gc.*hc)./f;
vv = sqrt(gv.*hv); %epsilon .* Q./Rc./hc;

T = Q .* ToQ;

for ii=1:20
    phys.f0 = f(ii);
    eddy.V0 = vv(ii);
    eddy.L0 = Rv(ii);
    eddy.Lz0 = hv(ii);

    flux.IntegrationDepth = hc(ii);
    flux.IsobathLocation = 0;
    flux.IsobathDepth = hc(ii);

    [Flux(ii), FluxError(ii)] = PredictFlux(phys, eddy, flux);
end

%figure;
clf;
hax(1) = subplot(121);
plot(Flux, T, '*');
ylabel('Lab measurement (m^3/s)');
xlabel('Prediction (m^3/s)');
title('all');
beautify;

hax(2) = subplot(122);
plot(Flux(epsilon > 1), T(epsilon > 1), '*');
ylabel('Lab measurement (m^3/s)');
xlabel('Prediction (m^3/s)');
title('\epsilon > 1');
linkaxes(hax, 'xy');
beautify;

export_fig images/cenedese2013.png