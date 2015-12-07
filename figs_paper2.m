% describe fluxes
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-34')
    ew = runs('../topoeddy/runew-34/');
end

fluxvec = ew.recalculateFlux(2*ew.bathy.hsb,3,3);
ew.calc_avgflux(fluxvec, 1);