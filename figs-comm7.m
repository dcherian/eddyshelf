if ~exist('ew4341', 'var')
    ew4341 = runs('../topoeddy/runew-4341/');
    ew4341.read_csdsurf;
    ew4341.read_eddsurf;
    ew4341.read_rhosurf;
end

opt.addcsdye = 1;
opt.csfluxplot = 1;
opt.csfluxIsobath = 1;
opt.csdcontours = [97.6 120]*1e3;
opt.limy = [0 350];

ew4341.makeVideo = 1;

ew4341.animate_field('rho', [], 1, Inf, opt);

%%
if ~exist('ew3341', 'var')
    ew3341 = runs('../topoeddy/runew-3341-2-big/');
    ew3341.read_csdsurf;
    ew3341.read_eddsurf;
    ew3341.read_rhosurf;
end

opt.addcsdye = 1;
opt.csfluxplot = 1;
opt.csfluxIsobath = 1;
opt.csdcontours = [97.6 120]*1e3;
opt.limy = [0 350];

ew3341.makeVideo = 1;

ew3341.animate_field('rho', [], 340, 1, opt);
