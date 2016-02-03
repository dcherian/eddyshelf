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
opt.limy = [0 300];

ew4341.makeVideo = 1;

ew4341.animate_field('eddye', [], 179, 1, opt);