% given a directory, calculates all required non-dimensional parameters

function [nondim] = calc_nondim(dir)
    
    params = read_params_from_ini(dir);
    if ~isfield(params,'phys') % run is from before I added physical parameters
       params.phys.f0   = 1.028e-4;
       params.phys.beta = 2e-11;
       params.phys.N2   = 1e-5;
    end
    
    phys = params.phys; bathy = params.bathy; eddy = params.eddy;
    
    N = sqrt(phys.N2);
    H = bathy.sl_shelf * bathy.L_shelf + ...
        bathy.sl_slope * bathy.L_slope;
    
    % estimate eddy parameters - need U
    if ~isfield(eddy,'U')
        fname = [dir '/' roms_find_file(dir,'ini')];
        N = length(ncread(fname,'s_rho'));
        u = ncread(fname,'u',[1 1 N 1],[Inf Inf 1 1]);
        eddy.U = max(u(:));
        clear u
        % save to ini file
        write_params_to_ini(fname,eddy.U,'eddy.U');
    end
    
    % topographic beta (total H)
    nondim.betatH = phys.f0/H * bathy.sl_slope;
    
    % topographic beta (H = fL/N)
    nondim.betafL = N/bathy.L_slope * bathy.sl_slope;
    
    % slope Burger number
    nondim.S = bathy.sl_slope * N / phys.f0;
    
    % Burger number
    nondim.Bu = N*eddy.depth/phys.f0 / (eddy.dia/2);
    
    % Rossby Number
    nondim.Ro = eddy.U/phys.f0 / (eddy.dia/2);
    
    % Rhines number? U/beta/L^2
    nondim.Rh = eddy.U/phys.beta / (eddy.dia/2)^2;
    
    % description
    nondim.comment = ['betatH = f0/H (dh/dx) - for slope, H = total water depth | ' ...
                     'betafL = " with H = fL/N | ' ...
                     'S = slope  burger number = slope_slope * N/f | ' ...
                     'Bu = eddy Burger number = N*eddy.depth/(f*eddy.radius) | ' ...
                     'Ro = Rossby Number = U/f0/eddy.radius | ' ...
                     'Rh = Rhines Number? = U/beta/(eddy.radius)^2'];
                 
    % save to file
    save([dir '/nondim.mat'],'nondim');
    