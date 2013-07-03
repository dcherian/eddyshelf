classdef runs < handle
    properties
        % dir & file names
        dir; out_file; ltrans_file; flt_file;
        % data
        zeta; temp;
        % grid & bathymetry
        rgrid; bathy
        % float data
        roms; ltrans;
        % eddy track data
        eddy;
        % initial params
        params
    end
    methods
        % constructor
        function [runs] = runs(dir)
            if isdir(dir)
                runs.dir = dir;
                runs.out_file = [dir '/ocean_avg.nc'];
            else
                runs.out_file = dir;
                dir = strrep(dir,'\','/');
                inds = strfind(dir,'/');
                dir = dir(1:inds(end));
            end
            runs.flt_file = [dir '/ocean_flt.nc'];
            runs.ltrans_file = [dir '/ltrans.nc'];
            runs.rgrid = roms_get_grid(runs.out_file,runs.out_file,0,1);
            runs.rgrid.xr = runs.rgrid.x_rho';
            runs.rgrid.yr = runs.rgrid.y_rho';
            try
                runs.zeta = roms_read_data(dir,'zeta');
            catch ME
                runs.zeta = double(ncread(runs.out_file,'zeta'));
            end
            
            if exist(runs.flt_file,'file')
                runs.roms = floats('roms',runs.flt_file,runs.rgrid);
            end
            
            if ~exist([dir '/eddytrack.mat'],'file')
                try
                    runs.eddy = track_eddy(dir);
                catch ME
                    disp('Couldn''t run track_eddy.m');
                end
            else
                edd = load([dir '/eddytrack.mat'],'eddy');
                runs.eddy = edd.eddy;
            end
            
            if exist(runs.ltrans_file,'file')
                runs.ltrans = floats('ltrans',runs.ltrans_file,runs.rgrid);
            end
            
            runs.params = read_params_from_ini(runs.dir);
            runs.bathy = runs.params.bathy;
            
            % fill bathy
            [runs.bathy.xsb,runs.bathy.isb,runs.bathy.hsb] = find_shelfbreak(runs.out_file);
            [runs.bathy.xsl,runs.bathy.isl,runs.bathy.hsl] = find_shelfbreak(runs.out_file,'slope');
            runs.bathy.h = runs.rgrid.h;
        end
        
        % create initial seed file for ltrans
        function [] = ltrans_create(runs)
            ltrans_create(runs.rgrid,runs.zeta,runs.eddy);
        end
        
        % create ltrans init file from roms out
        function [] = ltrans_create_from_roms(runs)
            ltrans_create_from_roms('ltrans_init_compare.txt',runs.flt_file,runs.rgrid);
        end
    end
end