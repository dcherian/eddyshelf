% creates LTRANS initial location csv file based on ROMS float output
%   [] = ltrans_create(fname,flt,rgrid)
function [] = ltrans_create(fname,flt,rgrid)

    fid = fopen(fname,'a+');
    if fid == -1
        error('Couldn''t open file to write.');
    end
    
    floats = read_floats('roms',flt,rgrid);
    
    fprintf(fid,'%.2f,%.2f,%.2f,%d \n',floats.init');
    fprintf('Added %d files to %s \n',size(floats.x,2),fname);

    fclose(fid);