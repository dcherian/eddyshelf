% creates LTRANS initial location csv file based on ROMS float output

function [] = ltrans_create(fname,flt)

    fid = fopen(fname,'a+');
    if fid == -1
        error('Couldn''t open file to write.');
    end
    
    floats = roms_read_floats(flt);
    
    fprintf(fid,'%.2f,%.2f,%.2f,%d \n',floats.init');
    fprintf('Added %d files to %s \n',size(floats.x,2),fname);
    
%     for ii=1:size(floats.x,2)
%         fprintf(fid,'
%     end
    
    fclose(fid);