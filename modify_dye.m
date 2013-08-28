%% read data and initialize
function [] = modify_dye(file)
    if isempty(file)
        file = 'runs/ocean_rst.nc';
    end

    [xrmat,yrmat,zrmat,~,~,~] = dc_roms_var_grid(file,'rho');
    [xsb,isb,hsb,ax] = find_shelfbreak(file);
    %[xsl,isl,hsl,ax] = find_shelfbreak(file,'slope');

    % set condition

    if ax == 'y'
        dye_01 = yrmat;
        dye_03 = xrmat;
    else
        dye_01 = xrmat;
        dye_03 = yrmat;
    end

    dye_02 = zrmat;
    % dye_01 = zeros(size(xrmat)); dye_02 = dye_01;
    % 
    % buffer = 6;
    % % decrease from 1 to 0
    % buffer_val_dec = repmat(cos( pi/2*(0:buffer-1)/(buffer-1) ), ...
    %                             [size(dye_01,1) 1 size(dye_01,3)]);
    % % increase from 0 to 1
    % buffer_val_inc = repmat(fliplr(cos( pi/2*(0:buffer-1)/(buffer-1) )), ...
    %                             [size(dye_01,1) 1 size(dye_01,3)]);
    %                         
    % dye_01(:,1:isb,:) = 1;
    % dye_01(:,isb+1:isb+buffer,:) = buffer_val_dec;
    % 
    % dye_02(:,isb-buffer:isb-1,:) = buffer_val_inc;
    % dye_02(:,isb:isl,:) = 1;
    % dye_02(:,isl+1:isl+buffer,:) = buffer_val_dec;
    %                         
    % 
    % subplot(211)
    % pcolorcen(dye_01(:,:,end)'); 
    % title('dye_{01}');
    % subplot(212)
    % pcolorcen(dye_02(:,:,end)');
    % title('dye_{02}');

    % write to file
    nt = length(ncread(file,'ocean_time'));
    ncwrite(file,'dye_01',repmat(dye_01,[1 1 1 2 nt]));
    ncwrite(file,'dye_02',repmat(dye_02,[1 1 1 2 nt]));
    try
    ncwrite(file,'dye_03',repmat(dye_03,[1 1 1 2 nt]));
    catch ME
        warning('dye_03 not found');
    end
    disp(['wrote to file ' file]);