% returns structure S used by roms_create
function [S] = extract_params(run)

    % if provided with runs object
    if isobject(run)
        fname = run.out_file;

         % gets all parameters I have saved
        params = run.params;

        grd = run.rgrid;
    else
        % if provided with folder (config/)
        dir = run;
        ininame = [dir '/' roms_find_file(dir, 'ini')];
        grdname = [dir '/' roms_find_file(dir, 'grd')];

        grd = roms_get_grid(grdname, ininame, 0);

        params = read_params_from_ini(ininame);

        % read other info from nc file
        fname = ininame;
    end

    % at this points, I expect to have the following structures
    %   - params, grd
    % and file 'fname' to read other stuff

    names = fieldnames(params);
    % read params into workspace
    for ii = 1:length(names)
        assignin('caller', names{ii}, params.(names{ii}));
    end

    % Now for the ones I don't save - grid params
    names = {'N', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', ...
             'Tcline'};
    for ii = 1:length(names)
        S.(names{ii}) = grd.(names{ii});
    end
    S.Lm =  size(grd.mask_rho, 2) - 2;
    S.Mm = size(grd.mask_rho, 1) - 2;

    ncid = netcdf.open(fname);
    [~, S.NT] = netcdf.inqDim(ncid, netcdf.inqDimID(ncid, 'tracer'));
    netcdf.close(ncid);
    S.NPT = S.NT - 2;
    S.spherical = ncread(fname, 'spherical');

    assignin('caller', 'calc_pv', 0);
end
