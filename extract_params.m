function [] = extract_params(run)

    name = inputname(1);

    % gets all parameters I have saved
    names = evalin('caller', ['fieldnames(' name '.params);']);
    for ii = 1:length(names)
        evalin('caller', [names{ii} ' = ' name '.params.' names{ii} ';']);
    end

    % Now for the ones I don't save - grid params
    names = {'N', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', ...
             'Tcline'};
    for ii = 1:length(names)
        evalin('caller', ['S.' names{ii} ' = ' name '.rgrid.' names{ii} ';']);
    end
    evalin('caller', ['S.Lm = size(' name '.rgrid.mask_rho, 2) - 2;']);
    evalin('caller', ['S.Mm = size(' name '.rgrid.mask_rho, 1) - ' ...
                        '2;']);

    ncid = netcdf.open(run.out_file);
    [~, NT] = netcdf.inqDim(ncid, netcdf.inqDimID(ncid, 'tracer'));
    netcdf.close(ncid);
    spherical = ncread(run.out_file, 'spherical');

    evalin('caller', ['S.NT = ' num2str(NT) ';']);
    evalin('caller', ['S.spherical = ' num2str(spherical) ';']);
    evalin('caller', 'S.NPT = S.NT - 2;');
    evalin('caller', 'calc_pv = 0;');
end
