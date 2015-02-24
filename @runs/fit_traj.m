% tanh fit to eddy.my or eddy.hcen
% save output in runs.traj
function [] = fit_traj(runs)

    debug = 0;

    tvec = runs.ndtime;
    % get unique timesteps
    [ut,uind,~] = unique(tvec);
    tvec = tvec(uind);

    use_my = 0;
    if use_my
        yvec = runs.eddy.my(uind);

        % locate origin in the middle of the
        % straight line section.
        iref = find_approx(yvec, ...
                           (yvec(1)+min(yvec))/2, 1);
    else
        yvec = runs.eddy.hcen(uind)';
        i0 = find(yvec == yvec(1), 1, 'last');

        % locate origin in the middle of the
        % straight line section.
        iref = find_approx(yvec, ...
                           (yvec(i0)+min(yvec))/2, 1);
    end

    yref = yvec(iref); tref = tvec(iref);

    tfit = tvec - tref;
    yfit = yvec - yref;
    [y0,T,y1] = runs.tanh_fit(tfit, yfit, 0);

    % location and water depth along idealized
    % trajectory.
    if use_my
        ytraj = y0*tanh(tfit./T) + y1*(tfit./T) + ...
                yref;
        itraj = vecfind(runs.rgrid.y_rho(:,1), ...
                        ytraj);
        htraj = runs.bathy.h(1,itraj);
        str = 'y = eddy.my';
    else
        % clamp to max water depth
        htraj = min(y0*tanh(tfit./T) + y1*(tfit./T) + ...
                    yref, max(runs.bathy.h(:)));
        itraj = vecfind(runs.bathy.h(1,:), ...
                        htraj);
        ytraj = runs.eddy.my(uind);
        str = 'y = hcen';
    end

    % sometimes T comes out as negative. not
    % sure why
    tcrit = 1.5;
    tind = find_approx(tfit, tcrit*abs(T), 1);

    H = htraj(tind);
    Y = ytraj(tind);

    if debug
        figure;
        subplot(211);
        plot(tvec, run.eddy.my(uind), '*', tvec, ytraj);
        liney(Y); title(runName);
        subplot(212);
        plot(tvec, run.eddy.hcen(uind), '*', tvec, ...
             htraj);
        liney(H);
    end
    % y,t scales (corrected for earlier reference shift)
    yscl = y0 + yref;
    tscl = tref + T;

    runs.traj.use_my = use_my;
    runs.traj.yscl = yscl;
    runs.traj.tscl = tscl;
    runs.traj.tcrit = tcrit;
    runs.traj.tind = tind;
    runs.traj.y0 = y0;
    runs.traj.y1 = y1;
    runs.traj.H = H;
    runs.traj.Y = Y;
    runs.traj.htraj = htraj;
    runs.traj.ytraj = ytraj;
    runs.traj.str = str;
    runs.traj.comment = ['Fit y = y0*tanh(t/T) + y1*(t/T) | ' ...
                        'Use y = water depth (hcen) if use_my =0, ' ...
                        'y = eddy.my otherwise | ' ...
                        '(H,Y) = interpolated location / water-depth ' ...
                        'at t = tcrit | (y0,y1,T) are parameters of ' ...
                        'the fit in a different reference frame | ' ...
                        '(yscl,tscl) are (y0,T) corrected for ' ...
                        'reference shift. | (htraj,ytraj) are actual ' ...
                        'fits. ytraj is meaningless when fitting to ' ...
                        'hcen i.e., with use_my = 0'];
    runs.traj.hash = githash([mfilename('fullpath') '.m']);

    traj = runs.traj;
    save([runs.dir '/traj.mat'], 'traj');
end