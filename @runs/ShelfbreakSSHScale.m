function [] = ShelfbreakSSHScale(runs)

    [varmean, maskmean,xivec,prof] = runs.avgProfile('zeta', runs.bathy.axis, 'sb', 1);

    % find interface. do fit of SSH till the interface
    ind = find(maskmean == 1, 1, 'first') - 1;
    if strcmpi(runs.name, 'ew-8352')
        prof = prof./prof(ind)
    end

    [y0,X,x0,y1,conf,fitobj] = tanh_fit(xivec(10:ind), prof(10:ind), 1);
    title(runs.name);

    sbssh.X = X;
    sbssh.conf = conf(:,2);
    sbssh.fitobj = fitobj;
    sbssh.xvec = xivec;
    sbssh.ssh = prof;

    sbssh.hash = githash([mfilename('fullpath') '.m']);

    save([runs.dir '/sbssh.mat'], 'sbssh');

    runs.sbssh = sbssh;
end
