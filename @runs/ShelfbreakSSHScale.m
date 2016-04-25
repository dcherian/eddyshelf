function [] = ShelfbreakSSHScale(runs)

% Use profile that is the time-average of (SSH .* mask)
    [~,~,xivec,prof] = runs.avgProfile('zeta', runs.bathy.axis, 'sb', 1);

    [y0,X,x0,y1,conf,fitobj] = tanh_fit(xivec, prof);

    sbssh.X = X;
    sbssh.conf = conf(:,2);
    sbssh.fitobj = fitobj;
    sbssh.xvec = xivec;
    sbssh.ssh = prof;

    sbssh.hash = githash([mfilename('fullpath') '.m']);

    save([runs.dir '/sbssh.mat'], 'sbssh');

    runs.sbssh = sbssh;
end
