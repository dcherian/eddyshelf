% Fit tanh to masked avgProfile of SSH at shelfbreak. Save in sbssh.
% 2*sbssh.X should be an along-shelf length scale of shelf water outflow at the shelfbreak
function [] = ShelfbreakSSHScale(runs)

    [varmean, maskmean,xivec,prof] = runs.avgProfile('zeta', runs.bathy.axis, 'sb', 1);

    % find interface. do fit of SSH till the interface
    ind = find(maskmean == 1, 1, 'first') - 1;
    if strcmpi(runs.name, 'ew-8352')
        prof = prof./prof(ind)
    end

    prof(ind:end) = prof(ind-1); % flat-line downstream of shelf-eddy water front
    [y0,X,x0,y1,conf,fitobj] = tanh_fit(xivec(10:end), prof(10:end), 1);
    title(runs.name);

    sbssh.X = X;
    sbssh.conf = conf(:,2);
    sbssh.fitobj = fitobj;
    sbssh.xvec = xivec;
    sbssh.ssh = prof;

    sbssh.hash = githash([mfilename('fullpath') '.m']);

    save([runs.dir '/sbssh.mat'], 'sbssh');

    disp(['Shelfbreak ΔSSH scale = ' num2str(abs(X/1000), '%.1f') ' km']);
    runs.sbssh = sbssh;
end
