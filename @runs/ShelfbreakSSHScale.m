% Fit tanh to masked avgProfile of SSH at shelfbreak. Save in sbssh.
% 2*sbssh.X should be an along-shelf length scale of shelf water outflow at the shelfbreak
function [] = ShelfbreakSSHScale(runs)

    [sshmean, maskmean, xivec, prof] = runs.avgProfile('zeta', runs.bathy.axis, 'sb', 1);

    % find interface. do fit of SSH till the interface
    ind = find(maskmean == 1, 1, 'first') - 1;
    if strcmpi(runs.name, 'ew-8352')
        prof = prof./prof(ind)
    end

    % For mean SSH, let's flatline at the maximum ssh point.
    sshmeanback = sshmean;
    [~,imax] = max(sshmean);
    sshmean(imax+1:end) = sshmean(imax);
    [y0,X,x0,y1,conf,fitobj] = ...
        tanh_fit(xivec(10:end), sshmean(10:end), 1);
    plot(xivec(10:end), sshmeanback(10:end), 'k--');
    title([runs.name ' | unmasked']);

    sbssh.nomask.X = X;
    sbssh.nomask.conf = conf(:,2);
    sbssh.nomask.fitobj = fitobj;
    sbssh.xvec = xivec;
    sbssh.nomask.ssh = prof;

    % for masked averaged SSH, flatline downstream of dye front
    prof(ind:end) = prof(ind-1);
    [y0,X,x0,y1,conf,fitobj] = ...
        tanh_fit(xivec(10:end), prof(10:end), 1);
    plot(xivec(10:end), sshmeanback(10:end), 'k--');
    title([runs.name ' | masked']);

    sbssh.X = X;
    sbssh.conf = conf(:,2);
    sbssh.fitobj = fitobj;
    sbssh.ssh = prof;

    sbssh.hash = githash([mfilename('fullpath') '.m']);

    save([runs.dir '/sbssh.mat'], 'sbssh');

    disp(['Shelfbreak ΔSSH scale (unmasked) = ' ...
          num2str(abs(sbssh.nomask.X/1000), '%.1f') ' km']);
    disp(['Shelfbreak ΔSSH scale (shelf-water) = ' ...
          num2str(abs(X/1000), '%.1f') ' km']);
    runs.sbssh = sbssh;
end
