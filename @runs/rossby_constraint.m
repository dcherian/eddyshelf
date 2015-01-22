function [] = rossby_contraint(runs)

    if isempty(runs.zeta)
        runs.read_zeta;
    end

    % define co-ordinate system with respect to eddy center
    % (x,y,t)
    x = bsxfun(@minus, runs.rgrid.x_rho', ...
               permute(runs.eddy.mx, [3 1 2]));
    y = bsxfun(@minus, runs.rgrid.y_rho', ...
               permute(runs.eddy.my, [3 1 2]));
    r = sqrt(x.^2 + y.^2);

    %    x = runs.rgrid.x_rho(1,:)';
    %y = runs.rgrid.y_rho(:,1)';

    r = sqrt(bsxfun(@plus, x(2:end-1,2:end-1,:).^2, y(2:end-1,2:end-1,:).^2));
    f = runs.rgrid.f(2:end-1, 2:end-1)';

    dzdx = bsxfun(@rdivide, diff(runs.zeta, 1, 1), diff(x,1,1));
    dzdy = bsxfun(@rdivide, diff(runs.zeta, 1, 2), diff(y,1,2));

    dzdr = bsxfun(@times, avg1(dzdx(:,2:end-1,:),1), ...
                  bsxfun(@rdivide, y(2:end-1,2:end-1,:), r)) + ...
           bsxfun(@times, avg1(dzdy(2:end-1,:,:),2), ...
                  bsxfun(@rdivide, x(2:end-1,2:end-1,:), r));

    diag = 4*9.81*bsxfun(@rdivide, dzdr./r,f.^2);


end