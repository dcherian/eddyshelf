function [] = taylorcolumn(runs)

% let's calculate β
% Assumptions:
% - Velocity decays as Gaussian
% - Initial values are good scales for eddy properties
    Ro = runs.eddy.Ro(1);
    if runs.bathy.axis == 'y'
        hvec = runs.bathy.h(1,:);
        axvec = runs.rgrid.y_rho(:,1);
    else
        hvec = runs.bathy.h(:,1)';
        axvec = runs.rgrid.x_rho(1,:)';
    end

    alpha = hvec./max(hvec(:));

    U = runs.eddy.V(1);
    Ubot = U .* exp(-(hvec./runs.eddy.Lgauss(1)).^2);

    beta = alpha./Ro .* Ubot./U;

    figure;
    plot(axvec/1000, beta);
    %linex((runs.bathy.xsb + runs.rrdeep)/1000);
    linex(runs.eddy.my(end)/1000, 'eddy center', 'r');
end