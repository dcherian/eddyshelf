function [cgout] = topowaves(runs)

    cgx = nan(1,6);
    f0 = runs.params.phys.f0;
    beta = runs.params.phys.beta;

    % scales
    U = runs.eddy.V(1);
    L = runs.rrdeep;
    D = 1200; runs.eddy.Lgauss(1);
    N = sqrt(runs.params.phys.N2);

    % bottom slope
    alpha = -1 * runs.bathy.sl_slope;

    % wavelengths
    lambda_x = runs.eddy.vor.lmaj(1);
    lambda_y = runs.eddy.vor.lmin(1);

    k0 = 2*pi ./ (lambda_x./L);
    l0 = 2*pi ./ (lambda_y./L);
    K = sqrt(k0^2 + l0^2);
    S0 = N * D ./ f0 ./ L;

    beta_t = f0./D .* alpha;

    disp(beta_t./beta);
% $$$     %{ QG SOLUTIONS
% $$$     %%%%%%%%% bottom trapped mode
% $$$     mu = [0.5:0.01:20];
% $$$     LHS = tanh(mu);
% $$$     RHS = -beta_t./beta .* (mu - S0.*K^2./mu);
% $$$
% $$$     f = @(x) tanh(x) + beta_t./beta .* (x - S0*K^2./x);
% $$$     options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
% $$$     mu0 = fzero(f, 12, options)
% $$$
% $$$     figure; hold all
% $$$     plot(mu, LHS, mu, RHS);
% $$$     linex(mu0);
% $$$     legend('LHS', 'RHS');
% $$$     beautify([16 18 20]);
% $$$
% $$$     Rh = U./beta./L^2;
% $$$     cgx(1) = - 1./Rh .* 1./(K^2 - mu0^2/S0) .* (1 + 2*k0^2./(K^2 - mu0^2/S0)) ...
% $$$           .* U
% $$$
% $$$     %%%%%%%%% surface intensified modes
% $$$     m = [0.5:0.01:20];
% $$$     LHS = tan(m);
% $$$     RHS = -beta_t./beta .* (m + S0.*K^2./m);
% $$$
% $$$     f = @(x) tan(x) + beta_t./beta .* (x + S0*K^2./x);
% $$$     options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
% $$$     for ii=1:5
% $$$         m0(ii) = fzero(f, ((ii-1)+1/2)*pi, options);
% $$$     end
% $$$
% $$$     cgx(2:end) = - 1./Rh .* 1./(K^2 + m0.^2/S0) .* (1 + 2*k0^2./(K^2 + m0.^2/S0)) ...
% $$$         .* U
% $$$
% $$$     figure; hold all
% $$$     plot(m, LHS, '.', m, RHS);
% $$$     linex(m0);
% $$$     legend('LHS', 'RHS');
% $$$     beautify([16 18 20]);
% $$$ %}
% $$$     %%%%%%%%%%%%%%%%%% group velocity for bottom trapped waves

    % need dimensional paramters here
    k0 = k0/L; l0 = l0/L;
    Sa0 = alpha*N/f0;

    syms k l Sa w

    % ω = σ/f (non-dimensional)
    ss = 1+Sa^2;
    kk = 1+k^2/l^2;

    for ii=1:2
        w(k,l,Sa) = (ss/2 + (-1)^(ii) * 1/2*sqrt(ss^2 - 4*Sa^2/kk))^(1/2);

        % m ≡ m(k,l,Sa) (non-dimensional vertical co-ordinate z')
        m = l/w .* (Sa^2 ./ (1-w^2 + Sa^2))^(1/2);

        % group velocity = ∂ω/∂l
        cgy = diff(w,l);

        w0 = double(w(k0,l0,Sa0));

        % non-dimensionalization for z
        R = N./f0 * 1./sqrt(1-w0);
        m0 = double(m(k0,l0,Sa0)*R);
        cgy0 = double(cgy(k0,l0,Sa0) .* f0 );

        % output
        if ii == 2
            disp('====== Positive root');
        else
            disp('====== Negative root');
            cgout = cgy0;
        end

        disp(['w0 = ', num2str(w0)]);
        disp(['vertical wavelength = ', num2str(2*pi./m0) ' m']);
        disp(['cg_y = ', num2str(cgy0) ' m/s']);
    end
end