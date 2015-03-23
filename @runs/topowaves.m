% topographic waves dispersion relation

function [cgout] = topowaves(runs)

% plot modal structures
%    z=[0:0.05:1];
%    j =1;
%    plot(cos(j*pi*(z-1));
%    hold all
%    plot(cos((j+1/2)*pi*(z-1)));

    plots = 0;

    cgx = nan(1,6);
    f0 = runs.params.phys.f0;
    beta = runs.params.phys.beta;

    % scales
    U = runs.eddy.V(1);
    Ld = runs.rrdeep;
    L = runs.eddy.vor.dia(runs.eddy.tscaleind)*2;
    D = runs.eddy.Lgauss(1);
    N = sqrt(runs.params.phys.N2);

    % bottom slope
    alpha = -1 * runs.bathy.sl_slope;

    % wavelengths
    lambda_x = runs.eddy.vor.lmaj(1);
    lambda_y = runs.eddy.vor.lmin(1);

    k0 = 2*pi ./ (lambda_x./L);
    l0 = 2*pi ./ (lambda_y./L);
    K = sqrt(k0^2 + l0^2); % NON DIMENSIONAL
    S0 = N * D ./ f0 ./ L;

    beta_t = f0./D .* alpha;

    disp(beta_t./beta);

    % QG SOLUTIONS - Pedlosky GFD pg. 411
    %%%%%%%%% bottom trapped mode
    mu = [0.5:0.01:20];
    LHS = tanh(mu);
    RHS = -beta_t./beta .* (mu - S0.*K^2./mu);

    f = @(x) tanh(x) + beta_t./beta .* (x - S0*K^2./x);
    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
    mu0 = fzero(f, 12, options)

    if plots
        figure; subplot(211); hold all
        plot(mu, LHS, mu, RHS);
        linex(mu0);
        legend('LHS', 'RHS');
        xlabel('\mu');
        title('\Psi = A cosh m(z-1)');
        beautify([16 18 20]);
    end

    m = [0.5:0.01:20];
    LHS = tan(m);
    RHS = -beta_t./beta .* (mu + S0.*K^2./m);

    f = @(x) tan(x) + beta_t./beta .* (x + S0*K^2./x);
    options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
    m0 = fzero(f, 12, options)

    if plots
        subplot(212); hold all
        plot(m, LHS, m, RHS);
        linex(m0);
        legend('LHS', 'RHS');
        xlabel('m');
        title('\Psi = A cos m(z-1)');
        beautify([16 18 20]);
    end

    % 1./β in Pedlosky's notation
    Rh = U./beta./L^2;

    % group velocity - μ
    cgx(1) = - 1./Rh .* 1./(K^2 - mu0^2/S0) .* (1 + 2*k0^2./(K^2 - mu0^2/S0)) ...
             .* U

    % topo wave frequency - μ - (non-dimensionalized by f?)
    omega = -(1./Rh) .* k0 ./ (K^2 - mu0^2/S0^2) .* f0;

    % assume wave packet gets to edge of slope after refraction.
    % Frequency is conserved.
    % assume long wavelength, then
    %      ω = - (β L_d^2) k_deep
    % Use this to determine deep water wavelength (all dimensional)
    k_deep = -1 * omega ./ (beta .* Ld^2);
    lambda_deep = 2*pi ./ k_deep;
    fprintf('\n Frequency = %.2e s^(-1)', omega);
    fprintf('\n Period = %f days', 2*pi/omega/86400);
    fprintf('\n Deep water wavelength = %.2f km \n\n', lambda_deep/1000);

    cgout = cgx(1);

    %%%%%%%%% surface intensified modes
    % m = [0.5:0.01:20];
    % LHS = tan(m);
    % RHS = -beta_t./beta .* (m + S0.*K^2./m);

    % f = @(x) tan(x) + beta_t./beta .* (x + S0*K^2./x);
    % options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
    % for ii=1:5
    %     m0(ii) = fzero(f, ((ii-1)+1/2)*pi, options);
    % end

    % cgx(2:end) = - 1./Rh .* 1./(K^2 + m0.^2/S0) .* (1 + 2*k0^2./(K^2 + m0.^2/S0)) ...
    %     .* U

    % figure; hold all
    % plot(m, LHS, '.', m, RHS);
    % linex(m0);
    % legend('LHS', 'RHS');
    % beautify([16 18 20]);

    %%%%%%%%%%%%%%%%%% group velocity for bottom trapped waves

    % % linear solution - Chapman, Rizzoli
    % % need dimensional paramters here
    % k0 = k0/L; l0 = l0/L;
    % Sa0 = alpha*N/f0;

    % syms k l Sa w

    % % ω = σ/f (non-dimensional)
    % ss = 1+Sa^2;
    % kk = 1+k^2/l^2;

    % for ii=1:2
    %     w(k,l,Sa) = (ss/2 + (-1)^(ii) * 1/2*sqrt(ss^2 - 4*Sa^2/kk))^(1/2);

    %     % m ≡ m(k,l,Sa) (non-dimensional vertical co-ordinate z')
    %     m = l/w .* (Sa^2 ./ (1-w^2 + Sa^2))^(1/2);

    %     % group velocity = ∂ω/∂l
    %     cgy = diff(w,l);

    %     w0 = double(w(k0,l0,Sa0));

    %     % non-dimensionalization for z
    %     R = N./f0 * 1./sqrt(1-w0);
    %     m0 = double(m(k0,l0,Sa0)*R);
    %     cgy0 = double(cgy(k0,l0,Sa0) .* f0 );

    %     % output
    %     if ii == 2
    %         disp('====== Positive root');
    %     else
    %         disp('====== Negative root');
    %         cgout = cgy0;
    %     end

    %     disp(['w0 = ', num2str(w0)]);
    %     disp(['vertical wavelength = ', num2str(2*pi./m0) ' m']);
    %     disp(['cg_y = ', num2str(cgy0) ' m/s']);
    % end
end