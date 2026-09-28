function price = bates_price_fft(x, S0, r, q, T, K)
% bates_price_fft   Carr–Madan + FFT pricing of European options under the Bates model
% x = [kappa,theta,sigma_v,rho,v0,lambdaJ,muJ,sigmaJ]
% S0 = underlying asset current price, r,q = risk-free rate, dividend yield
% T = time to maturity (annualized scalar), K = strike price (scalar or vector)

    alpha  = 1.5;
    N      = 2^12;
    eta    = 0.25;
    lambda = 2*pi/(N*eta);
    b      = lambda*N/2;

    % v grid & Simpson weights
    v = (0:N-1)' * eta;
    w = ones(N,1);
    w(2:2:end-1) = 4;
    w(3:2:end-2) = 2;
    w = w/3;

    % Carr–Madan integrand
    u   = v - 1i*(alpha+1);
    cf  = charfun_bates(u, x, S0, r, q, T);
    psi = exp(-r*T) .* cf ./ (alpha^2 + alpha - v.^2 + 1i*(2*alpha+1).*v);

    fft_in  = exp(1i*v*b) .* psi .* w * eta;
    fft_out = fft(fft_in);

    % log-strike grid
    k = -b + (0:N-1)' * lambda;

    % Call prices on the grid
    C = exp(-alpha*k)/pi .* real(fft_out);

    price = zeros(size(K));
    for j = 1:numel(K)
        [~, idx] = min(abs(k - log(K(j))));
        price(j) = C(idx);
    end
end
