function cf = charfun_bates(u, x, S0, r, q, T)
% charfun_bates  Compute characteristic function under the Bates model
% u - complex frequency vector; x = [kappa,theta,sigma_v,rho,v0,lambdaJ,muJ,sigmaJ]

    kappa   = x(1);
    theta   = x(2);
    sigma_v = x(3);
    rho     = x(4);
    v0      = x(5);
    lambdaJ = x(6);
    muJ     = x(7);
    sigmaJ  = x(8);

    % Jump drift compensation
    mJ = exp(muJ + 0.5*sigmaJ^2) - 1;

    % Complex unit multiplied by u
    iu = 1i .* u;

    % Intermediate variables for Heston part
    A = kappa - rho * sigma_v .* iu;
    d = sqrt((rho * sigma_v .* iu - kappa).^2 + sigma_v^2 .* u .* (u + 1i));
    g = (A - d) ./ (A + d);
    g(abs(g) > 1 - 1e-8) = 1 - 1e-8;

    % Exponent part of characteristic function in Carr–Madan formula
    C = (r - q - lambdaJ*mJ) .* (iu * T) ...
      + (theta*kappa/sigma_v^2) .* ( (A - d)*T - 2 .* log((1 - g .* exp(-d*T)) ./ (1 - g)) );

    % Variance term
    D = ((A - d)./sigma_v^2) .* ((1 - exp(-d*T)) ./ (1 - g .* exp(-d*T))) .* v0;

    % Jump part
    J = lambdaJ * T .* (exp(iu * muJ - 0.5 * sigmaJ^2 .* u.^2) - 1);

    % Final characteristic function
    C(isnan(C)) = 0; D(isnan(D)) = 0; J(isnan(J)) = 0;
    C(~isfinite(C)) = 0;
    D(~isfinite(D)) = 0;
    J(~isfinite(J)) = 0;

    cf = exp(C + D + J + iu * log(S0));

end
