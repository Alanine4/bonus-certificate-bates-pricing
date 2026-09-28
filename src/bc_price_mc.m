function price = bc_price_mc(params, S0, r, q, T, B, H, M, N)
% bc_price_mc  Monte Carlo price of a Bonus Certificate under the Bates model
% Euler scheme with M paths and N steps. Payoff: max(S_T, B) if the barrier H
% was never touched, S_T otherwise.

    % Unpack Bates parameters
    kappa   = params(1);
    theta   = params(2);
    sigma_v = params(3);
    rho     = params(4);
    v0      = params(5);
    lambdaJ = params(6);
    muJ     = params(7);
    sigmaJ  = params(8);
    mJ      = exp(muJ + 0.5 * sigmaJ^2) - 1;

    dt = T / N;
    sqrt_dt = sqrt(dt);
    payoffs = zeros(M,1);

    for m = 1:M
        S = S0;
        v = v0;
        touched = false;

        for i = 1:N
            dWv = sqrt_dt * randn;
            dWs = rho * dWv + sqrt(1 - rho^2) * sqrt_dt * randn;
            if rand < lambdaJ * dt
                jump = muJ + sigmaJ * randn;
            else
                jump = 0;
            end
            v = max(v + kappa*(theta - v)*dt + sigma_v*sqrt(max(v,0))*dWv, 0);
            S = S * exp((r - q - lambdaJ * mJ - 0.5 * v) * dt + sqrt(max(v,0)) * dWs + jump);

            if S < H
                touched = true;
            end
        end

        % Bonus Certificate payoff
        if ~touched
            if S > B
                payoffs(m) = S;
            else
                payoffs(m) = B;
            end
        else
            payoffs(m) = S;
        end
    end

    price = exp(-r * T) * mean(payoffs);
end
