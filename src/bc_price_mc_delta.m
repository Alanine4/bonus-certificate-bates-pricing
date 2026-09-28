function [price, delta] = bc_price_mc_delta(params, S0, r, q, T, B, H, M, N, eps)
% bc_price_mc_delta  Bonus Certificate price and delta under the Bates model
% Simulates paths from S0+eps and S0-eps with common random numbers and
% takes the central difference for the delta.
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
    
    payoff_up   = zeros(M,1);
    payoff_down = zeros(M,1);

    rng(123);

    for m = 1:M
        S_up   = S0 + eps;
        S_down = S0 - eps;
        v1 = v0;
        v2 = v0;
        touched_up = false;
        touched_down = false;

        for i = 1:N
            dWv = sqrt_dt * randn;
            dWs = rho * dWv + sqrt(1 - rho^2) * sqrt_dt * randn;
            jump_flag = rand < lambdaJ * dt;
            jump_val = jump_flag * (muJ + sigmaJ * randn);

            v1 = max(v1 + kappa * (theta - v1) * dt + sigma_v * sqrt(max(v1,0)) * dWv, 0);
            S_up = S_up * exp((r - q - lambdaJ * mJ - 0.5 * v1) * dt + sqrt(max(v1,0)) * dWs + jump_val);
            if S_up < H, touched_up = true; end

            v2 = max(v2 + kappa * (theta - v2) * dt + sigma_v * sqrt(max(v2,0)) * dWv, 0);
            S_down = S_down * exp((r - q - lambdaJ * mJ - 0.5 * v2) * dt + sqrt(max(v2,0)) * dWs + jump_val);
            if S_down < H, touched_down = true; end
        end

        if ~touched_up
            if S_up > B
                payoff_up(m) = S_up;
            else
                payoff_up(m) = B;
            end
        else
            payoff_up(m) = S_up;
        end

        if ~touched_down
            if S_down > B
                payoff_down(m) = S_down;
            else
                payoff_down(m) = B;
            end
        else
            payoff_down(m) = S_down;
        end
    end

    price = exp(-r * T) * mean( (payoff_up + payoff_down) / 2 );
    delta = exp(-r * T) * mean( (payoff_up - payoff_down) / (2 * eps) );
end