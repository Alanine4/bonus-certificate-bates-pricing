function [price, p_notouched] = do_put_mc(params, S0, r, q, T, K, H, M, N)
% do_put_mc  Monte Carlo price of a down-and-out put (strike K, barrier H)
% and the probability that the barrier is never touched.
    % Constant volatility approximation using the Heston initial variance v0
    v0    = params(5);
    sigma = sqrt(v0);

    dt = T / N;
    payoffs      = zeros(M,1);
    survived     = false(M,1);

    for m = 1:M
        S       = S0;
        touched = false;
        for i = 1:N
            Z = randn;
            S = S * exp((r - q - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z);
            if S < H
                touched = true;
                break;
            end
        end
        if ~touched
            survived(m)    = true;
            payoffs(m)     = max(K - S, 0);
        end
    end

    price       = exp(-r * T) * mean(payoffs);
    p_notouched = mean(survived);
end
