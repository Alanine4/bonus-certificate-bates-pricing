function err = calibration_err_iv(params, S0, r, q, T, K_vec, iv_market)
    % calibration_err_iv  Compute error between model and market IV using Bates FFT + implied volatility inversion
    %
    % params    - Bates model parameter vector
    % S0        - initial underlying price
    % r, q      - risk-free rate, dividend yield
    % T         - time to maturity (annualized)
    % K_vec     - strike price vector (combined Calls and Puts)
    % iv_market - corresponding market implied volatility vector

    n = numel(K_vec);
    iv_model = zeros(n,1);

    for j = 1:n
        K = K_vec(j);
        try
            % 1) Calculate price using Bates FFT
            price_j = bates_price_fft(params, S0, r, q, T, K);

            % 2) Theoretical arbitrage bounds (meaningful for Calls; same logic applicable for Puts)
            C_min = max(0, S0*exp(-q*T) - K*exp(-r*T));
            C_max = S0*exp(-q*T);

            % 3) Check price validity
            if isnan(price_j) || isinf(price_j) || price_j <= 0 || price_j < C_min || price_j > C_max
                fprintf('Price invalid at strike %.2f -> %.4f  (bounds [%.4f, %.4f])\n', ...
                        K, price_j, C_min, C_max);
                iv_model(j) = iv_market(j);
            else
                % 4) Invert implied volatility
                iv_try = blsimpv(S0, K, r, T, price_j, 'Yield', q);

                % 5) Secondary check of blsimpv output
                if isnan(iv_try) || isinf(iv_try) || iv_try <= 0
                    fprintf('blsimpv failed at strike %.2f (price=%.4f)\n', K, price_j);
                    iv_model(j) = iv_market(j);
                else
                    iv_model(j) = iv_try;
                end
            end

        catch ME
            % Catch any unexpected errors
            fprintf('EXCEPTION at strike %.2f: %s\n', K, ME.message);
            iv_model(j) = iv_market(j);
        end
    end

    % 6) Final cleanup to ensure no NaN/Inf
    bad = ~isfinite(iv_model);
    if any(bad)
        fprintf('Cleaning up %d bad iv_model entries\n', sum(bad));
        iv_model(bad) = iv_market(bad);
    end

    % 7) Return error vector
    err = iv_model - iv_market;
end
