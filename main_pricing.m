%% main_pricing.m
% Prices the Bonus Certificate with the calibrated Bates parameters and
% computes its delta for the hedge. Also plots the calibration fit.
% Requires data/calibrated_params.mat (from main_calibration.m).

clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
fig_dir = fullfile(root, 'figures');

%% 1. Calibrated parameters + market/product setup
load(fullfile(root, 'data', 'calibrated_params.mat'), 'x_opt');
params = x_opt;
S0 = 50.46;
r  = 0.04355;
q  = 0.0132;
T  = (datenum('19-Sep-2025') - datenum('15-May-2025'))/365;

%% 2. Market IV surface (calls + puts), strikes with IV = 0 removed
Kc = [20,23,25,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80];
iv_c = [0.8770,0.8555,0.7559,0.6880,0.6602,0.6689,0.6401,0.6167,0.6074,0.5962,0.5952,0.5837,0.5598,0.5491,0.5396,0.5339,0.5232,0.5105,0.5088,0.5071,0.5002,0.4915,0.4849,0.4804,0.4694,0.4618,0.4485,0.4392,0.4321,0.4269,0.4243,0.4292,0.4199,0.4209,0.4224,0.4531,0.4702]';

Kp = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50];
iv_p = [0.7988,0.7647,0.7324,0.7031,0.6934,0.6504,0.6299,0.6143,0.6270,0.5869,0.5889,0.5596,0.5542,0.5381,0.5313,0.5210,0.5112,0.5122,0.5095,0.4915,0.4805,0.4971,0.4695,0.4614,0.4622,0.4502,0.4475,0.4402,0.4351,0.4316,0.4258]';

%% 3. OTM quotes only: calls with K >= F, puts with K <= F, sorted by strike
Kc = Kc(:); iv_c = iv_c(:);
Kp = Kp(:); iv_p = iv_p(:);
F = S0 * exp((r - q) * T);
mask_c = (Kc >= F);
mask_p = (Kp <= F);
strike_vec = [ Kc(mask_c);   Kp(mask_p) ];
iv_mkt     = [ iv_c(mask_c); iv_p(mask_p) ];
[strike_vec, ix] = sort(strike_vec);
iv_mkt = iv_mkt(ix);

%% 4. Product and simulation parameters
B = 1.2 * S0;       % bonus level
H = 0.8 * S0;       % barrier
M = 2e5;            % Monte Carlo paths
N = 200;            % time steps
eps = 0.01 * S0;    % bump size for the delta

%% 5. Plain vanilla European put (put-call parity)
C_B = bates_price_fft(params, S0, r, q, T, B);
plain_put = C_B - S0*exp(-q*T) + B*exp(-r*T);

%% 6. Down-and-out put price and survival probability
rng(2025);
[barrier_put, p_notouched] = do_put_mc(params, S0, r, q, T, B, H, M, N);

%% 7. Bonus Certificate price & delta
[price_BC, delta_BC] = bc_price_mc_delta(params, S0, r, q, T, B, H, M, N, eps);

%% 8. Results
fprintf('\n===== Valuation Result =====\n');
fprintf('Vanilla Put Price        = %.4f\n', plain_put);
fprintf('Down-and-out Put Price   = %.4f\n', barrier_put);
fprintf('Bonus Certificate Price  = %.4f\n', price_BC);
fprintf('BC Delta (S0)            = %.6f\n', delta_BC);
fprintf('Shares needed to buy     = %.0f\n', delta_BC*1e6);
fprintf('Cash needed              = %.2f USD\n', 1e6*price_BC - delta_BC*1e6*S0);
fprintf('Survival Probability     = %.2f%%\n', 100*p_notouched);
fprintf('============================\n\n');

%% 9. Calibration plot: market vs model implied volatility
n = numel(strike_vec);
iv_model = zeros(n,1);
for j = 1:n
    p_j = bates_price_fft(params, S0, r, q, T, strike_vec(j));
    iv_model(j) = blsimpv(S0, strike_vec(j), r, T, p_j, 'Yield', q);
end

fig1 = figure; hold on;
plot(strike_vec, iv_mkt, 'o',  'MarkerSize',6,'LineWidth',1.2);
plot(strike_vec, iv_model,'o','MarkerSize',6,'LineWidth',1.2,'Color',[0.85,0.33,0.1]);
xlabel('Strike Price K'); ylabel('Implied Volatility');
legend('Market IV','Model IV','Location','Best');
title('Calibration Outcome: Market vs. Model IV');
grid on; hold off;
exportgraphics(fig1, fullfile(fig_dir, 'calibration_iv.png'), 'Resolution', 150);

%% 10. Vanilla option prices: model vs market
F = S0 * exp((r - q)*T);
mask_call = strike_vec >= F;
n = numel(strike_vec);

price_mkt   = zeros(n,1);
price_model = zeros(n,1);
for j = 1:n
    Kj = strike_vec(j);
    ivj = iv_mkt(j);

    % Market price: Black-Scholes with the quoted IV
    [c, p] = blsprice(S0, Kj, r, T, ivj, q);
    price_mkt(j) = mask_call(j)*c + (~mask_call(j))*p;

    % Model price: Bates FFT call, puts via put-call parity
    price_call = bates_price_fft(params, S0, r, q, T, Kj);
    if mask_call(j)
        price_model(j) = price_call;
    else
        price_model(j) = price_call ...
                         - S0*exp(-q*T) ...
                         + Kj*exp(-r*T);
    end
end

fig2 = figure; hold on;
plot(strike_vec, price_mkt,   'bo-','LineWidth',1.2, 'MarkerSize',6);
plot(strike_vec, price_model, 'r*-','LineWidth',1.2, 'MarkerSize',6);
xlabel('Strike K');   ylabel('Option Price');
legend('Market Price','Model Price','Location','Best');
title('Vanilla Option Prices: Model vs Market');
grid on; hold off;
exportgraphics(fig2, fullfile(fig_dir, 'vanilla_prices.png'), 'Resolution', 150);
