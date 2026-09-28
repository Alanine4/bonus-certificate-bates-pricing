%% main_calibration.m
% Two-step calibration of the Bates model (Heston + Merton jumps) to the
% implied-volatility smile of DAL options, pricing date 16-May-2025,
% expiry 19-Sep-2025.
%   Step 1: Heston-only fit (jump parameters fixed at zero)
%   Step 2: full 8-parameter Bates fit, started from the Step 1 solution
% Output: data/calibrated_params.mat (x_opt)

clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

%% 1. Market & maturity
S0 = 50.46;  r = 0.04321;  q = 0.0132;
T  = (datenum('19-Sep-2025') - datenum('16-May-2025'))/365;

%% 2. Market quotes (calls + puts)
K_call = [20;23;25;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;52.5;55;57.5;60;62.5;65;67.5;70;72.5;75;77.5;80];
P_call = [23.15;28.32;23.30;20.70;20.66;20.15;19.92;18.42;15.67;14.80;11.80;12.65;13.95;13.10;12.35;11.80;11.15;10.64;9.70;8.25;7.55;7.49;6.45;7.27;5.70;4.49;3.55;2.67;1.91;1.55;1.14;0.84;0.62;0.40;0.47;0.26;0.26];
K_put  = [20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;52.5;55;57.5;60;62.5;65;67.5;70;72.5;75;77.5];
P_put  = [0.11;0.55;0.12;0.28;0.24;0.23;0.64;0.48;0.39;0.30;0.30;0.35;0.67;0.53;0.55;0.73;0.81;0.93;1.08;1.30;1.50;1.60;1.87;2.18;2.43;2.93;3.00;3.40;3.75;3.92;4.95;6.05;7.35;9.35;10.60;12.78;17.25;23.30;31.30;27.18;10.43;23.90];

%% Market implied volatilities
strike_vec = [K_call; K_put];
iv_market  = zeros(size(strike_vec));
for i=1:length(K_call)
    iv_market(i) = blsimpv(S0,K_call(i),r,T,P_call(i),'Yield',q,'Class','call');
end
for i=1:length(K_put)
    iv_market(length(K_call)+i) = blsimpv(S0,K_put(i),r,T,P_put(i),'Yield',q,'Class','put');
end
% Keep strikes <= 60 with a finite, positive implied volatility
valid_mask = (strike_vec <= 60) & isfinite(iv_market) & (iv_market > 0);
strike_vec  = strike_vec(valid_mask);
iv_market   = iv_market(valid_mask);

%% 3. Step 1: Heston-only calibration (lambdaJ = muJ = sigmaJ = 0)
% x = [kappa, theta, sigma_v, rho, v0]
x0_h = [0.9, 0.12, 0.15, -0.4, 0.03];
lb_h = [0.1, 0.01, 0.01, -0.9, 0.01];
ub_h = [5.0, 1.0, 1.0, -0.1, 0.5];
opts = optimoptions('lsqnonlin','Display','iter','TolFun',1e-8,'UseParallel',true);

fun_h = @(xh) calibration_err_iv([xh,0,0,0], S0, r, q, T, strike_vec, iv_market);
[x_h, res_h] = lsqnonlin(fun_h, x0_h, lb_h, ub_h, opts);
fprintf('\n=> Heston-only parameters: [%s]\n', sprintf('%.4f ', x_h));

%% 4. Step 2: full Bates calibration
% Start from the Heston solution and add initial jump parameters [lambdaJ, muJ, sigmaJ]
x0_full = [ x_h, 0.1, -0.05, 0.1 ];
lb_full = [ lb_h, 0.00,   -0.2,   0.01 ];
ub_full = [ ub_h, 1.0,     0.2,   0.4  ];

fun_full = @(x) calibration_err_iv(x, S0, r, q, T, strike_vec, iv_market);
[x_opt, resnorm] = lsqnonlin(fun_full, x0_full, lb_full, ub_full, opts);

%% 5. Output final parameters
fprintf('\n===== Calibration Outcome =====\n');
fprintf('params = [ %s];\n', sprintf('%.6f, ', x_opt(1:end-1)) );
fprintf('          %.6f\n', x_opt(end));
fprintf('RSS = %.4e\n', resnorm);
save(fullfile(root, 'data', 'calibrated_params.mat'), 'x_opt');
