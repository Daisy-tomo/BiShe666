%% ========================================================================
% bayesian_inversion_single_condition.m
%
% 【程序说明】单工况贝叶斯反演 —— 双涵道涡扇发动机热力循环参数识别
%
% 本程序的核心任务：
%   - 观测量：单一飞行工况下的"比推力 R_ud"与"比油耗 C_ud"
%   - 未知量：13 个发动机修正参数（效率、总压恢复系数等）
%   - 方法：通过先验分布（均匀分布）+ 单工况高斯似然函数，
%           构造参数的后验分布，再用 Metropolis-Hastings MCMC 采样
%   - 目标：不是给出一个单一的"最优值"，而是得到参数的后验分布，
%           用于分析参数可辨识性、参数间相关性以及不确定性量化
%
% 【单工况设定】
%   增压比 pi_k      = 33
%   涡轮前温度 T_g   = 1800 K
%   飞行马赫数       = 0（地面静止）
%   飞行速度         = 0 m/s
%   环境温度 T_H     = 288 K
%   涵道比 m         = 8（与多工况程序一致）
%
% 使用"虚拟试验数据"验证：
%   - 指定真值参数 theta_true → 前向模型 → 生成带噪声观测
%   - MCMC 反演 → 后验分布应覆盖真值
%
% 【程序结构】
%   主程序 + 局部函数（单文件）
%   局部函数列表：
%     engine_forward, generate_virtual_data,
%     log_prior, log_likelihood, log_posterior_fn,
%     run_mcmc, reflect_to_unit_interval,
%     piecewise_kT, piecewise_RT, delta_cooling,
%     summarize_posterior, plot_results
%
% 【采样配置】
%   轻量模式（默认）：12000 次采样，3000 burn-in
%   高精度模式（注释）：可改为 30000~50000 次采样
%
% 作者：基于图片公式严格编码
% 日期：2026-03
% ========================================================================

clear; clc; close all;
rng(42);  % 固定随机种子，保证可重复

%% ========================================================================
%  第一部分：定义参数名称、上下界、真值、初值
%% ========================================================================

param_names = {
    'eta\_k',       ...  % 1  压气机效率
    'eta\_t',       ...  % 2  涡轮效率
    'eta\_m',       ...  % 3  机械效率
    'eta\_v',       ...  % 4  风扇效率
    'eta\_tv',      ...  % 5  风扇涡轮效率
    'eta\_c1',      ...  % 6  一次喷管效率（内涵）
    'eta\_c2',      ...  % 7  二次喷管效率（外涵）
    'sigma\_cc',    ...  % 8  进气道总压恢复系数
    'sigma\_kan',   ...  % 9  进气通道总压恢复系数
    'sigma\_kask',  ...  % 10 压气机级间总压恢复系数
    'sigma\_ks',    ...  % 11 燃烧室总压恢复系数
    'eta\_T',       ...  % 12 燃烧放热系数
    'lambda\_star'  ...  % 13 风扇涡轮热恢复系数
};

% 先验范围（均匀分布）
lb = [0.84,  0.86,  0.980, 0.85,  0.90,  0.94,  0.92,  0.98,  0.98,  0.98,  0.94,  0.97,  1.02];
ub = [0.92,  0.92,  0.995, 0.87,  0.92,  0.95,  0.94,  1.00,  0.99,  0.99,  0.96,  0.99,  1.04];

n_params = length(lb);

% 真值参数（非边界，偏中间取值，供虚拟数据生成用）
theta_true = [
    0.90,   ...  % eta_k
    0.89,   ...  % eta_t
    0.990,  ...  % eta_m
    0.860,  ...  % eta_v
    0.910,  ...  % eta_tv
    0.945,  ...  % eta_c1
    0.930,  ...  % eta_c2
    0.990,  ...  % sigma_cc
    0.985,  ...  % sigma_kan
    0.985,  ...  % sigma_kask
    0.950,  ...  % sigma_ks
    0.980,  ...  % eta_T
    1.030   ...  % lambda_star
];

% 初始点：先验中点（可替换为真值附近小扰动）
theta0 = (lb + ub) / 2;

fprintf('=== 参数初始化 ===\n');
fprintf('%-15s  %-8s  %-8s  %-8s  %-8s\n', '参数', '下界', '上界', '真值', '初值');
for i = 1:n_params
    fprintf('%-15s  %-8.4f  %-8.4f  %-8.4f  %-8.4f\n', ...
        param_names{i}, lb(i), ub(i), theta_true(i), theta0(i));
end

%% ========================================================================
%  第二部分：定义单工况
%  地面静止起飞工况：pi_k=33, T_g=1800K, M=0, T_H=288K, m=8
%% ========================================================================

cond.T_H       = 288;    % 环境温度 [K]
cond.M_flight  = 0.0;    % 飞行马赫数（地面静止）
cond.m         = 8;      % 涵道比
cond.pi_k      = 33;     % 增压比
cond.T_g       = 1800;   % 涡轮前温度 [K]

fprintf('\n=== 单工况参数 ===\n');
fprintf('  环境温度   T_H      = %.1f K\n',  cond.T_H);
fprintf('  飞行马赫数 M_flight = %.2f\n',    cond.M_flight);
fprintf('  涵道比     m        = %d\n',       cond.m);
fprintf('  增压比     pi_k     = %d\n',       cond.pi_k);
fprintf('  涡前温度   T_g      = %.0f K\n',  cond.T_g);

%% ========================================================================
%  第三部分：用真值参数验证前向模型，生成虚拟试验数据
%% ========================================================================

noise_level_R = 0.01;  % 比推力相对噪声标准差 1%
noise_level_C = 0.01;  % 比油耗相对噪声标准差 1%
rng_seed      = 2024;

data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed);

fprintf('\n=== 虚拟试验数据 ===\n');
fprintf('  R_true  = %.4f  N·s/kg\n',  data.R_true);
fprintf('  R_obs   = %.4f  N·s/kg  (加噪声后)\n', data.R_obs);
fprintf('  C_true  = %.6f kg/(N·h)\n', data.C_true);
fprintf('  C_obs   = %.6f kg/(N·h)  (加噪声后)\n', data.C_obs);
fprintf('  sigma_R = %.4f  N·s/kg\n',  data.sigma_R);
fprintf('  sigma_C = %.6f kg/(N·h)\n', data.sigma_C);

%% ========================================================================
%  第四部分：设置 MCMC 参数
%% ========================================================================

opts.n_samples          = 12000;  % 总采样步数（含 burn-in）
opts.burn_in            = 3000;   % 燃烧期长度
opts.proposal_sd        = 0.015;  % 归一化空间中的初始步长
opts.adapt_start        = 500;    % 自适应调节开始步
opts.adapt_end          = 3000;   % 自适应调节结束步
opts.adapt_interval     = 200;    % 每隔多少步调整一次
opts.target_accept_low  = 0.20;   % 目标接受率下限
opts.target_accept_high = 0.35;   % 目标接受率上限
opts.adapt_factor       = 1.10;   % 调节幅度因子

% ---- 高精度模式（取消注释以启用）----
% opts.n_samples     = 50000;
% opts.burn_in       = 10000;
% opts.proposal_sd   = 0.012;
% opts.adapt_start   = 1000;
% opts.adapt_end     = 10000;
% opts.adapt_interval= 500;

%% ========================================================================
%  第五部分：运行 MCMC
%% ========================================================================

fprintf('\n=== 开始 MCMC 采样 ===\n');
fprintf('总步数: %d，burn-in: %d，后验样本数: %d\n', ...
    opts.n_samples, opts.burn_in, opts.n_samples - opts.burn_in);

tic;
results = run_mcmc(data, cond, lb, ub, theta0, opts);
t_mcmc = toc;
fprintf('MCMC 用时: %.1f 秒\n', t_mcmc);

%% ========================================================================
%  第六部分：输出后验统计
%% ========================================================================

fprintf('\n=== 后验统计结果 ===\n');
fprintf('总体接受率: %.4f\n', results.accept_rate);
fprintf('\n%-15s  %-8s  %-8s  %-8s  %-8s  %-10s  %-10s\n', ...
    '参数', '真值', '后验均值', 'MAP', '95%下', '95%上', '后验标准差');
for i = 1:n_params
    fprintf('%-15s  %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-10.4f  %-10.6f\n', ...
        param_names{i}, theta_true(i), results.theta_mean(i), results.theta_map(i), ...
        results.theta_ci95(i,1), results.theta_ci95(i,2), results.theta_std(i));
end

%% ========================================================================
%  第七部分：后验预测验证（后验均值 + MAP 回代）
%% ========================================================================

fprintf('\n=== 后验预测验证 ===\n');

y_mean = engine_forward(results.theta_mean, cond);
y_map  = engine_forward(results.theta_map,  cond);

R_pred_mean = y_mean(1);
C_pred_mean = y_mean(2);
R_pred_map  = y_map(1);
C_pred_map  = y_map(2);

fprintf('%-12s  %-12s  %-14s  %-14s\n', '量', '观测值', '后验均值预测', 'MAP预测');
fprintf('%-12s  %-12.4f  %-14.4f  %-14.4f\n', 'R_ud [N·s/kg]',  data.R_obs, R_pred_mean, R_pred_map);
fprintf('%-12s  %-12.6f  %-14.6f  %-14.6f\n', 'C_ud [kg/N·h]',  data.C_obs, C_pred_mean, C_pred_map);

%% ========================================================================
%  第八部分：绘图
%% ========================================================================

plot_results(results, data, cond, param_names, theta_true, lb, ub, ...
    R_pred_mean, C_pred_mean, R_pred_map, C_pred_map);

fprintf('\n=== 程序运行完毕 ===\n');

%% ========================================================================
%%  局部函数区域
%% ========================================================================

% ------------------------------------------------------------------------
% 函数：piecewise_kT
% 用途：分段计算燃气多变指数 k_T(T_g)
% ------------------------------------------------------------------------
function kT = piecewise_kT(T_g)
    if T_g > 800 && T_g <= 1400
        kT = 1.33;
    elseif T_g > 1400 && T_g <= 1600
        kT = 1.30;
    elseif T_g > 1600
        kT = 1.25;
    else
        kT = 1.33;
    end
end

% ------------------------------------------------------------------------
% 函数：piecewise_RT
% 用途：分段计算燃气气体常数 R_T(T_g) [J/(kg·K)]
% ------------------------------------------------------------------------
function RT = piecewise_RT(T_g)
    if T_g > 800 && T_g <= 1400
        RT = 287.6;
    elseif T_g > 1400 && T_g <= 1600
        RT = 288.0;
    elseif T_g > 1600
        RT = 288.6;
    else
        RT = 287.6;
    end
end

% ------------------------------------------------------------------------
% 函数：delta_cooling
% 用途：计算涡轮冷却空气引气系数 delta(T_g)
% ------------------------------------------------------------------------
function d = delta_cooling(T_g)
    d = 0.02 + (T_g - 1200) / 100 * 0.02;
    d = max(0.0, min(0.10, d));
end

% ------------------------------------------------------------------------
% 函数：engine_forward
% 用途：双涵道涡扇发动机热力循环前向模型
% 输入：theta (13×1 参数向量), cond (单工况结构体)
% 输出：y = [R_ud, C_ud]
% ------------------------------------------------------------------------
function [y, aux] = engine_forward(theta, cond)

    eta_k      = theta(1);
    eta_t      = theta(2);
    eta_m      = theta(3);
    eta_v      = theta(4);
    eta_tv     = theta(5);
    eta_c1     = theta(6);
    eta_c2     = theta(7);
    sigma_cc   = theta(8);
    sigma_kan  = theta(9);
    sigma_kask = theta(10);
    sigma_ks   = theta(11);
    eta_T_comb = theta(12);
    lambda_star = theta(13);

    T_H      = cond.T_H;
    M_flight = cond.M_flight;
    m        = cond.m;
    pi_k     = cond.pi_k;
    T_g      = cond.T_g;

    y   = [NaN, NaN];
    aux = struct();

    k_air = 1.4;
    R_air = 287.3;

    a_sound  = sqrt(k_air * R_air * T_H);
    V_flight = a_sound * M_flight;

    kT    = piecewise_kT(T_g);
    RT    = piecewise_RT(T_g);
    delta = delta_cooling(T_g);

    inner_v = 1 + V_flight^2 / (2 * (k_air / (k_air - 1)) * R_air * T_H);
    if inner_v <= 0
        return;
    end
    tau_v = inner_v ^ (k_air / (k_air - 1));

    T_B = T_H * inner_v;

    pk_exp = (k_air - 1) / k_air;
    T_k = T_B * (1 + (pi_k ^ pk_exp - 1) / eta_k);
    if ~isfinite(T_k) || T_k <= 0
        return;
    end

    g_T = 3.293e-5 * T_g - 2.84e-5 * T_k - 0.004814;
    if g_T <= 0 || ~isfinite(g_T)
        return;
    end

    A_comp = (k_air / (k_air - 1)) * R_air * T_B * (pi_k ^ pk_exp - 1);
    denom_factor = (kT / (kT - 1)) * RT * T_g;

    num_lam   = 1 - A_comp / (denom_factor * eta_k);
    denom_lam = 1 - A_comp / (denom_factor * eta_k * eta_t);

    if abs(denom_lam) < 1e-10 || ~isfinite(num_lam) || ~isfinite(denom_lam)
        return;
    end
    lambda_heat = num_lam / denom_lam;

    if ~isfinite(lambda_heat) || lambda_heat <= 0
        return;
    end

    sigma_bx = sigma_cc * sigma_kan;

    exp_T = (kT - 1) / kT;
    press_prod = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
    if press_prod <= 0
        return;
    end
    turb_bracket = 1 - (1 / press_prod) ^ exp_T;
    turb_work = (kT / (kT - 1)) * RT * T_g * turb_bracket;

    comp_coeff = (1 + g_T) * eta_k * eta_t * eta_m * (1 - delta);
    if abs(comp_coeff) < 1e-10
        return;
    end
    comp_work = (k_air / (k_air - 1)) * R_air * T_B * (pi_k ^ pk_exp - 1) / comp_coeff;

    L_sv = lambda_heat * turb_work - comp_work;
    if ~isfinite(L_sv) || L_sv <= 0
        return;
    end

    num_xpc   = 1 + (m * V_flight^2) / (2 * L_sv * eta_tv * eta_v * eta_c2);
    denom_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda_star);

    if abs(denom_xpc) < 1e-10 || ~isfinite(num_xpc)
        return;
    end
    x_pc = num_xpc / denom_xpc;

    if x_pc <= 0 || x_pc >= 1 || ~isfinite(x_pc)
        return;
    end

    inner_c1 = 2 * eta_c1 * lambda_star * x_pc * L_sv;
    if inner_c1 < 0
        return;
    end
    inner_c2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 + V_flight^2;
    if inner_c2 < 0
        return;
    end

    R_ud = 1 / (1 + m) * ((1 + g_T) * sqrt(inner_c1) - V_flight) ...
         + m / (1 + m) * (sqrt(inner_c2) - V_flight);

    if ~isfinite(R_ud)
        return;
    end

    denom_C = R_ud * (1 + m);
    if abs(denom_C) < 1e-10 || ~isfinite(denom_C)
        return;
    end
    C_ud = 3600 * g_T * (1 - delta) / denom_C;

    if ~isfinite(C_ud) || C_ud <= 0
        return;
    end

    y = [R_ud, C_ud];

    aux.T_B         = T_B;
    aux.T_k         = T_k;
    aux.g_T         = g_T;
    aux.lambda_heat = lambda_heat;
    aux.L_sv        = L_sv;
    aux.x_pc        = x_pc;
    aux.delta       = delta;
    aux.kT          = kT;
    aux.RT          = RT;
    aux.tau_v       = tau_v;
    aux.sigma_bx    = sigma_bx;
end

% ------------------------------------------------------------------------
% 函数：generate_virtual_data
% 用途：生成虚拟试验数据（真值 + 噪声），单工况版本
% ------------------------------------------------------------------------
function data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed)

    rng(rng_seed);

    y = engine_forward(theta_true, cond);
    if any(isnan(y))
        error('真值参数在该工况下前向模型返回 NaN，请检查参数设置');
    end

    R_true = y(1);
    C_true = y(2);

    sigma_R = noise_level_R * abs(R_true);
    sigma_C = noise_level_C * abs(C_true);

    R_obs = R_true + sigma_R * randn;
    C_obs = C_true + sigma_C * randn;

    data.R_true  = R_true;
    data.C_true  = C_true;
    data.R_obs   = R_obs;
    data.C_obs   = C_obs;
    data.sigma_R = sigma_R;
    data.sigma_C = sigma_C;
end

% ------------------------------------------------------------------------
% 函数：log_prior
% 用途：计算对数先验（均匀先验）
% ------------------------------------------------------------------------
function lp = log_prior(theta, lb, ub)
    if all(theta >= lb) && all(theta <= ub)
        lp = 0;
    else
        lp = -Inf;
    end
end

% ------------------------------------------------------------------------
% 函数：log_likelihood
% 用途：计算单工况对数似然（独立高斯误差）
% ------------------------------------------------------------------------
function ll = log_likelihood(theta, data, cond)

    y = engine_forward(theta, cond);

    if any(isnan(y)) || any(isinf(y))
        ll = -Inf;
        return;
    end

    R_pred = y(1);
    C_pred = y(2);

    sR = data.sigma_R;
    sC = data.sigma_C;

    if sR <= 0 || sC <= 0
        ll = -Inf;
        return;
    end

    % 对数似然（高斯）
    ll = - 0.5 * ( ((data.R_obs - R_pred) / sR)^2 ...
                 + ((data.C_obs - C_pred) / sC)^2 ...
                 + log(2 * pi * sR^2) ...
                 + log(2 * pi * sC^2) );
end

% ------------------------------------------------------------------------
% 函数：log_posterior_fn
% 用途：对数后验 = 对数先验 + 对数似然
% ------------------------------------------------------------------------
function lpost = log_posterior_fn(theta, data, cond, lb, ub)
    lp = log_prior(theta, lb, ub);
    if isinf(lp)
        lpost = -Inf;
        return;
    end
    ll    = log_likelihood(theta, data, cond);
    lpost = lp + ll;
end

% ------------------------------------------------------------------------
% 函数：reflect_to_unit_interval
% 用途：将提议步骤中越界的归一化坐标 z 做反射处理
% ------------------------------------------------------------------------
function z_ref = reflect_to_unit_interval(z)
    z_ref = z;
    for j = 1:length(z)
        zj = z(j);
        max_iter = 20;
        iter = 0;
        while (zj < 0 || zj > 1) && iter < max_iter
            if zj < 0
                zj = -zj;
            end
            if zj > 1
                zj = 2 - zj;
            end
            iter = iter + 1;
        end
        zj = max(0, min(1, zj));
        z_ref(j) = zj;
    end
end

% ------------------------------------------------------------------------
% 函数：run_mcmc
% 用途：随机游走 Metropolis-Hastings 采样（归一化参数空间）
% ------------------------------------------------------------------------
function results = run_mcmc(data, cond, lb, ub, theta0, opts)

    n_params  = length(lb);
    N         = opts.n_samples;
    burn_in   = opts.burn_in;
    prop_sd   = opts.proposal_sd;

    chain_full   = zeros(N, n_params);
    logpost_full = zeros(N, 1);

    z_curr     = (theta0 - lb) ./ (ub - lb);
    theta_curr = lb + z_curr .* (ub - lb);
    lpost_curr = log_posterior_fn(theta_curr, data, cond, lb, ub);

    if isinf(lpost_curr)
        error('初始点对数后验为 -Inf，请检查 theta0 是否在先验范围内');
    end

    n_accept = 0;
    accept_count_window = 0;
    window_start = 1;

    fprintf('进度: ');

    for s = 1:N

        if mod(s, round(N/10)) == 0
            fprintf('%d%% ', round(s/N*100));
        end

        z_prop = z_curr + prop_sd * randn(1, n_params);
        z_prop = reflect_to_unit_interval(z_prop);
        theta_prop = lb + z_prop .* (ub - lb);

        lpost_prop = log_posterior_fn(theta_prop, data, cond, lb, ub);

        log_alpha = lpost_prop - lpost_curr;
        if log(rand) < log_alpha
            z_curr     = z_prop;
            theta_curr = theta_prop;
            lpost_curr = lpost_prop;
            n_accept   = n_accept + 1;
            accept_count_window = accept_count_window + 1;
        end

        chain_full(s, :)   = theta_curr;
        logpost_full(s)    = lpost_curr;

        if s >= opts.adapt_start && s <= opts.adapt_end
            if mod(s - opts.adapt_start, opts.adapt_interval) == 0 && s > opts.adapt_start
                window_len = s - window_start + 1;
                local_rate = accept_count_window / window_len;

                if local_rate < opts.target_accept_low
                    prop_sd = prop_sd / opts.adapt_factor;
                elseif local_rate > opts.target_accept_high
                    prop_sd = prop_sd * opts.adapt_factor;
                end

                accept_count_window = 0;
                window_start = s + 1;
            end
        end
    end

    fprintf('\n');

    chain_post   = chain_full(burn_in+1:end, :);
    logpost_post = logpost_full(burn_in+1:end);

    theta_mean = mean(chain_post, 1);
    theta_std  = std(chain_post, 0, 1);

    [~, best_idx] = max(logpost_post);
    theta_map = chain_post(best_idx, :);

    theta_ci95 = zeros(n_params, 2);
    for j = 1:n_params
        theta_ci95(j, :) = quantile(chain_post(:, j), [0.025, 0.975]);
    end

    results.chain_full    = chain_full;
    results.chain_post    = chain_post;
    results.logpost_full  = logpost_full;
    results.logpost_post  = logpost_post;
    results.accept_rate   = n_accept / N;
    results.theta_mean    = theta_mean;
    results.theta_map     = theta_map;
    results.theta_ci95    = theta_ci95;
    results.theta_std     = theta_std;
    results.best_idx      = best_idx;
    results.final_prop_sd = prop_sd;
end

% ------------------------------------------------------------------------
% 函数：plot_results
% 用途：绘制所有可视化图形（单工况版本）
% ------------------------------------------------------------------------
function plot_results(results, data, cond, param_names, theta_true, lb, ub, ...
                      R_pred_mean, C_pred_mean, R_pred_map, C_pred_map)

    chain_post = results.chain_post;
    chain_full = results.chain_full;
    n_params   = size(chain_post, 2);

    % --- 图1：链轨迹图（Trace Plots）---
    figure('Name', 'MCMC链轨迹图', 'Position', [50, 50, 1400, 900]);
    n_rows = ceil(n_params / 3);
    for j = 1:n_params
        subplot(n_rows, 3, j);
        plot(chain_full(:, j), 'Color', [0.3 0.5 0.8], 'LineWidth', 0.5);
        hold on;
        xline(size(chain_full,1) - size(chain_post,1), 'r--', 'LineWidth', 1.5);
        yline(theta_true(j), 'g-', 'LineWidth', 1.5);
        xlabel('迭代步');
        ylabel(param_names{j});
        title(param_names{j});
        legend({'链', 'burn-in边界', '真值'}, 'FontSize', 6, 'Location', 'best');
        set(gca, 'FontSize', 8);
    end
    sgtitle(sprintf('MCMC 链轨迹图 | 单工况: pi\_k=%d, T\_g=%.0fK, M=%.1f, T\_H=%.0fK', ...
        cond.pi_k, cond.T_g, cond.M_flight, cond.T_H), 'FontSize', 11);

    % --- 图2：边缘后验直方图 ---
    figure('Name', '边缘后验分布', 'Position', [50, 50, 1400, 900]);
    for j = 1:n_params
        subplot(n_rows, 3, j);
        histogram(chain_post(:, j), 40, 'FaceColor', [0.4 0.6 0.8], ...
                  'EdgeColor', 'none', 'Normalization', 'pdf');
        hold on;
        xline(theta_true(j),          'g-',  'LineWidth', 2);
        xline(results.theta_mean(j),  'b--', 'LineWidth', 2);
        xline(results.theta_map(j),   'r:',  'LineWidth', 2);
        xlabel(param_names{j});
        ylabel('PDF');
        title(param_names{j});
        legend({'后验', '真值', '均值', 'MAP'}, 'FontSize', 6, 'Location', 'best');
        set(gca, 'FontSize', 8);
    end
    sgtitle('参数边缘后验分布（单工况）', 'FontSize', 12);

    % --- 图3：参数相关系数热图 ---
    figure('Name', '参数相关性热图', 'Position', [100, 100, 700, 600]);
    C_corr = corrcoef(chain_post);
    imagesc(C_corr, [-1, 1]);
    colorbar;
    colormap(redblue_cmap());
    set(gca, 'XTick', 1:n_params, 'XTickLabel', param_names, ...
             'YTick', 1:n_params, 'YTickLabel', param_names, ...
             'XTickLabelRotation', 45, 'FontSize', 8);
    title('后验样本参数相关系数矩阵（单工况）');
    for r = 1:n_params
        for c = 1:n_params
            text(c, r, sprintf('%.2f', C_corr(r,c)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 6, 'Color', 'k');
        end
    end

    % --- 图4：观测与预测对比（单工况柱状图）---
    figure('Name', '观测与预测对比', 'Position', [150, 150, 800, 400]);

    subplot(1, 2, 1);
    vals_R = [data.R_obs, R_pred_mean, R_pred_map];
    b = bar(vals_R, 'FaceColor', 'flat');
    b.CData = [0.2 0.2 0.2; 0.3 0.5 0.8; 0.8 0.3 0.3];
    hold on;
    errorbar(1, data.R_obs, 2*data.sigma_R, 'k', 'LineWidth', 2);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'观测值', '后验均值', 'MAP'}, 'FontSize', 10);
    ylabel('比推力 R_{ud} [N·s/kg]');
    title('比推力对比');
    grid on;

    subplot(1, 2, 2);
    vals_C = [data.C_obs, C_pred_mean, C_pred_map];
    b2 = bar(vals_C, 'FaceColor', 'flat');
    b2.CData = [0.2 0.2 0.2; 0.3 0.5 0.8; 0.8 0.3 0.3];
    hold on;
    errorbar(1, data.C_obs, 2*data.sigma_C, 'k', 'LineWidth', 2);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'观测值', '后验均值', 'MAP'}, 'FontSize', 10);
    ylabel('比油耗 C_{ud} [kg/(N·h)]');
    title('比油耗对比');
    grid on;

    sgtitle(sprintf('单工况观测值与预测值对比 | pi\_k=%d, T\_g=%.0fK', ...
        cond.pi_k, cond.T_g), 'FontSize', 12);

    % --- 图5：可辨识性指标（后验标准差 / 先验标准差）---
    figure('Name', '参数可辨识性', 'Position', [200, 200, 800, 400]);
    prior_std = (ub - lb) / sqrt(12);
    shrinkage = results.theta_std ./ prior_std;

    bar(1:n_params, shrinkage, 'FaceColor', [0.3 0.7 0.5], 'EdgeColor', 'none');
    hold on;
    yline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', '无信息基准线');
    yline(0.5, 'b:',  'LineWidth', 1.5, 'DisplayName', '50%收缩阈值');
    set(gca, 'XTick', 1:n_params, 'XTickLabel', param_names, ...
             'XTickLabelRotation', 45, 'FontSize', 9);
    ylabel('后验标准差 / 先验标准差');
    title('参数可辨识性指标（单工况）—— 比值越小 → 可辨识性越强');
    legend('Location', 'northeast', 'FontSize', 9);
    ylim([0, 1.2]);
    grid on;

    for j = 1:n_params
        text(j, shrinkage(j) + 0.02, sprintf('%.2f', shrinkage(j)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end

    fprintf('\n图形已绘制完成（共5张）\n');
end

% ------------------------------------------------------------------------
% 辅助函数：redblue_cmap
% ------------------------------------------------------------------------
function cmap = redblue_cmap()
    n = 64;
    r2 = [linspace(0.17, 1, n/2), linspace(1, 0.84, n/2)];
    g2 = [linspace(0.51, 1, n/2), linspace(1, 0.19, n/2)];
    b2 = [linspace(0.73, 1, n/2), linspace(1, 0.15, n/2)];
    cmap = [r2', g2', b2'];
end
