%% ========================================================================
% bayesian_inversion_multicondition.m
%
% 【程序说明】多工况联合贝叶斯反演 —— 双涵道涡扇发动机热力循环参数识别
%
% 本程序的核心任务：
%   - 观测量：多个飞行工况下的"比推力 R_ud"与"比油耗 C_ud"
%   - 未知量：13 个发动机修正参数（效率、总压恢复系数等）
%   - 方法：通过先验分布（均匀分布）+ 多工况联合高斯似然函数，
%           构造参数的后验分布，再用 Metropolis-Hastings MCMC 采样
%   - 目标：不是给出一个单一的"最优值"，而是得到参数的后验分布，
%           用于分析参数可辨识性、参数间相关性以及不确定性量化
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
%   轻量模式（默认）：8 工况，12000 次采样，3000 burn-in
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
%  第二部分：构造多工况（8个工况）
%  工况变量：T_H(环境温度), M_flight(飞行马赫数), m(涵道比), pi_k(压比), T_g(涡前温度)
%% ========================================================================

% conds(i) 结构体，每个工况定义已知输入
conds = struct('T_H', {}, 'M_flight', {}, 'm', {}, 'pi_k', {}, 'T_g', {});

% 工况 1：起飞（低速、高推力、高温）
conds(1).T_H       = 288;
conds(1).M_flight  = 0.0;
conds(1).m         = 8;
conds(1).pi_k      = 25;
conds(1).T_g       = 1650;

% 工况 2：低速爬升
conds(2).T_H       = 280;
conds(2).M_flight  = 0.3;
conds(2).m         = 9;
conds(2).pi_k      = 22;
conds(2).T_g       = 1550;

% 工况 3：中速爬升
conds(3).T_H       = 265;
conds(3).M_flight  = 0.5;
conds(3).m         = 10;
conds(3).pi_k      = 28;
conds(3).T_g       = 1500;

% 工况 4：巡航低速（高空）
conds(4).T_H       = 250;
conds(4).M_flight  = 0.7;
conds(4).m         = 11;
conds(4).pi_k      = 30;
conds(4).T_g       = 1450;

% 工况 5：巡航标准（主工况）
conds(5).T_H       = 245;
conds(5).M_flight  = 0.8;
conds(5).m         = 10;
conds(5).pi_k      = 32;
conds(5).T_g       = 1400;

% 工况 6：巡航高速
conds(6).T_H       = 240;
conds(6).M_flight  = 0.85;
conds(6).m         = 9;
conds(6).pi_k      = 35;
conds(6).T_g       = 1480;

% 工况 7：下降低推力
conds(7).T_H       = 260;
conds(7).M_flight  = 0.6;
conds(7).m         = 12;
conds(7).pi_k      = 18;
conds(7).T_g       = 1300;

% 工况 8：地面慢车（静止、低温比）
conds(8).T_H       = 288;
conds(8).M_flight  = 0.0;
conds(8).m         = 8;
conds(8).pi_k      = 12;
conds(8).T_g       = 1250;

n_conds = length(conds);

fprintf('\n=== 工况列表 ===\n');
fprintf('%-4s  %-6s  %-8s  %-6s  %-6s  %-6s\n', ...
    '工况','T_H','M_flight','m','pi_k','T_g');
for i = 1:n_conds
    fprintf('%-4d  %-6.1f  %-8.2f  %-6.0f  %-6.0f  %-6.0f\n', ...
        i, conds(i).T_H, conds(i).M_flight, conds(i).m, conds(i).pi_k, conds(i).T_g);
end

%% ========================================================================
%  第三部分：用真值参数验证前向模型，生成虚拟试验数据
%% ========================================================================

noise_level_R = 0.01;  % 比推力相对噪声标准差 1%
noise_level_C = 0.01;  % 比油耗相对噪声标准差 1%
rng_seed      = 2024;

data = generate_virtual_data(theta_true, conds, noise_level_R, noise_level_C, rng_seed);

fprintf('\n=== 虚拟试验数据 ===\n');
fprintf('%-4s  %-10s  %-10s  %-10s  %-10s\n', ...
    '工况','R_true','R_obs','C_true','C_obs');
for i = 1:n_conds
    fprintf('%-4d  %-10.4f  %-10.4f  %-10.6f  %-10.6f\n', ...
        i, data.R_true(i), data.R_obs(i), data.C_true(i), data.C_obs(i));
end

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
results = run_mcmc(data, conds, lb, ub, theta0, opts);
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
R_pred_mean = zeros(1, n_conds);
C_pred_mean = zeros(1, n_conds);
R_pred_map  = zeros(1, n_conds);
C_pred_map  = zeros(1, n_conds);

for i = 1:n_conds
    y_mean = engine_forward(results.theta_mean, conds(i));
    y_map  = engine_forward(results.theta_map,  conds(i));
    R_pred_mean(i) = y_mean(1);
    C_pred_mean(i) = y_mean(2);
    R_pred_map(i)  = y_map(1);
    C_pred_map(i)  = y_map(2);
end

fprintf('%-4s  %-10s  %-12s  %-12s  %-10s  %-12s  %-12s\n', ...
    '工况','R_obs','R_pred_mean','R_pred_map','C_obs','C_pred_mean','C_pred_map');
for i = 1:n_conds
    fprintf('%-4d  %-10.4f  %-12.4f  %-12.4f  %-10.6f  %-12.6f  %-12.6f\n', ...
        i, data.R_obs(i), R_pred_mean(i), R_pred_map(i), ...
        data.C_obs(i), C_pred_mean(i), C_pred_map(i));
end

%% ========================================================================
%  第八部分：绘图
%% ========================================================================

plot_results(results, data, conds, param_names, theta_true, lb, ub, ...
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
        % T_g <= 800，超出图片定义范围，外推用 1.33
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
% 按图片公式：delta = 0.02 + (T_g - 1200)/100 * 0.02
% 注：此处按图片最接近的数学结构实现，即线性插值形式
% ------------------------------------------------------------------------
function d = delta_cooling(T_g)
    d = 0.02 + (T_g - 1200) / 100 * 0.02;
    % 物理约束：引气系数应在合理范围内
    d = max(0.0, min(0.10, d));
end

% ------------------------------------------------------------------------
% 函数：engine_forward
% 用途：双涵道涡扇发动机热力循环前向模型
% 输入：theta (13×1 参数向量), cond (单工况结构体)
% 输出：y = [R_ud, C_ud]，aux 为中间量结构体
% ------------------------------------------------------------------------
function [y, aux] = engine_forward(theta, cond)

    % --- 解包参数 ---
    eta_k      = theta(1);   % 压气机效率
    eta_t      = theta(2);   % 涡轮效率
    eta_m      = theta(3);   % 机械效率
    eta_v      = theta(4);   % 风扇效率
    eta_tv     = theta(5);   % 风扇涡轮效率
    eta_c1     = theta(6);   % 一次喷管效率（内涵）
    eta_c2     = theta(7);   % 二次喷管效率（外涵）
    sigma_cc   = theta(8);   % 进气道总压恢复系数
    sigma_kan  = theta(9);   % 进气通道总压恢复系数
    sigma_kask = theta(10);  % 压气机级间总压恢复系数
    sigma_ks   = theta(11);  % 燃烧室总压恢复系数
    eta_T_comb = theta(12);  % 燃烧放热系数（避免与温度变量名冲突）
    lambda_star = theta(13); % 风扇涡轮热恢复系数

    % --- 解包工况 ---
    T_H      = cond.T_H;
    M_flight = cond.M_flight;
    m        = cond.m;
    pi_k     = cond.pi_k;
    T_g      = cond.T_g;

    % 默认返回 NaN（异常时使用）
    y   = [NaN, NaN];
    aux = struct();

    % --- (1) 基本常数 ---
    k_air = 1.4;
    R_air = 287.3;

    % --- (2) 飞行速度 [m/s] ---
    a_sound  = sqrt(k_air * R_air * T_H);
    V_flight = a_sound * M_flight;

    % --- (3) 分段函数 ---
    kT    = piecewise_kT(T_g);
    RT    = piecewise_RT(T_g);
    delta = delta_cooling(T_g);

    % --- (4) 进口总压比 tau_v（进气道总压恢复）---
    % 按图片：tau_v = (1 + V^2 / (2*(k/(k-1))*R*T_H))^(k/(k-1))
    % 注：此处按图片最接近的数学结构实现
    inner_v = 1 + V_flight^2 / (2 * (k_air / (k_air - 1)) * R_air * T_H);
    if inner_v <= 0
        return;
    end
    tau_v = inner_v ^ (k_air / (k_air - 1));

    % --- (5) 压气机入口总温 T_B [K] ---
    % 按图片：T_B = T_H * (1 + V^2 / (2*(k/(k-1))*R*T_H))
    T_B = T_H * inner_v;

    % --- (6) 压气机出口温度 T_k [K] ---
    % T_k = T_B * (1 + (pi_k^((k-1)/k) - 1) / eta_k)
    pk_exp = (k_air - 1) / k_air;
    T_k = T_B * (1 + (pi_k ^ pk_exp - 1) / eta_k);
    if ~isfinite(T_k) || T_k <= 0
        return;
    end

    % --- (7) 相对耗油量 g_T ---
    % g_T = 3.293e-5 * T_g - 2.84e-5 * T_k - 0.004814
    g_T = 3.293e-5 * T_g - 2.84e-5 * T_k - 0.004814;
    if g_T <= 0 || ~isfinite(g_T)
        return;
    end

    % --- (8) 热恢复系数 lambda_heat(pi_k, T_g) ---
    % 按图片中分式结构实现：
    % 分子中的压气机功：
    %   W_k = (k/(k-1)) * R_air * T_B * (pi_k^((k-1)/k) - 1)
    % 分母需要用到涡轮膨胀功
    %
    % 按图片公式（lambda 分子分母）：
    %   公共项 A = (k/(k-1)) * R_air * T_B * (pi_k^((k-1)/k) - 1)
    %
    %   分子 = 1 - A / ( (kT/(kT-1)) * RT * T_g * eta_k )  [注：上分支]
    %   分母 = 1 - A / ( (kT/(kT-1)) * RT * T_g * eta_k * eta_t )
    %
    % 注：此处按图片最接近的数学结构实现
    %     图片中 lambda 分子分母的区别在于分母多了一个 eta_t

    A_comp = (k_air / (k_air - 1)) * R_air * T_B * (pi_k ^ pk_exp - 1);

    denom_factor = (kT / (kT - 1)) * RT * T_g;

    % 分子（不含 eta_t）
    num_lam   = 1 - A_comp / (denom_factor * eta_k);
    % 分母（含 eta_t）
    denom_lam = 1 - A_comp / (denom_factor * eta_k * eta_t);

    if abs(denom_lam) < 1e-10 || ~isfinite(num_lam) || ~isfinite(denom_lam)
        return;
    end
    lambda_heat = num_lam / denom_lam;

    if ~isfinite(lambda_heat) || lambda_heat <= 0
        return;
    end

    % --- (9) 进口总压恢复系数 sigma_bx ---
    sigma_bx = sigma_cc * sigma_kan;

    % --- (10) 单位自由能 L_sv [J/kg] ---
    % 按图片 L_cb 公式实现，命名为 L_sv
    % L_sv = lambda_heat * [ (kT/(kT-1)) * RT * T_g *
    %           (1 - (1/(tau_v*sigma_bx*pi_k*sigma_kask*sigma_ks))^((kT-1)/kT)) ]
    %        - (k/(k-1)) * R_air * T_B * (pi_k^((k-1)/k) - 1)
    %          * 1 / ((1 + g_T) * eta_k * eta_t * eta_m * (1 - delta))

    % 膨胀比指数
    exp_T = (kT - 1) / kT;

    % 涡轮膨胀压比乘积
    press_prod = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
    if press_prod <= 0
        return;
    end
    turb_bracket = 1 - (1 / press_prod) ^ exp_T;

    % 涡轮做功项
    turb_work = (kT / (kT - 1)) * RT * T_g * turb_bracket;

    % 压气机耗功项（含损失修正）
    comp_coeff = (1 + g_T) * eta_k * eta_t * eta_m * (1 - delta);
    if abs(comp_coeff) < 1e-10
        return;
    end
    comp_work = (k_air / (k_air - 1)) * R_air * T_B * (pi_k ^ pk_exp - 1) / comp_coeff;

    L_sv = lambda_heat * turb_work - comp_work;
    if ~isfinite(L_sv) || L_sv <= 0
        return;
    end

    % --- (11) 最优自由能分配系数 x_pc ---
    % 按图片公式：
    % x_pc = (1 + m*V^2 / (2*L_sv*eta_tv*eta_v*eta_c2))
    %        / (1 + m*eta_tv*eta_v*eta_c2 / (eta_c1*lambda_star))
    %
    % 注：此处按图片最接近的数学结构实现（分子分母分别含 m 项）

    num_xpc   = 1 + (m * V_flight^2) / (2 * L_sv * eta_tv * eta_v * eta_c2);
    denom_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda_star);

    if abs(denom_xpc) < 1e-10 || ~isfinite(num_xpc)
        return;
    end
    x_pc = num_xpc / denom_xpc;

    % 物理约束：0 < x_pc < 1
    if x_pc <= 0 || x_pc >= 1 || ~isfinite(x_pc)
        return;
    end

    % --- (12) 最优比推力 R_ud [N·s/kg] ---
    % 按图片公式（内涵 + 外涵两部分）：
    % R_ud = 1/(1+m) * [(1+g_T)*sqrt(2*eta_c1*lambda_star*x_pc*L_sv) - V]
    %      + m/(1+m) * [sqrt(2*(1-x_pc)/m * L_sv*eta_tv*eta_v*eta_c2 + V^2) - V]
    %
    % 注：此处按图片最接近的数学结构实现

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

    % --- (13) 比油耗 C_ud [kg/(N·h)] ---
    % C_ud = 3600 * g_T * (1 - delta) / (R_ud * (1 + m))
    denom_C = R_ud * (1 + m);
    if abs(denom_C) < 1e-10 || ~isfinite(denom_C)
        return;
    end
    C_ud = 3600 * g_T * (1 - delta) / denom_C;

    if ~isfinite(C_ud) || C_ud <= 0
        return;
    end

    % --- 输出 ---
    y = [R_ud, C_ud];

    % 中间量输出
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
% 用途：生成虚拟试验数据（真值 + 噪声）
% ------------------------------------------------------------------------
function data = generate_virtual_data(theta_true, conds, noise_level_R, noise_level_C, rng_seed)

    rng(rng_seed);
    n = length(conds);

    R_true  = zeros(1, n);
    C_true  = zeros(1, n);
    R_obs   = zeros(1, n);
    C_obs   = zeros(1, n);
    sigma_R = zeros(1, n);
    sigma_C = zeros(1, n);

    for i = 1:n
        y = engine_forward(theta_true, conds(i));
        if any(isnan(y))
            error('真值参数在工况 %d 下前向模型返回 NaN，请检查参数设置', i);
        end
        R_true(i) = y(1);
        C_true(i) = y(2);

        sigma_R(i) = noise_level_R * abs(R_true(i));
        sigma_C(i) = noise_level_C * abs(C_true(i));

        R_obs(i) = R_true(i) + sigma_R(i) * randn;
        C_obs(i) = C_true(i) + sigma_C(i) * randn;
    end

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
        lp = 0;  % log(constant) = 0
    else
        lp = -Inf;
    end
end

% ------------------------------------------------------------------------
% 函数：log_likelihood
% 用途：计算多工况联合对数似然（独立高斯误差）
% ------------------------------------------------------------------------
function ll = log_likelihood(theta, data, conds)

    n = length(conds);
    ll = 0;

    for i = 1:n
        y = engine_forward(theta, conds(i));

        % 前向模型异常检查
        if any(isnan(y)) || any(isinf(y))
            ll = -Inf;
            return;
        end

        R_pred = y(1);
        C_pred = y(2);

        sR = data.sigma_R(i);
        sC = data.sigma_C(i);

        if sR <= 0 || sC <= 0
            ll = -Inf;
            return;
        end

        % 对数似然累加（高斯）
        ll = ll - 0.5 * ( ((data.R_obs(i) - R_pred) / sR)^2 ...
                        + ((data.C_obs(i) - C_pred) / sC)^2 ...
                        + log(2 * pi * sR^2) ...
                        + log(2 * pi * sC^2) );
    end
end

% ------------------------------------------------------------------------
% 函数：log_posterior_fn
% 用途：对数后验 = 对数先验 + 对数似然
% ------------------------------------------------------------------------
function lpost = log_posterior_fn(theta, data, conds, lb, ub)
    lp = log_prior(theta, lb, ub);
    if isinf(lp)
        lpost = -Inf;
        return;
    end
    ll    = log_likelihood(theta, data, conds);
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
        % 反射直到落入 [0,1]
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
        % 保险：若仍越界则裁剪
        zj = max(0, min(1, zj));
        z_ref(j) = zj;
    end
end

% ------------------------------------------------------------------------
% 函数：run_mcmc
% 用途：随机游走 Metropolis-Hastings 采样（归一化参数空间）
% ------------------------------------------------------------------------
function results = run_mcmc(data, conds, lb, ub, theta0, opts)

    n_params  = length(lb);
    N         = opts.n_samples;
    burn_in   = opts.burn_in;
    prop_sd   = opts.proposal_sd;

    % 初始化链
    chain_full   = zeros(N, n_params);
    logpost_full = zeros(N, 1);

    % 初始点（归一化）
    z_curr     = (theta0 - lb) ./ (ub - lb);
    theta_curr = lb + z_curr .* (ub - lb);
    lpost_curr = log_posterior_fn(theta_curr, data, conds, lb, ub);

    if isinf(lpost_curr)
        error('初始点对数后验为 -Inf，请检查 theta0 是否在先验范围内');
    end

    n_accept = 0;
    accept_count_window = 0;
    window_start = 1;

    fprintf('进度: ');

    for s = 1:N

        % 打印进度
        if mod(s, round(N/10)) == 0
            fprintf('%d%% ', round(s/N*100));
        end

        % --- 提议步骤（归一化空间随机游走）---
        z_prop = z_curr + prop_sd * randn(1, n_params);

        % 反射边界处理
        z_prop = reflect_to_unit_interval(z_prop);

        % 映射回参数空间
        theta_prop = lb + z_prop .* (ub - lb);

        % --- 计算对数后验 ---
        lpost_prop = log_posterior_fn(theta_prop, data, conds, lb, ub);

        % --- M-H 接受准则 ---
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

        % --- 自适应调节提议步长 ---
        if s >= opts.adapt_start && s <= opts.adapt_end
            if mod(s - opts.adapt_start, opts.adapt_interval) == 0 && s > opts.adapt_start
                window_len = s - window_start + 1;
                local_rate = accept_count_window / window_len;

                if local_rate < opts.target_accept_low
                    prop_sd = prop_sd / opts.adapt_factor;
                elseif local_rate > opts.target_accept_high
                    prop_sd = prop_sd * opts.adapt_factor;
                end

                % 重置窗口
                accept_count_window = 0;
                window_start = s + 1;
            end
        end
    end

    fprintf('\n');

    % --- 提取后验样本（去掉 burn-in）---
    chain_post   = chain_full(burn_in+1:end, :);
    logpost_post = logpost_full(burn_in+1:end);

    % --- 计算后验统计 ---
    theta_mean = mean(chain_post, 1);
    theta_std  = std(chain_post, 0, 1);

    % MAP（后验样本中对数后验最大的点）
    [~, best_idx] = max(logpost_post);
    theta_map = chain_post(best_idx, :);

    % 95% 置信区间
    theta_ci95 = zeros(n_params, 2);
    for j = 1:n_params
        theta_ci95(j, :) = quantile(chain_post(:, j), [0.025, 0.975]);
    end

    % --- 输出 ---
    results.chain_full   = chain_full;
    results.chain_post   = chain_post;
    results.logpost_full = logpost_full;
    results.logpost_post = logpost_post;
    results.accept_rate  = n_accept / N;
    results.theta_mean   = theta_mean;
    results.theta_map    = theta_map;
    results.theta_ci95   = theta_ci95;
    results.theta_std    = theta_std;
    results.best_idx     = best_idx;
    results.final_prop_sd = prop_sd;
end

% ------------------------------------------------------------------------
% 函数：plot_results
% 用途：绘制所有可视化图形
% ------------------------------------------------------------------------
function plot_results(results, data, conds, param_names, theta_true, lb, ub, ...
                      R_pred_mean, C_pred_mean, R_pred_map, C_pred_map)

    chain_post = results.chain_post;
    chain_full = results.chain_full;
    n_params   = size(chain_post, 2);
    n_conds    = length(conds);

    % --- 图1：链轨迹图（Trace Plots）---
    fig1 = figure('Name', 'MCMC链轨迹图', 'Position', [50, 50, 1400, 900]);
    n_rows = ceil(n_params / 3);
    for j = 1:n_params
        subplot(n_rows, 3, j);
        plot(chain_full(:, j), 'Color', [0.3 0.5 0.8], 'LineWidth', 0.5);
        hold on;
        xline(size(results.chain_full,1) - size(chain_post,1), 'r--', 'LineWidth', 1.5);
        yline(theta_true(j), 'g-', 'LineWidth', 1.5);
        xlabel('迭代步');
        ylabel(param_names{j});
        title(param_names{j});
        legend({'链', 'burn-in边界', '真值'}, 'FontSize', 6, 'Location', 'best');
        set(gca, 'FontSize', 8);
    end
    sgtitle('MCMC 链轨迹图（红线=burn-in边界，绿线=真值）', 'FontSize', 12);

    % --- 图2：边缘后验直方图 ---
    fig2 = figure('Name', '边缘后验分布', 'Position', [50, 50, 1400, 900]);
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
    sgtitle('参数边缘后验分布', 'FontSize', 12);

    % --- 图3：参数相关系数热图 ---
    fig3 = figure('Name', '参数相关性热图', 'Position', [100, 100, 700, 600]);
    C_corr = corrcoef(chain_post);
    imagesc(C_corr, [-1, 1]);
    colorbar;
    colormap(redblue_cmap());
    set(gca, 'XTick', 1:n_params, 'XTickLabel', param_names, ...
             'YTick', 1:n_params, 'YTickLabel', param_names, ...
             'XTickLabelRotation', 45, 'FontSize', 8);
    title('后验样本参数相关系数矩阵');
    % 在格中标注数值
    for r = 1:n_params
        for c = 1:n_params
            text(c, r, sprintf('%.2f', C_corr(r,c)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 6, 'Color', 'k');
        end
    end

    % --- 图4：观测与预测对比图 ---
    fig4 = figure('Name', '观测与预测对比', 'Position', [150, 150, 1000, 450]);

    subplot(1, 2, 1);
    idx = 1:n_conds;
    errorbar(idx, data.R_obs, 2*data.sigma_R, 'ko', 'LineWidth', 1.5, ...
             'MarkerFaceColor', 'k', 'DisplayName', '观测值±2\sigma');
    hold on;
    plot(idx, R_pred_mean, 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '后验均值预测');
    plot(idx, R_pred_map,  'r^--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'MAP预测');
    xlabel('工况编号');
    ylabel('比推力 R_{ud} [N·s/kg]');
    title('比推力：观测 vs 预测');
    legend('Location', 'best', 'FontSize', 9);
    grid on;

    subplot(1, 2, 2);
    errorbar(idx, data.C_obs, 2*data.sigma_C, 'ko', 'LineWidth', 1.5, ...
             'MarkerFaceColor', 'k', 'DisplayName', '观测值±2\sigma');
    hold on;
    plot(idx, C_pred_mean, 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '后验均值预测');
    plot(idx, C_pred_map,  'r^--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'MAP预测');
    xlabel('工况编号');
    ylabel('比油耗 C_{ud} [kg/(N·h)]');
    title('比油耗：观测 vs 预测');
    legend('Location', 'best', 'FontSize', 9);
    grid on;

    sgtitle('各工况观测值与预测值对比', 'FontSize', 13);

    % --- 图5：可辨识性指标（后验标准差 / 先验标准差）---
    fig5 = figure('Name', '参数可辨识性', 'Position', [200, 200, 800, 400]);
    prior_std = (ub - lb) / sqrt(12);  % 均匀分布标准差
    shrinkage = results.theta_std ./ prior_std;

    bar(1:n_params, shrinkage, 'FaceColor', [0.3 0.7 0.5], 'EdgeColor', 'none');
    hold on;
    yline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', '无信息基准线');
    yline(0.5, 'b:',  'LineWidth', 1.5, 'DisplayName', '50%收缩阈值');
    set(gca, 'XTick', 1:n_params, 'XTickLabel', param_names, ...
             'XTickLabelRotation', 45, 'FontSize', 9);
    ylabel('后验标准差 / 先验标准差');
    title('参数可辨识性指标（比值越小 → 后验收缩越显著 → 可辨识性越强）');
    legend('Location', 'northeast', 'FontSize', 9);
    ylim([0, 1.2]);
    grid on;

    % 在柱上标注数值
    for j = 1:n_params
        text(j, shrinkage(j) + 0.02, sprintf('%.2f', shrinkage(j)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end

    fprintf('\n图形已绘制完成（共5张）\n');
end

% ------------------------------------------------------------------------
% 辅助函数：redblue_cmap
% 用途：生成红-白-蓝颜色图（用于相关性热图）
% ------------------------------------------------------------------------
function cmap = redblue_cmap()
    n = 64;
    r = [linspace(0.7, 1, n/2), linspace(1, 1, n/2)];
    g = [linspace(0.2, 1, n/2), linspace(1, 0.2, n/2)];
    b = [linspace(1, 1, n/2), linspace(1, 0.7, n/2)];
    % 更标准：蓝->白->红
    r2 = [linspace(0.17, 1, n/2), linspace(1, 0.84, n/2)];
    g2 = [linspace(0.51, 1, n/2), linspace(1, 0.19, n/2)];
    b2 = [linspace(0.73, 1, n/2), linspace(1, 0.15, n/2)];
    cmap = [r2', g2', b2'];
end
