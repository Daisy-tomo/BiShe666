%% =========================================================================
% bayesian_inversion_single_condition.m
%
% 【研究说明】
% 1. 本程序研究的是单工况条件下涡扇发动机热力循环的贝叶斯后验反演；
% 2. 观测量只有比推力 R_ud 和比油耗 C_ud 两个标量；
% 3. 未知量是 13 个修正参数（效率、总压恢复系数、热恢复系数等）；
% 4. 由于观测少于参数数目（2 < 13），问题是严重欠定的；
% 5. 因此程序目标不是给出唯一确定解，而是：
%       - 分析参数后验分布的收缩程度（可辨识性）；
%       - 分析参数间的后验相关性；
%       - 定量描述参数的不确定性区间；
% 6. 程序先用虚拟试验数据验证反演流程，保证前向模型与反演代码正确联通。
%
% 使用方法：直接在 MATLAB 中运行本文件即可。
% 建议环境：MATLAB R2019b 及以上版本。
% 运行时间：普通笔记本电脑约 1~3 分钟（10000 次采样）。
% =========================================================================

clear; clc; close all;
rng(42);  % 固定随机种子，保证结果可重复

%% =========================================================================
% 第一部分：参数定义
% =========================================================================

% 13 个未知修正参数的名称
param_names = {
    'eta\_k',       ...  % 1  压气机绝热效率
    'eta\_t',       ...  % 2  涡轮绝热效率
    'eta\_m',       ...  % 3  机械效率
    'eta\_v',       ...  % 4  风扇效率
    'eta\_tv',      ...  % 5  风扇涡轮效率
    'eta\_c1',      ...  % 6  一次喷管效率
    'eta\_c2',      ...  % 7  二次喷管效率
    'sigma\_cc',    ...  % 8  进气道激波总压恢复系数
    'sigma\_kan',   ...  % 9  进气通道总压恢复系数
    'sigma\_kask',  ...  % 10 压气机级间总压恢复系数
    'sigma\_ks',    ...  % 11 燃烧室总压恢复系数
    'eta\_T',       ...  % 12 燃烧放热系数
    'lambda\_star'  ...  % 13 风扇涡轮热恢复系数（lambda'）
};

% 参数下界
lb = [0.84, 0.86, 0.980, 0.85, 0.90, 0.94, 0.92, 0.98, 0.98, 0.98, 0.94, 0.97, 1.02];

% 参数上界
ub = [0.92, 0.92, 0.995, 0.87, 0.92, 0.95, 0.94, 1.00, 0.99, 0.99, 0.96, 0.99, 1.04];

n_params = length(lb);

%% =========================================================================
% 第二部分：单工况设置
% =========================================================================
% 下面设置一组典型的起飞/低速飞行工况（飞行马赫数接近 0，低空起飞条件）
% T_g 选取在 1400~1600 K 区间，以激活中段分段函数

cond.T_H      = 288.15;   % [K]  环境温度（标准大气，海平面）
cond.M_flight = 0.0;      % [-]  飞行马赫数（起飞状态，地面静止）
cond.m        = 10.0;     % [-]  涵道比（大涵道比涡扇，典型值 8~13，取中间值）
cond.pi_k     = 30.0;     % [-]  压气机总压比（典型高涵道比涡扇，取中段值）
cond.T_g      = 1500.0;   % [K]  涡轮前总温（落在 1400<T_g<=1600 分段区间）

fprintf('========================================================\n');
fprintf('  涡扇发动机参数贝叶斯反演 - 单工况版本\n');
fprintf('========================================================\n');
fprintf('工况设置：\n');
fprintf('  T_H      = %.2f K\n', cond.T_H);
fprintf('  M_flight = %.2f\n',   cond.M_flight);
fprintf('  m (涵道比) = %.1f\n', cond.m);
fprintf('  pi_k     = %.1f\n',   cond.pi_k);
fprintf('  T_g      = %.1f K\n', cond.T_g);
fprintf('\n');

%% =========================================================================
% 第三部分：设置真值参数与初始值
% =========================================================================
% theta_true：程序内置的"真实"参数值，落在参数范围内部（不贴边界）
% 各参数取范围中点附近，加小扰动，确保不在边界上

theta_true = [
    0.88,    ...  % eta_k
    0.89,    ...  % eta_t
    0.988,   ...  % eta_m
    0.860,   ...  % eta_v
    0.910,   ...  % eta_tv
    0.945,   ...  % eta_c1
    0.930,   ...  % eta_c2
    0.990,   ...  % sigma_cc
    0.985,   ...  % sigma_kan
    0.985,   ...  % sigma_kask
    0.950,   ...  % sigma_ks
    0.980,   ...  % eta_T
    1.030    ...  % lambda_star
];

% 验证 theta_true 落在范围内
assert(all(theta_true >= lb) && all(theta_true <= ub), ...
    '错误：theta_true 超出参数范围！');

% theta0：MCMC 初始点，取先验中点（稳健起点）
theta0 = 0.5 * (lb + ub);

fprintf('真值参数 theta_true：\n');
for i = 1:n_params
    fprintf('  %-15s = %.4f  [范围: %.4f ~ %.4f]\n', ...
        param_names{i}, theta_true(i), lb(i), ub(i));
end
fprintf('\n初始值 theta0（先验中点）：\n');
for i = 1:n_params
    fprintf('  %-15s = %.4f\n', param_names{i}, theta0(i));
end
fprintf('\n');

%% =========================================================================
% 第四部分：验证前向模型
% =========================================================================
fprintf('--- 前向模型验证 ---\n');
[y_true, aux_true] = engine_forward(theta_true, cond);
if ~all(isfinite(y_true))
    error('前向模型在 theta_true 处输出非有限值，请检查参数设置！');
end
fprintf('theta_true 前向计算结果：\n');
fprintf('  R_ud_true = %.4f [N·s/kg]\n', y_true(1));
fprintf('  C_ud_true = %.6f [kg/(N·h)]\n', y_true(2));
fprintf('中间变量：\n');
fprintf('  T_B = %.2f K, T_k = %.2f K\n', aux_true.T_B, aux_true.T_k);
fprintf('  g_T = %.6f\n', aux_true.g_T);
fprintf('  lambda_heat = %.4f\n', aux_true.lambda_heat);
fprintf('  L_sv = %.2f J/kg\n', aux_true.L_sv);
fprintf('  x_pc = %.4f\n', aux_true.x_pc);
fprintf('  tau_v = %.4f\n', aux_true.tau_v);
fprintf('\n');

%% =========================================================================
% 第五部分：生成虚拟试验数据
% =========================================================================
noise_level_R = 0.01;   % R_ud 的相对噪声水平（1%）
noise_level_C = 0.01;   % C_ud 的相对噪声水平（1%）
rng_seed      = 123;    % 噪声随机种子

data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed);

fprintf('--- 虚拟试验数据 ---\n');
fprintf('  R_ud_true = %.4f,  R_ud_obs = %.4f,  sigma_R = %.4f\n', ...
    data.R_true, data.R_obs, data.sigma_R);
fprintf('  C_ud_true = %.6f,  C_ud_obs = %.6f,  sigma_C = %.6f\n', ...
    data.C_true, data.C_obs, data.sigma_C);
fprintf('\n');

%% =========================================================================
% 第六部分：设置 MCMC 参数
% =========================================================================
opts.n_samples          = 10000;  % 总采样步数（轻量配置）
opts.burn_in            = 2000;   % 烧入期步数
opts.proposal_sd        = 0.012;  % 归一化空间中的提案标准差
opts.adapt_start        = 500;    % 自适应开始步
opts.adapt_end          = 2500;   % 自适应结束步
opts.adapt_interval     = 100;    % 自适应调整间隔
opts.target_accept_low  = 0.18;   % 目标接受率下限
opts.target_accept_high = 0.35;   % 目标接受率上限
opts.print_interval     = 1000;   % 打印进度间隔

fprintf('--- MCMC 配置 ---\n');
fprintf('  总采样步数:   %d\n', opts.n_samples);
fprintf('  烧入期:       %d\n', opts.burn_in);
fprintf('  后验链长度:   %d\n', opts.n_samples - opts.burn_in);
fprintf('  初始提案 SD:  %.4f（归一化空间）\n', opts.proposal_sd);
fprintf('\n');

%% =========================================================================
% 第七部分：运行 MCMC
% =========================================================================
fprintf('开始 MCMC 采样...\n');
tic;
results = run_mcmc(data, cond, lb, ub, theta0, opts);
t_elapsed = toc;
fprintf('MCMC 完成，用时 %.1f 秒\n\n', t_elapsed);

%% =========================================================================
% 第八部分：输出后验统计
% =========================================================================
fprintf('========================================================\n');
fprintf('  后验统计汇总\n');
fprintf('========================================================\n');
fprintf('总体接受率: %.3f\n\n', results.accept_rate);

fprintf('%-15s  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
    '参数名', '真值', 'Post均值', 'MAP', 'CI95低', 'CI95高', '先验中点');
fprintf('%s\n', repmat('-', 1, 75));
for i = 1:n_params
    fprintf('%-15s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
        param_names{i}, ...
        theta_true(i), ...
        results.theta_mean(i), ...
        results.theta_map(i), ...
        results.theta_ci95(i,1), ...
        results.theta_ci95(i,2), ...
        theta0(i));
end
fprintf('\n');

%% =========================================================================
% 第九部分：后验预测
% =========================================================================
[y_mean, ~] = engine_forward(results.theta_mean, cond);
[y_map,  ~] = engine_forward(results.theta_map,  cond);

R_pred_mean = y_mean(1);  C_pred_mean = y_mean(2);
R_pred_map  = y_map(1);   C_pred_map  = y_map(2);

fprintf('--- 后验预测对比 ---\n');
fprintf('  %-20s  %10s  %10s\n', '量', 'R_ud [N·s/kg]', 'C_ud [kg/(N·h)]');
fprintf('  %-20s  %10.4f  %10.6f\n', '真值',            data.R_true,   data.C_true);
fprintf('  %-20s  %10.4f  %10.6f\n', '观测值（含噪声）', data.R_obs,    data.C_obs);
fprintf('  %-20s  %10.4f  %10.6f\n', '后验均值预测',     R_pred_mean,   C_pred_mean);
fprintf('  %-20s  %10.4f  %10.6f\n', 'MAP预测',          R_pred_map,    C_pred_map);
fprintf('\n');

% 相对误差
err_mean_R = abs(R_pred_mean - data.R_true) / abs(data.R_true) * 100;
err_mean_C = abs(C_pred_mean - data.C_true) / abs(data.C_true) * 100;
err_map_R  = abs(R_pred_map  - data.R_true) / abs(data.R_true) * 100;
err_map_C  = abs(C_pred_map  - data.C_true) / abs(data.C_true) * 100;
fprintf('  后验均值预测相对误差: R_ud=%.3f%%, C_ud=%.3f%%\n', err_mean_R, err_mean_C);
fprintf('  MAP预测相对误差:       R_ud=%.3f%%, C_ud=%.3f%%\n', err_map_R,  err_map_C);
fprintf('\n');

%% =========================================================================
% 第十部分：绘图
% =========================================================================
plot_results(results, theta_true, lb, ub, param_names);

fprintf('所有图形已生成。程序运行完毕。\n');

%% =========================================================================
%  ↓↓↓  以下为所有局部函数  ↓↓↓
%% =========================================================================

%% -------------------------------------------------------------------------
%  分段函数：k_T(T_g) — 燃气绝热指数
%% -------------------------------------------------------------------------
function kT = piecewise_kT(T_g)
% piecewise_kT  按图片公式精确实现燃气绝热指数分段函数
%   800 < T_g <= 1400  ->  kT = 1.33
%   1400 < T_g <= 1600 ->  kT = 1.30
%   T_g > 1600         ->  kT = 1.25
if T_g > 800 && T_g <= 1400
    kT = 1.33;
elseif T_g > 1400 && T_g <= 1600
    kT = 1.30;
elseif T_g > 1600
    kT = 1.25;
else
    % T_g <= 800，外推用最小值（保护性处理）
    kT = 1.33;
end
end

%% -------------------------------------------------------------------------
%  分段函数：R_T(T_g) — 燃气气体常数 [J/(kg·K)]
%% -------------------------------------------------------------------------
function RT = piecewise_RT(T_g)
% piecewise_RT  按图片公式精确实现燃气气体常数分段函数
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

%% -------------------------------------------------------------------------
%  分段函数：delta(T_g) — 涡轮冷却引气系数
%% -------------------------------------------------------------------------
function d = delta_cooling(T_g)
% delta_cooling  按图片公式实现冷却引气系数
%   delta(T_g) = 0.02 + (T_g - 1200)/100 * 0.02
%   此处按图片最接近的数学结构实现（线性插值形式）
d = 0.02 + (T_g - 1200) / 100 * 0.02;
% 物理约束：delta 不应为负，也不宜过大
d = max(0.0, min(d, 0.15));
end

%% -------------------------------------------------------------------------
%  前向模型
%% -------------------------------------------------------------------------
function [y, aux] = engine_forward(theta, cond)
% engine_forward  双涵道涡扇发动机热力循环前向模型
%
%   输入：
%     theta  - 13 个修正参数向量
%     cond   - 工况结构体（T_H, M_flight, m, pi_k, T_g）
%
%   输出：
%     y      - [R_ud, C_ud]，比推力和比油耗
%     aux    - 中间变量结构体
%
%   参数解包（严格按定义顺序）
eta_k      = theta(1);   % 压气机效率
eta_t      = theta(2);   % 涡轮效率
eta_m      = theta(3);   % 机械效率
eta_v      = theta(4);   % 风扇效率
eta_tv     = theta(5);   % 风扇涡轮效率
eta_c1     = theta(6);   % 一次喷管效率
eta_c2     = theta(7);   % 二次喷管效率
sigma_cc   = theta(8);   % 进气道激波总压恢复系数
sigma_kan  = theta(9);   % 进气通道总压恢复系数
sigma_kask = theta(10);  % 压气机级间总压恢复系数
sigma_ks   = theta(11);  % 燃烧室总压恢复系数
eta_T      = theta(12);  % 燃烧放热系数
lambda_star= theta(13);  % 风扇涡轮热恢复系数 lambda'

%   工况解包
T_H      = cond.T_H;
M_flight = cond.M_flight;
m        = cond.m;
pi_k     = cond.pi_k;
T_g      = cond.T_g;

% 默认返回无效值（在检测到问题时使用）
y   = [NaN, NaN];
aux = struct();

try
    %% ----------------------------------------------------------------
    % (1) 基本常数
    k_air = 1.4;       % 空气绝热指数
    R_air = 287.3;     % [J/(kg·K)] 空气气体常数

    % 声速
    a = sqrt(k_air * R_air * T_H);
    if ~isfinite(a) || a <= 0, return; end

    % 飞行速度
    V_flight = a * M_flight;   % 起飞时 M_flight=0 -> V_flight=0

    %% ----------------------------------------------------------------
    % (2) 分段函数值
    kT = piecewise_kT(T_g);    % 燃气绝热指数
    RT = piecewise_RT(T_g);    % 燃气气体常数
    d  = delta_cooling(T_g);   % 冷却引气系数

    %% ----------------------------------------------------------------
    % (3) 进口总压比 tau_v 与压气机入口温度 T_B
    %   tau_v = (1 + V^2 / (2*(k/(k-1))*R*T_H))^(k/(k-1))
    %   此处按图片最接近的数学结构实现
    inner = 1 + V_flight^2 / (2 * (k_air/(k_air-1)) * R_air * T_H);
    if inner <= 0, return; end
    exp1  = k_air / (k_air - 1);   % = 3.5 for air
    tau_v = inner^exp1;

    % 压气机入口温度（滞止温度）
    T_B = T_H * inner;   % = T_H * (1 + V^2/(2*(k/(k-1))*R*T_H))
    if ~isfinite(T_B) || T_B <= 0, return; end

    %% ----------------------------------------------------------------
    % (4) 压气机出口温度 T_k
    %   T_k = T_B * (1 + (pi_k^((k-1)/k) - 1) / eta_k)
    pi_k_ratio = pi_k^((k_air-1)/k_air);   % pi_k^(0.2857...)
    if ~isfinite(pi_k_ratio) || pi_k_ratio < 1, return; end
    T_k = T_B * (1 + (pi_k_ratio - 1) / eta_k);
    if ~isfinite(T_k) || T_k <= 0, return; end

    %% ----------------------------------------------------------------
    % (5) 相对耗油量 g_T
    %   g_T = 3.293e-5 * T_g - 2.84e-5 * T_k - 0.004814
    g_T = 3.293e-5 * T_g - 2.84e-5 * T_k - 0.004814;
    if ~isfinite(g_T) || g_T <= 0
        % g_T 应为正值；若参数组合导致 g_T<=0，视为无效
        return;
    end

    %% ----------------------------------------------------------------
    % (6) 热恢复系数 lambda_heat（即图中 λ(π_к', T_г')）
    %
    % 严格按图片公式实现（第4张图）：
    %
    %              1 - W_c / (H_g · η_к)
    %   λ = ─────────────────────────────────────────
    %              1 - W_c / (H_g · η_к · η_т · η_г)
    %
    % 其中：
    %   W_c = k/(k-1) · R · T_В · (π_к^((k-1)/k) - 1)   [压气机压缩功参数]
    %   H_g = k_T/(k_T-1) · R_T · T_г                    [燃气焓参数]
    %   η_к  → eta_k  （压气机效率）
    %   η_т  → eta_t  （涡轮效率）
    %   η_г  → eta_T  （燃烧放热系数）
    %
    % 注意：分子分母仅有 η_к，分母有 η_к · η_т · η_г，两者不同。

    % W_c：压气机压缩功参数
    W_c = (k_air / (k_air - 1)) * R_air * T_B * (pi_k_ratio - 1);

    % H_g：燃气焓参数  kT/(kT-1) · R_T · T_g
    H_g = (kT / (kT - 1)) * RT * T_g;
    if abs(H_g) < 1e-6, return; end

    % λ 分子：1 - W_c / (H_g · η_к)
    %   ——分子分母只含 eta_k，不含 eta_t、eta_T
    num_lambda = 1 - W_c / (H_g * eta_k);

    % λ 分母：1 - W_c / (H_g · η_к · η_т · η_г)
    %   ——分母分母含 eta_k · eta_t · eta_T（三个效率都在此处）
    den_lambda = 1 - W_c / (H_g * eta_k * eta_t * eta_T);

    if abs(den_lambda) < 1e-10, return; end
    lambda_heat = num_lambda / den_lambda;
    if ~isfinite(lambda_heat), return; end

    %% ----------------------------------------------------------------
    % (7) 进口总压恢复系数 σ_вх
    %   σ_вх = σ_сс · σ_кан
    sigma_bx = sigma_cc * sigma_kan;

    %% ----------------------------------------------------------------
    % (8) 单位自由能 L_св（即图中 L_св(π_к', T_г')）
    %
    % 严格按图片公式（第5张图）实现：
    %
    %   L_св = λ · [k_T/(k_T-1) · R_T · T_г · (1 - (1/Π)^((k_T-1)/k_T))]
    %        - k/(k-1) · R · T_В · (π_к^((k-1)/k) - 1)
    %          · 1 / [(1 + g_T) · η_к · η_т · η_М · (1 - δ)]
    %
    % 其中：
    %   Π = τ_v · σ_вх · π_к · σ_каск · σ_кс   [膨胀端总压比乘积]
    %   第一项：燃气膨胀可用功（由 λ 修正）
    %   第二项：驱动压气机消耗的功（除以效率链）
    %   第二项分母含 η_к · η_т · η_М（注意：此处是 η_т 涡轮效率，不是 η_г）

    % 膨胀端总压比乘积 Π = τ_v · σ_вх · π_к · σ_каск · σ_кс
    Pi_exp = tau_v * sigma_bx * pi_k * sigma_kask * sigma_ks;
    if Pi_exp <= 0, return; end

    % 膨胀指数 (k_T - 1)/k_T
    exp_kT = (kT - 1) / kT;

    % 膨胀项 (1/Π)^((k_T-1)/k_T)
    expansion_term = (1.0 / Pi_exp)^exp_kT;
    if ~isfinite(expansion_term), return; end

    % 第一项：λ · H_g · (1 - expansion_term)
    %   = λ · k_T/(k_T-1) · R_T · T_г · (1 - (1/Π)^((k_T-1)/k_T))
    Lsv_term1 = lambda_heat * H_g * (1 - expansion_term);

    % 第二项分母：(1 + g_T) · η_к · η_т · η_М · (1 - δ)
    %   注：此处为 eta_t（涡轮效率），与 lambda 分母的 eta_t 相同
    Lsv_denom2 = (1 + g_T) * eta_k * eta_t * eta_m * (1 - d);
    if abs(Lsv_denom2) < 1e-10, return; end

    % 第二项：W_c / Lsv_denom2
    %   = k/(k-1) · R · T_В · (π_к^((k-1)/k)-1) / [(1+g_T)·η_к·η_т·η_М·(1-δ)]
    Lsv_term2 = W_c / Lsv_denom2;

    % L_св = 第一项 - 第二项
    L_sv = Lsv_term1 - Lsv_term2;
    if ~isfinite(L_sv) || L_sv <= 0
        % L_sv 须为正（有效自由能为正才能产生推力）
        return;
    end

    %% ----------------------------------------------------------------
    % (9) 最优自由能分配系数 x_pc（即图中 x_pc(pi_k, T_g, m)）
    %
    %   x_pc = (1 + m*V^2 / (2*L_sv*eta_tv*eta_v*eta_c2))
    %        / (1 + m*eta_tv*eta_v*eta_c2 / (eta_c1*lambda_star))
    %
    %   此处按图片最接近的数学结构实现

    % 分子
    if L_sv <= 0, return; end
    V2_term = m * V_flight^2;   % = 0 when V_flight=0
    num_xpc = 1 + V2_term / (2 * L_sv * eta_tv * eta_v * eta_c2);

    % 分母
    den_xpc = 1 + (m * eta_tv * eta_v * eta_c2) / (eta_c1 * lambda_star);
    if abs(den_xpc) < 1e-10, return; end

    x_pc = num_xpc / den_xpc;
    if ~isfinite(x_pc) || x_pc <= 0 || x_pc >= 1
        % x_pc 应在 (0,1) 之间（自由能分配系数）
        % 物理要求宽松，允许稍大于 1 的值，但不允许负值
        if x_pc <= 0, return; end
    end

    %% ----------------------------------------------------------------
    % (10) 最优比推力 R_ud
    %
    %   R_ud = 1/(1+m) * [(1+g_T)*sqrt(2*eta_c1*lambda_star*x_pc*L_sv) - V_flight]
    %        + m/(1+m) * [sqrt(2*(1-x_pc)/m*L_sv*eta_tv*eta_v*eta_c2 + V_flight^2) - V_flight]
    %
    %   此处按图片最接近的数学结构实现

    % 一次喷管（内涵）速度分量
    inner_sq1 = 2 * eta_c1 * lambda_star * x_pc * L_sv;
    if inner_sq1 < 0, return; end
    V_j1 = (1 + g_T) * sqrt(inner_sq1) - V_flight;

    % 二次喷管（外涵）速度分量
    inner_sq2 = 2 * (1 - x_pc) / m * L_sv * eta_tv * eta_v * eta_c2 + V_flight^2;
    if inner_sq2 < 0, return; end
    V_j2 = sqrt(inner_sq2) - V_flight;

    R_ud = (1/(1+m)) * V_j1 + (m/(1+m)) * V_j2;
    if ~isfinite(R_ud) || R_ud <= 0, return; end

    %% ----------------------------------------------------------------
    % (11) 比油耗 C_ud
    %   C_ud = 3600 * g_T * (1 - delta) / (R_ud * (1 + m))
    %   单位：kg/(N·h)
    denom_C = R_ud * (1 + m);
    if abs(denom_C) < 1e-10, return; end
    C_ud = 3600 * g_T * (1 - d) / denom_C;
    if ~isfinite(C_ud) || C_ud <= 0, return; end

    %% ----------------------------------------------------------------
    % 输出
    y = [R_ud, C_ud];

    aux.T_B         = T_B;
    aux.T_k         = T_k;
    aux.tau_v       = tau_v;
    aux.g_T         = g_T;
    aux.lambda_heat = lambda_heat;
    aux.sigma_bx    = sigma_bx;
    aux.L_sv        = L_sv;
    aux.x_pc        = x_pc;
    aux.kT          = kT;
    aux.RT          = RT;
    aux.delta       = d;

catch ME
    % 捕获任何数值异常，返回 NaN
    % （调试时可取消注释下面这行）
    % warning('engine_forward caught: %s', ME.message);
end
end

%% -------------------------------------------------------------------------
%  生成虚拟试验数据
%% -------------------------------------------------------------------------
function data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed)
% generate_virtual_data  基于真值参数生成含噪声的虚拟观测
%
%   R_obs = R_true + N(0, sigma_R^2)，  sigma_R = noise_level_R * |R_true|
%   C_obs = C_true + N(0, sigma_C^2)，  sigma_C = noise_level_C * |C_true|

rng(rng_seed);  % 固定噪声种子

[y_true, ~] = engine_forward(theta_true, cond);
if ~all(isfinite(y_true))
    error('generate_virtual_data: 真值前向模型返回非有限值！');
end
R_true = y_true(1);
C_true = y_true(2);

sigma_R = noise_level_R * abs(R_true);
sigma_C = noise_level_C * abs(C_true);

R_obs = R_true + sigma_R * randn();
C_obs = C_true + sigma_C * randn();

data.R_true  = R_true;
data.C_true  = C_true;
data.R_obs   = R_obs;
data.C_obs   = C_obs;
data.sigma_R = sigma_R;
data.sigma_C = sigma_C;
end

%% -------------------------------------------------------------------------
%  对数先验
%% -------------------------------------------------------------------------
function lp = log_prior(theta, lb, ub)
% log_prior  均匀先验对数值
%   各参数独立均匀分布于 [lb_i, ub_i]
%   若 theta 落在范围内，返回 0（对数尺度下的常数）
%   否则返回 -Inf

if all(theta >= lb) && all(theta <= ub)
    lp = 0.0;
else
    lp = -Inf;
end
end

%% -------------------------------------------------------------------------
%  对数似然
%% -------------------------------------------------------------------------
function ll = log_likelihood(theta, data, cond)
% log_likelihood  高斯观测噪声下的对数似然
%
%   logL = -0.5*[((R_obs-R_pred)/sigma_R)^2 + log(2*pi*sigma_R^2)
%               +((C_obs-C_pred)/sigma_C)^2 + log(2*pi*sigma_C^2)]

ll = -Inf;  % 默认无效

[y_pred, ~] = engine_forward(theta, cond);
if ~all(isfinite(y_pred))
    return;
end
R_pred = y_pred(1);
C_pred = y_pred(2);

R_obs   = data.R_obs;
C_obs   = data.C_obs;
sigma_R = data.sigma_R;
sigma_C = data.sigma_C;

if sigma_R <= 0 || sigma_C <= 0
    return;
end

% 对数似然（高斯噪声，两观测量独立）
res_R = (R_obs - R_pred) / sigma_R;
res_C = (C_obs - C_pred) / sigma_C;

ll = -0.5 * (res_R^2 + log(2*pi*sigma_R^2) ...
           + res_C^2 + log(2*pi*sigma_C^2));

if ~isfinite(ll)
    ll = -Inf;
end
end

%% -------------------------------------------------------------------------
%  对数后验
%% -------------------------------------------------------------------------
function lpost = log_posterior_fn(theta, data, cond, lb, ub)
% log_posterior_fn  对数后验 = 对数先验 + 对数似然
lp = log_prior(theta, lb, ub);
if isinf(lp)
    lpost = -Inf;
    return;
end
ll    = log_likelihood(theta, data, cond);
lpost = lp + ll;
end

%% -------------------------------------------------------------------------
%  反射边界处理
%% -------------------------------------------------------------------------
function z_ref = reflect_to_unit_interval(z)
% reflect_to_unit_interval  将归一化参数反射回 [0,1]
%   对越界分量使用镜像反射，而不是简单截断
%   这样可以避免采样在边界处堆积

z_ref = z;
for i = 1:length(z)
    zi = z(i);
    % 反复折叠，直到落入 [0,1]
    max_iter = 10;
    iter = 0;
    while (zi < 0 || zi > 1) && iter < max_iter
        if zi < 0
            zi = -zi;           % 关于 0 反射
        elseif zi > 1
            zi = 2 - zi;        % 关于 1 反射
        end
        iter = iter + 1;
    end
    % 若多次反射后仍越界（极端情况），强制截断
    z_ref(i) = max(0, min(1, zi));
end
end

%% -------------------------------------------------------------------------
%  MCMC 采样器（随机游走 Metropolis-Hastings）
%% -------------------------------------------------------------------------
function results = run_mcmc(data, cond, lb, ub, theta0, opts)
% run_mcmc  归一化空间中的随机游走 MH 采样
%
%   在归一化空间 z = (theta-lb)/(ub-lb) ∈ [0,1]^n 中进行提案，
%   越界分量用反射边界处理，然后映射回物理空间计算对数后验。

n_params   = length(lb);
n_samples  = opts.n_samples;
burn_in    = opts.burn_in;
prop_sd    = opts.proposal_sd;

% 链存储（全部步数）
chain_full   = zeros(n_samples, n_params);
logpost_full = zeros(n_samples, 1);

% 初始化
theta_curr = theta0;
% 确保初始值在范围内
theta_curr = max(lb, min(ub, theta_curr));
z_curr = (theta_curr - lb) ./ (ub - lb);

lpost_curr = log_posterior_fn(theta_curr, data, cond, lb, ub);
if ~isfinite(lpost_curr)
    % 如果初始点无效，尝试先验中点
    warning('MCMC: 初始点对数后验无效，改用先验中点...');
    theta_curr = 0.5*(lb+ub);
    z_curr     = 0.5 * ones(1, n_params);
    lpost_curr = log_posterior_fn(theta_curr, data, cond, lb, ub);
    if ~isfinite(lpost_curr)
        error('MCMC: 先验中点也无效，请检查前向模型与参数设置！');
    end
end

n_accept_total = 0;
n_accept_window = 0;  % 自适应窗口内的接受次数

fprintf('  进度: ');
for s = 1:n_samples
    % 打印进度
    if mod(s, opts.print_interval) == 0
        fprintf('%d/%d (accept=%.2f, prop_sd=%.4f)\n', ...
            s, n_samples, n_accept_total/s, prop_sd);
    end

    % 提案（归一化空间中高斯游走）
    z_prop = z_curr + prop_sd * randn(1, n_params);

    % 反射边界
    z_prop = reflect_to_unit_interval(z_prop);

    % 映射回物理空间
    theta_prop = lb + z_prop .* (ub - lb);

    % 计算对数后验
    lpost_prop = log_posterior_fn(theta_prop, data, cond, lb, ub);

    % MH 接受准则
    log_alpha = lpost_prop - lpost_curr;
    if log(rand()) < log_alpha
        % 接受
        z_curr     = z_prop;
        theta_curr = theta_prop;
        lpost_curr = lpost_prop;
        n_accept_total  = n_accept_total + 1;
        n_accept_window = n_accept_window + 1;
    end

    % 存储
    chain_full(s, :) = theta_curr;
    logpost_full(s)  = lpost_curr;

    % 自适应调整提案方差（在 adapt_start 到 adapt_end 之间）
    if s >= opts.adapt_start && s <= opts.adapt_end
        if mod(s - opts.adapt_start, opts.adapt_interval) == 0 && s > opts.adapt_start
            local_rate = n_accept_window / opts.adapt_interval;
            if local_rate < opts.target_accept_low
                prop_sd = prop_sd * 0.9;
            elseif local_rate > opts.target_accept_high
                prop_sd = prop_sd * 1.1;
            end
            % 限制 prop_sd 范围，防止过大或过小
            prop_sd = max(1e-4, min(prop_sd, 0.5));
            n_accept_window = 0;  % 重置窗口计数
        end
    end
end
fprintf('\n');

% 丢弃 burn-in
chain_post   = chain_full(burn_in+1:end, :);
logpost_post = logpost_full(burn_in+1:end);

% 后验统计
theta_mean = mean(chain_post, 1);
theta_std  = std(chain_post, 0, 1);

% MAP（后验链中对数后验最大的点）
[~, best_idx_post] = max(logpost_post);
theta_map = chain_post(best_idx_post, :);

% 95% 后验区间
theta_ci95 = zeros(n_params, 2);
for i = 1:n_params
    theta_ci95(i, :) = quantile(chain_post(:,i), [0.025, 0.975]);
end

% 输出结构体
results.chain_full   = chain_full;
results.chain_post   = chain_post;
results.logpost_full = logpost_full;
results.logpost_post = logpost_post;
results.accept_rate  = n_accept_total / n_samples;
results.theta_mean   = theta_mean;
results.theta_std    = theta_std;
results.theta_map    = theta_map;
results.theta_ci95   = theta_ci95;
results.best_idx     = best_idx_post;
results.prop_sd_final= prop_sd;
end

%% -------------------------------------------------------------------------
%  绘图函数
%% -------------------------------------------------------------------------
function plot_results(results, theta_true, lb, ub, param_names)
% plot_results  生成后验分析图形（共3张：链轨迹、边缘后验、相关热图）

chain_post = results.chain_post;
chain_full = results.chain_full;
n_params   = size(chain_post, 2);
n_post     = size(chain_post, 1);
n_full     = size(chain_full, 1);

%% --- 图1：链轨迹图（全链，含 burn-in）---
fig1 = figure('Name','链轨迹图','Position',[50 50 1400 900]);
n_cols = 3;
n_rows = ceil(n_params / n_cols);
for i = 1:n_params
    subplot(n_rows, n_cols, i);
    plot(1:n_full, chain_full(:,i), 'b-', 'LineWidth', 0.3);
    hold on;
    xline(size(chain_full,1) - n_post, 'r--', 'LineWidth', 1.5, ...
        'Label','burn-in结束');
    yline(theta_true(i), 'g-', 'LineWidth', 1.5);
    xlabel('迭代步数');
    ylabel(param_names{i});
    title(sprintf('%s 链轨迹', param_names{i}));
    legend('链','burn-in边界','真值','Location','best','FontSize',6);
    ylim([lb(i)-0.005*(ub(i)-lb(i)), ub(i)+0.005*(ub(i)-lb(i))]);
end
sgtitle('MCMC 链轨迹图（蓝线=全链，绿线=真值，红虚线=burn-in结束）');

%% --- 图2：边缘后验直方图 ---
fig2 = figure('Name','边缘后验直方图','Position',[100 50 1400 900]);
for i = 1:n_params
    subplot(n_rows, n_cols, i);
    histogram(chain_post(:,i), 40, 'Normalization', 'pdf', ...
        'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    xline(theta_true(i),          'g-',  'LineWidth', 2);
    xline(results.theta_mean(i),  'r--', 'LineWidth', 2);
    xline(results.theta_map(i),   'm:',  'LineWidth', 2);
    xlabel(param_names{i});
    ylabel('概率密度');
    title(sprintf('%s\n真值=%.4f, 均值=%.4f, MAP=%.4f', ...
        param_names{i}, theta_true(i), results.theta_mean(i), results.theta_map(i)));
    legend('后验','真值','后验均值','MAP','Location','best','FontSize',6);
    xlim([lb(i), ub(i)]);
end
sgtitle('边缘后验分布直方图');

%% --- 图3：参数相关系数热图 ---
fig3 = figure('Name','后验相关性热图','Position',[150 50 700 600]);
C_mat = corrcoef(chain_post);
imagesc(C_mat, [-1, 1]);
colorbar;
colormap(redblue_colormap());
set(gca, 'XTick', 1:n_params, 'XTickLabel', param_names, ...
         'XTickLabelRotation', 45, ...
         'YTick', 1:n_params, 'YTickLabel', param_names, ...
         'FontSize', 8);
title('后验参数相关系数热图（-1=负相关，+1=正相关）');
% 在格子中显示数值
for i = 1:n_params
    for j = 1:n_params
        text(j, i, sprintf('%.2f', C_mat(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 6, ...
            'Color', txt_color(abs(C_mat(i,j))));
    end
end

%% --- 保存所有图形 ---
try
    saveas(fig1, 'fig1_trace_plots.png');
    saveas(fig2, 'fig2_marginal_posteriors.png');
    saveas(fig3, 'fig3_correlation_heatmap.png');
    fprintf('图形已保存至当前目录（fig1~fig3）。\n');
catch
    fprintf('图形保存失败（可能是无显示环境），但已正常绘制。\n');
end
end

%% -------------------------------------------------------------------------
%  辅助：文字颜色选择（根据背景深浅）
%% -------------------------------------------------------------------------
function c = txt_color(abs_val)
% txt_color  若绝对值 > 0.5 则用白色，否则用黑色（用于热图数字标注）
if abs_val > 0.5
    c = 'white';
else
    c = 'black';
end
end

%% -------------------------------------------------------------------------
%  辅助：生成红蓝配色图
%% -------------------------------------------------------------------------
function cmap = redblue_colormap()
% redblue_colormap  生成 64 级蓝-白-红色图，用于相关矩阵热图
n = 32;
% 前半段：蓝 (0,0,1) -> 白 (1,1,1)
half1 = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)];
% 后半段：白 (1,1,1) -> 红 (1,0,0)
half2 = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];
cmap  = [half1; half2];
end
