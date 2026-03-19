# 涡扇发动机性能参数 MCMC 贝叶斯反演程序说明

**文件**：`turbine_mcmc_inversion.m`
**语言**：MATLAB
**方法**：Metropolis-Hastings 马尔可夫链蒙特卡洛（MCMC）贝叶斯反演

---

## 目录

1. [程序概述](#1-程序概述)
2. [问题背景与物理模型](#2-问题背景与物理模型)
3. [参数体系](#3-参数体系)
4. [程序结构与执行流程](#4-程序结构与执行流程)
5. [子函数详解](#5-子函数详解)
6. [贝叶斯推断框架](#6-贝叶斯推断框架)
7. [MCMC 算法说明](#7-mcmc-算法说明)
8. [输出结果说明](#8-输出结果说明)
9. [参数配置指南](#9-参数配置指南)
10. [运行方法](#10-运行方法)

---

## 1. 程序概述

本程序针对**涡扇发动机**，以比推力 $R_{ud}$（N·s/kg）和比油耗 $C_{ud}$（kg/(N·h)）为可观测量，对4个核心效率参数进行**贝叶斯统计反演**。

**核心思路：**

$$
p(\boldsymbol{\theta} \mid \mathbf{y}_\text{obs}) \propto p(\mathbf{y}_\text{obs} \mid \boldsymbol{\theta}) \cdot p(\boldsymbol{\theta})
$$

即：**后验分布 ∝ 似然函数 × 先验分布**

通过 MCMC 采样近似后验分布，从而以概率形式估计未知参数的取值及其不确定性。

---

## 2. 问题背景与物理模型

### 2.1 发动机类型

大涵道比涡扇发动机，设计工况为**地面静止状态**（飞行马赫数 $M=0$）。

### 2.2 工况参数

| 参数 | 符号 | 值 | 单位 |
|------|------|----|------|
| 大气总温 | $T_H$ | 288 | K |
| 飞行马赫数 | $M_\text{flight}$ | 0.0 | — |
| 涵道比 | $m$ | 10.0 | — |
| 压气机压比 | $\pi_k$ | 33.0 | — |
| 涡轮前总温 | $T_g$ | 1700 | K |

### 2.3 前向模型计算链

前向模型 `engine_forward` 按以下顺序计算：

```
输入参数 θ、工况 cond
     │
     ▼
(1) 声速 a = √(k_air · R_air · T_H)，飞行速度 V_flight = a · M
     │
     ▼
(2) 分段查表：燃气绝热指数 kT，燃气气体常数 RT，冷却引气系数 δ
     │
     ▼
(3) 进口总压比 τ_v，压气机入口温度 T_B
     │
     ▼
(4) 压气机出口温度 T_k = T_B · [1 + (π_k^((k-1)/k) - 1) / η_k]
     │
     ▼
(5) 相对耗油量 g_T = 3×10⁻⁵ · T_g - 2.69×10⁻⁵ · T_k - 0.003
     │
     ▼
(6) 热恢复系数 λ_heat（分子/分母形式）
     │
     ▼
(7) 综合总压恢复系数 σ_bx = σ_cc · σ_kan
     │
     ▼
(8) 单位自由能 L_sv = λ_heat · (term1 - term2)
     │
     ▼
(9) 最优自由能分配系数 x_pc
     │
     ▼
(10) 主流喷射速度 V_j1，风扇喷射速度 V_j2 → 比推力 R_ud
     │
     ▼
(11) 比油耗 C_ud = 3600 · g_T · (1 - δ) / [R_ud · (1 + m)]
     │
     ▼
输出：y = [R_ud, C_ud]
```

**各步关键公式：**

| 步骤 | 公式 |
|------|------|
| 压气机出口温度 | $T_k = T_B\left(1 + \dfrac{\pi_k^{(k-1)/k}-1}{\eta_k}\right)$ |
| 热恢复系数 | $\lambda_\text{heat} = \dfrac{1 - W_k/(h_g \eta_k)}{1 - W_k/(h_g \eta_k \eta_t)}$ |
| 单位自由能 | $L_{sv} = \lambda_\text{heat}(L_1 - L_2)$ |
| 比推力 | $R_{ud} = \dfrac{1}{1+m}V_{j1} + \dfrac{m}{1+m}V_{j2}$ |
| 比油耗 | $C_{ud} = \dfrac{3600\,g_T(1-\delta)}{R_{ud}(1+m)}$ |

---

## 3. 参数体系

### 3.1 待反演参数（4个敏感参数）

| 序号 | 符号 | 物理含义 | 真值 | 先验下界 | 先验上界 | 初始值 |
|------|------|----------|------|----------|----------|--------|
| 1 | `eta_k` ($\eta_k$) | 压气机绝热效率 | 0.850 | 0.80 | 0.90 | 0.850 |
| 2 | `eta_t` ($\eta_t$) | 涡轮绝热效率 | 0.890 | 0.83 | 0.95 | 0.890 |
| 3 | `eta_T` ($\eta_T$) | 燃烧放热系数 | 0.980 | 0.95 | 1.01 | 0.980 |
| 4 | `lambda` ($\lambda$) | 风扇涡轮热恢复系数 | 1.030 | 1.00 | 1.06 | 1.030 |

> 初始值为先验区间中点 `theta0 = 0.5 * (lb + ub)`

### 3.2 固定参数（9个，不参与反演）

| 序号 | 符号 | 物理含义 | 固定值 |
|------|------|----------|--------|
| 1 | `eta_m` | 机械效率 | 0.988 |
| 2 | `eta_v` | 风扇效率 | 0.860 |
| 3 | `eta_tv` | 风扇涡轮效率 | 0.910 |
| 4 | `eta_c1` | 一次喷管效率 | 0.945 |
| 5 | `eta_c2` | 二次喷管效率 | 0.930 |
| 6 | `sigma_cc` | 进气道激波总压恢复 | 0.990 |
| 7 | `sigma_kan` | 进气通道总压恢复 | 0.985 |
| 8 | `sigma_kask` | 压气机级间总压恢复 | 0.985 |
| 9 | `sigma_ks` | 燃烧室总压恢复 | 0.950 |

### 3.3 分段物性参数

| 涡轮前温度 $T_g$ 范围 | 燃气绝热指数 $k_T$ | 燃气气体常数 $R_T$ (J/kg·K) |
|----------------------|-------------------|------------------------------|
| $800 < T_g \leq 1400$ K | 1.33 | 287.6 |
| $1400 < T_g \leq 1600$ K | 1.30 | 288.0 |
| $T_g > 1600$ K | 1.25 | 288.6 |

**冷却引气系数：**

$$
\delta = \text{clip}\left(0.02 + \frac{T_g - 1200}{100} \times 0.02,\ 0,\ 0.15\right)
$$

当 $T_g = 1700$ K 时：$\delta = 0.02 + 5 \times 0.02 = 0.12$

---

## 4. 程序结构与执行流程

```
主程序 turbine_mcmc_inversion.m
├── §1  参数定义
│     ├── 敏感参数名、先验上下界
│     ├── 工况参数 cond
│     ├── 真值 theta_true，固定参数 theta_fixed
│     └── 初始值 theta0（先验中点）
│
├── §2  前向模型验证
│     └── engine_forward(theta_true) → 打印 R_ud, C_ud
│
├── §3  生成虚拟观测数据
│     └── generate_virtual_data → 真值 + 1% 高斯噪声
│
├── §4  MCMC 配置
│     └── 步数、烧入期、提议标准差、自适应参数
│
├── §5  运行 MCMC
│     └── run_mcmc → results 结构体
│
├── §6  后验统计打印
│     └── 均值、MAP、95%置信区间
│
├── §7  后验预测对比
│     └── 用后验均值/MAP前向计算，比较与真值的相对误差
│
└── §8  绘图
      └── plot_results → 链轨迹图 + 边缘后验直方图 → 保存PNG
```

---

## 5. 子函数详解

### 5.1 `piecewise_kT(T_g)` — 燃气绝热指数

```matlab
function kT = piecewise_kT(T_g)
```

根据涡轮前温度 $T_g$ 分段返回燃气绝热指数 $k_T$，反映高温燃气热容比随温度变化的规律。

---

### 5.2 `piecewise_RT(T_g)` — 燃气气体常数

```matlab
function RT = piecewise_RT(T_g)
```

分段返回燃气气体常数 $R_T$（J/kg·K），高温时燃气组分变化导致气体常数略有不同。

---

### 5.3 `delta_cooling(T_g)` — 冷却引气系数

```matlab
function d = delta_cooling(T_g)
```

线性插值计算涡轮冷却引气占总气流的比例，范围限制在 $[0, 0.15]$。温度越高，需要引更多冷却气。

---

### 5.4 `engine_forward(theta_sensitive, theta_fixed, cond)` — 前向模型

```matlab
function [y, aux] = engine_forward(theta_sensitive, theta_fixed, cond)
```

**输入：**
- `theta_sensitive`：4个敏感参数向量 $[\eta_k, \eta_t, \eta_T, \lambda]$
- `theta_fixed`：9个固定参数向量
- `cond`：工况结构体

**输出：**
- `y = [R_ud, C_ud]`：比推力和比油耗
- `aux`：中间变量结构体（$T_B, T_k, g_T, \lambda_\text{heat}, L_{sv}, x_{pc}$ 等）

**数值保护：** 全程对每个中间量进行有限性检验（`isfinite`），出现 NaN/Inf 或物理不合理值（≤0）时立即返回 `[NaN, NaN]`，保证 MCMC 采样不因数值异常崩溃。

---

### 5.5 `generate_virtual_data(...)` — 生成虚拟观测数据

```matlab
function data = generate_virtual_data(theta_true, theta_fixed, cond, noise_level_R, noise_level_C)
```

用真值参数调用前向模型，得到"真实"性能值，再叠加高斯白噪声模拟观测误差：

$$
R_\text{obs} = R_\text{true} + \sigma_R \cdot \varepsilon_R, \quad \varepsilon_R \sim \mathcal{N}(0,1)
$$

$$
\sigma_R = \text{noise\_level\_R} \times |R_\text{true}|, \quad \sigma_C = \text{noise\_level\_C} \times |C_\text{true}|
$$

当前噪声水平：**1%**（`noise_level_R = noise_level_C = 0.01`）

---

### 5.6 `log_prior(theta, lb, ub)` — 对数先验

```matlab
function lp = log_prior(theta, lb, ub)
```

采用**均匀先验**（Uniform Prior）：

$$
p(\boldsymbol{\theta}) = \begin{cases} 1 & \text{若 } \boldsymbol{\theta} \in [\mathbf{lb},\, \mathbf{ub}] \\ 0 & \text{其他} \end{cases}
$$

对数形式：参数在边界内返回 0，否则返回 $-\infty$。

---

### 5.7 `log_likelihood(theta, theta_fixed, data, cond)` — 对数似然

```matlab
function ll = log_likelihood(theta, theta_fixed, data, cond)
```

假设观测误差服从独立高斯分布，对数似然为：

$$
\ln p(\mathbf{y}_\text{obs} \mid \boldsymbol{\theta}) = -\frac{1}{2}\left[\left(\frac{R_\text{obs}-R_\text{pred}}{\sigma_R}\right)^2 + \ln(2\pi\sigma_R^2) + \left(\frac{C_\text{obs}-C_\text{pred}}{\sigma_C}\right)^2 + \ln(2\pi\sigma_C^2)\right]
$$

---

### 5.8 `log_posterior_fn(...)` — 对数后验

```matlab
function lpost = log_posterior_fn(theta, theta_fixed, data, cond, lb, ub)
```

贝叶斯定理的对数形式：

$$
\ln p(\boldsymbol{\theta} \mid \mathbf{y}_\text{obs}) = \ln p(\boldsymbol{\theta}) + \ln p(\mathbf{y}_\text{obs} \mid \boldsymbol{\theta}) + \text{const}
$$

若先验为 $-\infty$（越界），直接短路返回，不再调用前向模型（提升效率）。

---

### 5.9 `run_mcmc(...)` — MCMC 主采样函数

```matlab
function results = run_mcmc(data, theta_fixed, cond, lb, ub, theta0, opts)
```

详见 [§7 MCMC 算法说明](#7-mcmc-算法说明)。

**返回的 `results` 结构体字段：**

| 字段 | 含义 |
|------|------|
| `chain_full` | 完整 MCMC 链（含烧入期），形状 `[n_samples, 4]` |
| `chain_post` | 后验链（去除烧入期），形状 `[n_samples-burn_in, 4]` |
| `logpost_full` | 每步对数后验值 |
| `logpost_post` | 后验段对数后验值 |
| `accept_rate` | 总体接受率 |
| `theta_mean` | 后验均值（1×4） |
| `theta_std` | 后验标准差（1×4） |
| `theta_map` | MAP 估计（后验概率最大点，1×4） |
| `theta_ci95` | 95% 可信区间（4×2，每行 [2.5%, 97.5%]） |
| `prop_sd_final` | 最终自适应步长 |

---

### 5.10 `plot_results(...)` — 可视化函数

```matlab
function plot_results(results, theta_true, lb, ub, param_names)
```

生成两张图并保存为 PNG：

**图1：MCMC 链轨迹图** (`链轨迹图_4参数.png`)
- 蓝线：参数链的历史轨迹
- 红虚线：烧入期结束边界
- 绿线：参数真值参考线

**图2：边缘后验分布直方图** (`边缘后验图_4参数.png`)
- 蓝色直方图：后验概率密度（PDF 归一化）
- 绿线：真值
- 红虚线：后验均值
- 紫点线：MAP 估计

---

## 6. 贝叶斯推断框架

### 6.1 模型设定

$$
\mathbf{y}_\text{obs} = \mathbf{f}(\boldsymbol{\theta}) + \boldsymbol{\varepsilon}, \quad \boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma})
$$

其中：
- $\mathbf{f}(\boldsymbol{\theta}) = [R_{ud}(\boldsymbol{\theta}),\ C_{ud}(\boldsymbol{\theta})]$：前向模型
- $\boldsymbol{\Sigma} = \text{diag}(\sigma_R^2,\ \sigma_C^2)$：观测噪声协方差（已知）

### 6.2 先验

均匀先验（无信息先验），等价于在参数可行域内纯由数据驱动：

$$
\boldsymbol{\theta} \sim \mathcal{U}(\mathbf{lb},\ \mathbf{ub})
$$

### 6.3 后验估计指标

| 指标 | 定义 | 含义 |
|------|------|------|
| 后验均值 | $\bar{\boldsymbol{\theta}} = \mathbb{E}[\boldsymbol{\theta} \mid \mathbf{y}]$ | MMSE 估计，最优期望意义 |
| MAP | $\hat{\boldsymbol{\theta}} = \arg\max_{\boldsymbol{\theta}} p(\boldsymbol{\theta}\mid\mathbf{y})$ | 后验概率最大点 |
| 95% 可信区间 | $[Q_{2.5\%},\ Q_{97.5\%}]$ | 参数真值有95%概率落在此区间内 |

---

## 7. MCMC 算法说明

### 7.1 算法：Metropolis-Hastings（MH）

**提议分布：** 在标准化参数空间（$[0,1]^4$）中使用各向同性高斯随机游走：

$$
\mathbf{z}^\star = \mathbf{z}^{(s)} + \sigma_\text{prop} \cdot \boldsymbol{\xi}, \quad \boldsymbol{\xi} \sim \mathcal{N}(\mathbf{0}, \mathbf{I})
$$

$$
\boldsymbol{\theta}^\star = \mathbf{lb} + \mathbf{z}^\star \odot (\mathbf{ub} - \mathbf{lb})
$$

> 在标准化空间中提议，使各参数量纲统一，步长更易调控。

**接受准则：**

$$
\alpha = \min\left(1,\ \frac{p(\boldsymbol{\theta}^\star \mid \mathbf{y})}{p(\boldsymbol{\theta}^{(s)} \mid \mathbf{y})}\right)
$$

以概率 $\alpha$ 接受新样本，否则保留当前样本。

### 7.2 自适应步长

在前 `adapt_end = 5000` 步内，每 `adapt_interval = 200` 步检查局部接受率，自动调整提议标准差：

```
若接受率 < 20%  →  σ_prop × 0.9   （减小步长，提高接受率）
若接受率 > 35%  →  σ_prop × 1.1   （增大步长，提升探索效率）
σ_prop 限制在 [1e-5, 0.5] 内
```

**目标接受率：20%～35%**（高维 MH 的理论最优约为 23.4%）

### 7.3 烧入期（Burn-in）

前 10000 步为烧入期，链从初始值向高后验概率区域收敛，该阶段样本不用于统计。

### 7.4 MCMC 配置参数

| 参数 | 值 | 说明 |
|------|----|------|
| `n_samples` | 40000 | 总采样步数 |
| `burn_in` | 10000 | 烧入期步数（占25%） |
| 后验链长度 | 30000 | 用于统计的有效样本数 |
| `proposal_sd` | 0.0005 | 初始提议标准差（标准化空间） |
| `adapt_start` | 0 | 自适应起始步 |
| `adapt_end` | 5000 | 自适应结束步 |
| `adapt_interval` | 200 | 自适应调整间隔 |

---

## 8. 输出结果说明

### 8.1 控制台输出顺序

```
===== 工况参数 =====
===== 4个敏感参数信息 =====
===== 前向模型验证 =====
===== 观测数据（噪声1%）=====
===== MCMC 配置 =====
===== 运行 MCMC =====
  已完成 5000/40000 步，当前接受率: x.xxx
  ...
MCMC 完成！总体接受率: x.xxx
===== 后验结果统计 =====
===== 后验预测对比 =====
图片已保存：...
```

### 8.2 后验结果统计表

```
参数名           真值    后验均值    MAP    CI95低   CI95高
----------------------------------------------------------------------
eta_k         0.8500   0.xxxx   0.xxxx   0.xxxx   0.xxxx
eta_t         0.8900   0.xxxx   0.xxxx   0.xxxx   0.xxxx
eta_T         0.9800   0.xxxx   0.xxxx   0.xxxx   0.xxxx
lambda        1.0300   0.xxxx   0.xxxx   0.xxxx   0.xxxx
```

### 8.3 后验预测对比表

```
性能指标              R_ud [N·s/kg]   C_ud [kg/(N·h)]
真值                     xxx.xxxx        0.xxxxxx
观测值                   xxx.xxxx        0.xxxxxx
后验均值预测             xxx.xxxx        0.xxxxxx
MAP 预测                 xxx.xxxx        0.xxxxxx
```

### 8.4 输出图像

| 文件名 | 内容 |
|--------|------|
| `链轨迹图_4参数.png` | 4个参数的 MCMC 采样轨迹，含烧入期边界和真值参考线 |
| `边缘后验图_4参数.png` | 4个参数的边缘后验概率密度直方图，含真值/均值/MAP标注 |

---

## 9. 参数配置指南

### 9.1 噪声水平

```matlab
noise_level_R = 0.01;  % R_ud 观测噪声标准差 = 1% × 真值
noise_level_C = 0.01;  % C_ud 观测噪声标准差 = 1% × 真值
```

噪声越小，后验分布越窄，参数估计越精确。

### 9.2 MCMC 调参建议

| 问题现象 | 建议调整 |
|----------|----------|
| 接受率过低（< 10%） | 减小 `proposal_sd` |
| 接受率过高（> 50%） | 增大 `proposal_sd` |
| 链未收敛 | 增大 `n_samples`，增大 `burn_in` |
| 后验分布过宽 | 降低噪声水平，或收窄先验范围 |
| 结果可复现 | 取消 `rng(42)` 注释 |

### 9.3 先验范围调整

```matlab
lb = [0.80, 0.83, 0.95, 1.00];  % 下界
ub = [0.90, 0.95, 1.01, 1.06];  % 上界
```

缩窄先验范围可加速收敛，但需确保真值在范围内（程序有 `assert` 检验）。

---

## 10. 运行方法

### 10.1 环境要求

- MATLAB R2019b 或更高版本
- 无需额外工具箱（仅用 `quantile` 函数，Statistics Toolbox）

### 10.2 运行步骤

1. 将 `turbine_mcmc_inversion.m` 置于 MATLAB 工作目录
2. 在命令行窗口运行：
   ```matlab
   turbine_mcmc_inversion
   ```
3. 等待 MCMC 采样完成（约数十秒至数分钟，取决于机器性能）
4. 查看控制台输出和自动保存的两张 PNG 图片

### 10.3 修改真值或工况

修改以下变量后重新运行即可：

```matlab
% 修改工况
cond.T_g = 1600.0;   % 更改涡轮前温度

% 修改真值（验证反演能力）
theta_true = [0.83, 0.87, 0.975, 1.025];

% 修改固定参数
theta_fixed(1) = 0.985;  % 修改机械效率
```

---

*本文档由程序代码自动分析生成，2025年*
