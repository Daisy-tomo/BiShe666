# 涡扇发动机参数贝叶斯反演程序说明

**文件名：** `bayesian_inversion_single_condition.m`
**适用环境：** MATLAB R2019b 及以上
**运行时间：** 普通笔记本电脑约 1～3 分钟

---

## 目录

1. [研究背景与问题定义](#1-研究背景与问题定义)
2. [程序整体结构](#2-程序整体结构)
3. [13个未知修正参数](#3-13个未知修正参数)
4. [单工况设置](#4-单工况设置)
5. [前向热力循环模型](#5-前向热力循环模型)
6. [贝叶斯反演框架](#6-贝叶斯反演框架)
7. [MCMC采样器设计](#7-mcmc采样器设计)
8. [虚拟试验数据生成](#8-虚拟试验数据生成)
9. [输出结果与可视化](#9-输出结果与可视化)
10. [鲁棒性设计](#10-鲁棒性设计)
11. [局部函数索引](#11-局部函数索引)
12. [常见问题](#12-常见问题)

---

## 1. 研究背景与问题定义

### 1.1 问题描述

本程序针对**大涵道比双涵道涡扇发动机**，在**单个工作条件**下，利用两个可观测的宏观性能指标：

- **比推力** $R_{ud}$（单位：N·s/kg）—— 单位空气流量产生的推力
- **比油耗** $C_{ud}$（单位：kg/(N·h)）—— 产生单位推力每小时消耗的燃油量

对发动机热力循环中的 **13 个修正参数**进行**贝叶斯后验反演**。

### 1.2 问题的欠定性

| 量 | 数目 |
|---|---|
| 观测量（$R_{ud}$，$C_{ud}$） | 2 |
| 未知参数 | 13 |
| 欠定程度 | 严重欠定（13 − 2 = 11 个自由度） |

由于观测量远少于未知参数，**不存在唯一确定解**。程序的目标是：

- 给出每个参数的**后验分布**（而非单点估计）
- 量化参数间的**后验相关性**
- 分析每个参数的**可辨识性**（先验信息能被数据更新多少）

### 1.3 贝叶斯视角

$$
\underbrace{p(\boldsymbol{\theta} \mid \mathbf{y})}_{\text{后验}} \propto \underbrace{p(\mathbf{y} \mid \boldsymbol{\theta})}_{\text{似然}} \cdot \underbrace{p(\boldsymbol{\theta})}_{\text{先验}}
$$

其中：
- $\boldsymbol{\theta}$：13 维参数向量
- $\mathbf{y} = [R_{ud}^{obs},\; C_{ud}^{obs}]$：含噪声的观测值
- 先验取均匀分布（参数范围内等可能）
- 似然取高斯噪声模型

---

## 2. 程序整体结构

```
bayesian_inversion_single_condition.m
│
├─ 主程序（script 部分）
│   ├─ Part 1: 参数定义（param_names, lb, ub）
│   ├─ Part 2: 单工况设置（cond 结构体）
│   ├─ Part 3: 真值参数 theta_true + 初始值 theta0
│   ├─ Part 4: 前向模型验证
│   ├─ Part 5: 虚拟试验数据生成
│   ├─ Part 6: MCMC 参数设置（opts 结构体）
│   ├─ Part 7: 运行 MCMC 采样
│   ├─ Part 8: 打印后验统计汇总
│   ├─ Part 9: 后验预测与误差对比
│   └─ Part 10: 调用绘图函数
│
└─ 局部函数（functions 部分）
    ├─ piecewise_kT      分段燃气绝热指数
    ├─ piecewise_RT      分段燃气气体常数
    ├─ delta_cooling     冷却引气系数
    ├─ engine_forward    完整前向热力循环模型
    ├─ generate_virtual_data   虚拟数据生成
    ├─ log_prior         对数先验
    ├─ log_likelihood    对数似然
    ├─ log_posterior_fn  对数后验
    ├─ reflect_to_unit_interval  反射边界处理
    ├─ run_mcmc          MH-MCMC 采样器
    ├─ plot_results      5张图形绘制
    ├─ txt_color         热图文字颜色辅助
    └─ redblue_colormap  蓝-白-红色图生成
```

---

## 3. 13个未知修正参数

参数向量按以下顺序排列：

$$
\boldsymbol{\theta} = [\eta_k,\; \eta_t,\; \eta_m,\; \eta_v,\; \eta_{tv},\; \eta_{c1},\; \eta_{c2},\; \sigma_{cc},\; \sigma_{kan},\; \sigma_{kask},\; \sigma_{ks},\; \eta_T,\; \lambda^*]
$$

| 序号 | 参数名 | 物理含义 | 范围下界 | 范围上界 | 程序内真值 |
|:---:|:---:|---|:---:|:---:|:---:|
| 1 | $\eta_k$ | 压气机绝热效率 | 0.84 | 0.92 | 0.88 |
| 2 | $\eta_t$ | 涡轮绝热效率 | 0.86 | 0.92 | 0.89 |
| 3 | $\eta_m$ | 机械效率（传动损失） | 0.980 | 0.995 | 0.988 |
| 4 | $\eta_v$ | 风扇效率 | 0.85 | 0.87 | 0.860 |
| 5 | $\eta_{tv}$ | 风扇涡轮效率 | 0.90 | 0.92 | 0.910 |
| 6 | $\eta_{c1}$ | 内涵（一次）喷管效率 | 0.94 | 0.95 | 0.945 |
| 7 | $\eta_{c2}$ | 外涵（二次）喷管效率 | 0.92 | 0.94 | 0.930 |
| 8 | $\sigma_{cc}$ | 进气道激波总压恢复系数 | 0.98 | 1.00 | 0.990 |
| 9 | $\sigma_{kan}$ | 进气通道总压恢复系数 | 0.98 | 0.99 | 0.985 |
| 10 | $\sigma_{kask}$ | 压气机级间总压恢复系数 | 0.98 | 0.99 | 0.985 |
| 11 | $\sigma_{ks}$ | 燃烧室总压恢复系数 | 0.94 | 0.96 | 0.950 |
| 12 | $\eta_T$ | 燃烧放热系数 | 0.97 | 0.99 | 0.980 |
| 13 | $\lambda^*$ | 风扇涡轮热恢复系数 | 1.02 | 1.04 | 1.030 |

> **先验分布**：各参数独立，均匀分布于各自范围内，即 $p(\theta_i) \propto \mathbf{1}_{[lb_i,\, ub_i]}(\theta_i)$。

---

## 4. 单工况设置

程序内置一组典型的**起飞地面工况**：

| 工况参数 | 符号 | 值 | 单位 | 说明 |
|---|---|---|---|---|
| 环境温度 | $T_H$ | 288.15 | K | 标准大气，海平面 |
| 飞行马赫数 | $M_\pi$ | 0.0 | — | 地面静止起飞 |
| 涵道比 | $m$ | 10.0 | — | 大涵道比涡扇典型值 |
| 压气机总压比 | $\pi_k$ | 30.0 | — | 高涵道比涡扇典型值 |
| 涡轮前总温 | $T_g$ | 1500.0 | K | 落在 1400＜$T_g$≤1600 分段区间 |

$T_g = 1500\,\text{K}$ 的选取，确保分段函数 $k_T$、$R_T$ 激活中段取值（$k_T=1.30$，$R_T=288.0$）。

---

## 5. 前向热力循环模型

函数签名：`[y, aux] = engine_forward(theta, cond)`

### 5.1 公式链条总览

```
基本常数 (k_air, R_air)
    │
    ▼
声速 a, 飞行速度 V_flight
    │
    ├──► 分段函数：k_T(T_g), R_T(T_g), δ(T_g)
    │
    ▼
进口总压比 τ_v  →  压气机入口温度 T_B
    │
    ▼
压气机出口温度 T_k
    │
    ▼
相对耗油量 g_T
    │
    ▼
热恢复系数 λ_heat  ←  进口总压恢复系数 σ_bx
    │
    ▼
单位自由能 L_sv
    │
    ▼
最优自由能分配系数 x_pc
    │
    ▼
比推力 R_ud   →   比油耗 C_ud
```

### 5.2 各步公式详解

#### (1) 基本常数与飞行速度

$$
a = \sqrt{k \cdot R \cdot T_H}, \quad V_\pi = a \cdot M_\pi
$$

其中 $k=1.4$，$R=287.3\;\text{J/(kg·K)}$。

#### (2) 分段函数

**冷却引气系数：**
$$
\delta(T_g) = 0.02 + \frac{T_g - 1200}{100} \times 0.02
$$

**燃气绝热指数 $k_T(T_g)$：**

| 条件 | $k_T$ |
|---|---|
| $800 < T_g \leq 1400$ | 1.33 |
| $1400 < T_g \leq 1600$ | 1.30 |
| $T_g > 1600$ | 1.25 |

**燃气气体常数 $R_T(T_g)$（J/(kg·K)）：**

| 条件 | $R_T$ |
|---|---|
| $800 < T_g \leq 1400$ | 287.6 |
| $1400 < T_g \leq 1600$ | 288.0 |
| $T_g > 1600$ | 288.6 |

#### (3) 进口总压比与压气机入口温度

$$
\tau_v = \left(1 + \frac{V_\pi^2}{2 \cdot \frac{k}{k-1} \cdot R \cdot T_H}\right)^{\!\frac{k}{k-1}}
$$

$$
T_B = T_H \cdot \left(1 + \frac{V_\pi^2}{2 \cdot \frac{k}{k-1} \cdot R \cdot T_H}\right)
$$

起飞工况 $M_\pi=0$ 时，$\tau_v = 1$，$T_B = T_H$。

#### (4) 压气机出口温度

$$
T_k = T_B \cdot \left(1 + \frac{\pi_k^{\frac{k-1}{k}} - 1}{\eta_k}\right)
$$

#### (5) 相对耗油量

$$
g_T = 3.293\times10^{-5} \cdot T_g - 2.84\times10^{-5} \cdot T_k - 0.004814
$$

#### (6) 热恢复系数

$$
\lambda_\text{heat} = \frac{
  1 - \dfrac{\frac{k}{k-1} R T_B (\pi_k^{\frac{k-1}{k}}-1)}{\frac{k_T}{k_T-1} R_T T_g \cdot \eta_k \eta_T}
}{
  1 - \dfrac{\frac{k}{k-1} R T_B (\pi_k^{\frac{k-1}{k}}-1)}{\frac{k_T}{k_T-1} R_T T_g \cdot \eta_k \eta_T \eta_t}
}
$$

> 分子对应不含涡轮效率的理想情形，分母含涡轮效率；比值反映了热量回收程度。

#### (7) 进口总压恢复系数

$$
\sigma_{bx} = \sigma_{cc} \cdot \sigma_{kan}
$$

#### (8) 单位自由能

$$
L_{sv} = \lambda_\text{heat} \cdot \frac{k_T}{k_T-1} R_T T_g \left[1 - \left(\frac{1}{\tau_v \sigma_{bx} \pi_k \sigma_{kask} \sigma_{ks}}\right)^{\!\frac{k_T-1}{k_T}}\right]
- \frac{\frac{k}{k-1} R T_B (\pi_k^{\frac{k-1}{k}}-1)}{(1+g_T)\,\eta_k\,\eta_t\,\eta_m\,(1-\delta)}
$$

> 第一项为燃气膨胀可用功，第二项为驱动压气机消耗的功。

#### (9) 最优自由能分配系数

$$
x_{pc} = \frac{1 + \dfrac{m V_\pi^2}{2 L_{sv} \eta_{tv} \eta_v \eta_{c2}}}{1 + \dfrac{m\,\eta_{tv}\,\eta_v\,\eta_{c2}}{\eta_{c1}\,\lambda^*}}
$$

> $x_{pc}$ 决定了内涵/外涵之间的最优能量分配比例。

#### (10) 比推力

$$
R_{ud} = \frac{1}{1+m}\left[(1+g_T)\sqrt{2\eta_{c1}\lambda^* x_{pc} L_{sv}} - V_\pi\right]
+ \frac{m}{1+m}\left[\sqrt{\frac{2(1-x_{pc})}{m} L_{sv}\,\eta_{tv}\,\eta_v\,\eta_{c2} + V_\pi^2} - V_\pi\right]
$$

> 两项分别对应内涵（高温燃气）和外涵（旁路冷空气）的动量增量贡献。

#### (11) 比油耗

$$
C_{ud} = \frac{3600\; g_T\,(1-\delta)}{R_{ud}\,(1+m)}
$$

> 系数 3600 将单位从 kg/(N·s) 换算为 kg/(N·h)。

---

## 6. 贝叶斯反演框架

### 6.1 先验

均匀先验，各参数独立：

$$
\log p(\boldsymbol{\theta}) = \begin{cases} 0 & \text{若 } lb_i \le \theta_i \le ub_i \;\forall i \\ -\infty & \text{否则} \end{cases}
$$

对应 MATLAB 函数：`log_prior(theta, lb, ub)`

### 6.2 似然

观测误差服从独立高斯分布：

$$
R_{ud}^{obs} = R_{ud}(\boldsymbol{\theta}) + \varepsilon_R, \quad \varepsilon_R \sim \mathcal{N}(0,\,\sigma_R^2)
$$

$$
C_{ud}^{obs} = C_{ud}(\boldsymbol{\theta}) + \varepsilon_C, \quad \varepsilon_C \sim \mathcal{N}(0,\,\sigma_C^2)
$$

对数似然：

$$
\log \mathcal{L}(\boldsymbol{\theta}) = -\frac{1}{2}\left[\left(\frac{R^{obs}-R^{pred}}{\sigma_R}\right)^{\!2} + \ln(2\pi\sigma_R^2) + \left(\frac{C^{obs}-C^{pred}}{\sigma_C}\right)^{\!2} + \ln(2\pi\sigma_C^2)\right]
$$

噪声标准差由观测值的相对比例确定：

$$
\sigma_R = 1\% \times |R_{true}|, \quad \sigma_C = 1\% \times |C_{true}|
$$

对应 MATLAB 函数：`log_likelihood(theta, data, cond)`

### 6.3 后验

$$
\log p(\boldsymbol{\theta}\mid\mathbf{y}) = \log p(\boldsymbol{\theta}) + \log\mathcal{L}(\boldsymbol{\theta})
$$

对应 MATLAB 函数：`log_posterior_fn(theta, data, cond, lb, ub)`

---

## 7. MCMC采样器设计

函数签名：`results = run_mcmc(data, cond, lb, ub, theta0, opts)`

### 7.1 算法：随机游走 Metropolis-Hastings

在归一化空间 $\mathbf{z} = (\boldsymbol{\theta} - lb) \oslash (ub - lb) \in [0,1]^{13}$ 中进行采样：

1. **提案**：$\mathbf{z}' = \mathbf{z}^{(t)} + \sigma_{prop} \cdot \boldsymbol{\xi}$，$\boldsymbol{\xi} \sim \mathcal{N}(\mathbf{0}, \mathbf{I})$
2. **边界处理**：对越界分量使用**反射边界**（不截断）
3. **映射**：$\boldsymbol{\theta}' = lb + \mathbf{z}' \odot (ub - lb)$
4. **接受概率**：$\alpha = \min\!\left(1,\; e^{\log p(\boldsymbol{\theta}'|\mathbf{y}) - \log p(\boldsymbol{\theta}^{(t)}|\mathbf{y})}\right)$
5. **更新**：以概率 $\alpha$ 接受提案

在归一化空间中采样的好处：不同参数的量纲和范围差异被消除，单一标量 $\sigma_{prop}$ 即可适配所有参数。

### 7.2 反射边界

```
若 z < 0：z ← -z         （关于 0 镜像反射）
若 z > 1：z ← 2 - z      （关于 1 镜像反射）
```

反复折叠直至落入 $[0,1]$，保证采样在物理有效域内，且不破坏细致平衡条件。

### 7.3 自适应提案方差

在步数 `adapt_start`（500）到 `adapt_end`（2500）之间，每 100 步统计局部接受率，自动调整 $\sigma_{prop}$：

| 局部接受率 | 调整 |
|---|---|
| $< 0.18$（过低，步子太大） | $\sigma_{prop} \leftarrow 0.9\,\sigma_{prop}$ |
| $> 0.35$（过高，步子太小） | $\sigma_{prop} \leftarrow 1.1\,\sigma_{prop}$ |
| 在 $[0.18,\, 0.35]$ 之间 | 不调整 |

### 7.4 采样配置（轻量版）

| 配置项 | 值 | 说明 |
|---|---|---|
| 总采样步数 | 10,000 | 适合笔记本电脑 |
| 烧入期（burn-in） | 2,000 | 丢弃预热阶段 |
| 后验链长度 | 8,000 | 用于统计分析 |
| 初始提案 SD | 0.012 | 归一化空间中 |
| 自适应区间 | 步 500～2500 | 自动调节步长 |

### 7.5 输出结构体

| 字段 | 含义 |
|---|---|
| `chain_full` | 全链（含 burn-in），$10000 \times 13$ |
| `chain_post` | 后验链，$8000 \times 13$ |
| `logpost_full` | 全链对数后验值 |
| `accept_rate` | 总接受率 |
| `theta_mean` | 后验均值 |
| `theta_std` | 后验标准差 |
| `theta_map` | MAP 估计（后验链最大值对应点）|
| `theta_ci95` | 95% 后验区间，$13 \times 2$ |

---

## 8. 虚拟试验数据生成

函数签名：`data = generate_virtual_data(theta_true, cond, noise_level_R, noise_level_C, rng_seed)`

**流程：**
1. 以 `theta_true` 调用 `engine_forward`，得到 $R_{true}$，$C_{true}$
2. 计算噪声标准差：$\sigma_R = 0.01 \times |R_{true}|$，$\sigma_C = 0.01 \times |C_{true}|$
3. 加入高斯白噪声：$R_{obs} = R_{true} + \sigma_R \cdot \xi_R$，$\xi_R \sim \mathcal{N}(0,1)$

**目的：** 在真实参数已知的情况下，验证反演流程的正确性——若反演结果的后验均值（或 MAP）接近 `theta_true`，说明代码无误。

---

## 9. 输出结果与可视化

程序共输出 5 张图，自动保存为 PNG 文件：

### 图1：链轨迹图（`fig1_trace_plots.png`）

- 13 个参数各一子图，横轴为迭代步数
- 蓝线：全部链（含 burn-in）
- 绿线：真值 `theta_true`
- 红虚线：burn-in 结束位置
- **用途**：判断链是否收敛、是否混合良好

### 图2：边缘后验直方图（`fig2_marginal_posteriors.png`）

- 对 burn-in 后的 8000 个样本绘制直方图（概率密度归一化）
- 同时标注：真值（绿）、后验均值（红虚线）、MAP（品红点线）
- **用途**：直观观察每个参数的后验分布形态与展宽

### 图3：后验相关系数热图（`fig3_correlation_heatmap.png`）

- $13 \times 13$ 相关矩阵，由 `corrcoef(chain_post)` 计算
- 蓝-白-红色图：蓝色=负相关，红色=正相关
- 每格标注相关系数数值
- **用途**：识别参数间的强相关结构（高度相关的参数对难以同时被辨识）

### 图4：观测-预测对比图（`fig4_obs_pred_comparison.png`）

对比以下 4 组数值（分别对 $R_{ud}$ 和 $C_{ud}$）：

| 类别 | 说明 |
|---|---|
| 真值 | 调用前向模型 + `theta_true` |
| 观测值 | 含 1% 高斯噪声的虚拟观测（带 2σ 误差棒）|
| 后验均值预测 | 调用前向模型 + `theta_mean` |
| MAP 预测 | 调用前向模型 + `theta_map` |

### 图5：可辨识性指标（`fig5_identifiability.png`）

定义每个参数的可辨识性指标：

$$
\rho_i = \frac{\sigma_i^{post}}{\sigma_i^{prior}}
$$

其中 $\sigma_i^{prior} = (ub_i - lb_i)/\sqrt{12}$（均匀分布标准差）。

| 指标值 | 含义 |
|---|---|
| $\rho_i \approx 1$ | 数据对该参数几乎没有约束，后验≈先验 |
| $\rho_i \ll 1$ | 后验显著收缩，参数可辨识性强 |

---

## 10. 鲁棒性设计

### 10.1 前向模型防护

`engine_forward` 中对以下情况进行检测，发现问题则立即返回 `[NaN, NaN]`：

| 检测点 | 保护条件 |
|---|---|
| 声速 $a$ | `a > 0` |
| 内层计算项 `inner` | `inner > 0`（防止虚数）|
| $T_B$、$T_k$ | `> 0`（热力学温度为正）|
| $g_T$ | `> 0`（耗油量须为正）|
| 分母项 | `abs(denom) > 1e-10`（防止除零）|
| $L_{sv}$ | `> 0`（可用功须为正）|
| 开方项 | `>= 0`（防止虚数）|
| $R_{ud}$、$C_{ud}$ | `> 0`（物理正值）|

所有防护均用 `try-catch` 包裹，任何 MATLAB 运行时异常都被捕获并返回 `NaN`。

### 10.2 似然函数防护

`log_likelihood` 中：
- 检测前向模型返回值是否 `isfinite`
- 若参数样本无效，返回 `log_likelihood = -Inf`（对应后验概率为 0，MH 算法将自动拒绝该样本）

### 10.3 MCMC 防护

- 若初始点 `theta0` 的后验为 `-Inf`，自动切换到先验中点重试
- 若先验中点也无效，程序报错并给出诊断信息

---

## 11. 局部函数索引

| 函数名 | 输入 | 输出 | 位置 |
|---|---|---|---|
| `piecewise_kT(T_g)` | 涡轮前温度 | 燃气绝热指数 $k_T$ | 约第244行 |
| `piecewise_RT(T_g)` | 涡轮前温度 | 燃气气体常数 $R_T$ | 约第264行 |
| `delta_cooling(T_g)` | 涡轮前温度 | 冷却引气系数 $\delta$ | 约第280行 |
| `engine_forward(theta, cond)` | 参数向量，工况 | $[R_{ud},C_{ud}]$，中间量 | 约第296行 |
| `generate_virtual_data(...)` | 真值参数，工况，噪声 | 观测数据结构体 | 约第530行 |
| `log_prior(theta, lb, ub)` | 参数，下界，上界 | 对数先验 | 约第562行 |
| `log_likelihood(theta, data, cond)` | 参数，数据，工况 | 对数似然 | 约第578行 |
| `log_posterior_fn(...)` | 同上+先验边界 | 对数后验 | 约第617行 |
| `reflect_to_unit_interval(z)` | 归一化向量 | 反射后向量 | 约第631行 |
| `run_mcmc(...)` | 数据，工况，先验，初值，配置 | results 结构体 | 约第658行 |
| `plot_results(...)` | results 等 | 5张图 | 约第781行 |
| `txt_color(abs_val)` | 绝对相关系数值 | 颜色字符串 | 约第915行 |
| `redblue_colormap()` | — | 64级色图矩阵 | 约第927行 |

---

## 12. 常见问题

**Q1：程序报 "前向模型在 theta_true 处输出非有限值" 怎么办？**

检查 `cond.pi_k`（压比）和 `cond.T_g`（涡轮前温度）是否合理。`pi_k` 建议在 15～40 之间，`T_g` 须大于 800K。

**Q2：MCMC 接受率非常低（< 5%）怎么优化？**

减小 `opts.proposal_sd`，例如从 0.012 改为 0.005。由于问题高度欠定，接受率偏低（15%～25%）是正常的。

**Q3：后验分布与先验几乎相同（可辨识性指标≈1），是否正常？**

是正常的。单工况只有两个观测量，对 13 个参数的约束能力非常有限。部分参数（尤其是各效率参数间高度相关的组合）无法被单独辨识，这正是欠定问题的本质特征。

**Q4：如何增强约束（提高可辨识性）？**

扩展为多工况版本（多个 $\pi_k$、$T_g$ 组合）或增加更多观测量（如效率、温度分布），参数数目相同而方程数增多，欠定程度随之降低。

**Q5：增加采样量会明显改善后验质量吗？**

对于欠定问题，增加采样量（如 50000 步）能提高后验分布估计的精度，但不会从根本上改变参数不可辨识的情况——后验分布仍然宽泛，这是信息不足的物理限制，而非数值问题。

---

*本文档由程序自动说明文本生成，与 `bayesian_inversion_single_condition.m` 配套使用。*
