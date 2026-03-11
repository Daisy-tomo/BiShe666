# 单工况贝叶斯反演代码说明

**文件名**：`bayesian_inversion_single_condition.m`
**对应多工况版本**：`bayesian_inversion_multicondition.m`
**日期**：2026-03

---

## 一、程序目的

本程序在**单一飞行工况**下，对双涵道涡扇发动机的 **13 个修正参数**进行贝叶斯反演，用于分析：

- 单工况下各参数的**可辨识性**（后验是否收缩）
- 参数之间的**后验相关性**
- 与多工况版本的**对比基准**

---

## 二、单工况设定

| 工况变量 | 符号 | 数值 | 说明 |
|----------|------|------|------|
| 增压比 | $\pi_k$ | **33** | 高压压气机总压比 |
| 涡轮前温度 | $T_g$ | **1800 K** | 燃烧室出口总温 |
| 飞行马赫数 | $M$ | **0** | 地面静止状态 |
| 飞行速度 | $V$ | **0 m/s** | 由马赫数计算得 |
| 环境温度 | $T_H$ | **288 K** | 标准海平面大气 |
| 涵道比 | $m$ | **8** | 外涵/内涵质量流量比 |

> **与多工况版本的区别**：多工况版本使用 8 个不同工况联合约束参数，本版本仅使用 1 个工况，观测信息量更少，预期参数可辨识性更弱。

---

## 三、13 个待识别参数

| 序号 | 参数名 | 物理含义 | 先验下界 | 先验上界 | 真值 |
|------|--------|----------|---------|---------|------|
| 1 | $\eta_k$ | 压气机等熵效率 | 0.84 | 0.92 | 0.90 |
| 2 | $\eta_t$ | 涡轮等熵效率 | 0.86 | 0.92 | 0.89 |
| 3 | $\eta_m$ | 机械传动效率 | 0.980 | 0.995 | 0.990 |
| 4 | $\eta_v$ | 风扇等熵效率 | 0.85 | 0.87 | 0.860 |
| 5 | $\eta_{tv}$ | 风扇涡轮效率 | 0.90 | 0.92 | 0.910 |
| 6 | $\eta_{c1}$ | 内涵喷管效率 | 0.94 | 0.95 | 0.945 |
| 7 | $\eta_{c2}$ | 外涵喷管效率 | 0.92 | 0.94 | 0.930 |
| 8 | $\sigma_{cc}$ | 进气道总压恢复系数 | 0.98 | 1.00 | 0.990 |
| 9 | $\sigma_{kan}$ | 进气通道总压恢复 | 0.98 | 0.99 | 0.985 |
| 10 | $\sigma_{kask}$ | 压气机级间总压恢复 | 0.98 | 0.99 | 0.985 |
| 11 | $\sigma_{ks}$ | 燃烧室总压恢复 | 0.94 | 0.96 | 0.950 |
| 12 | $\eta_T$ | 燃烧放热系数 | 0.97 | 0.99 | 0.980 |
| 13 | $\lambda^*$ | 风扇涡轮热恢复系数 | 1.02 | 1.04 | 1.030 |

先验分布均为**均匀分布** $\boldsymbol{\theta} \sim \mathcal{U}(\text{lb}, \text{ub})$。

---

## 四、贝叶斯框架

### 4.1 贝叶斯公式

$$
\underbrace{p(\boldsymbol{\theta} \mid \text{obs})}_{\text{后验}} = Z^{-1} \underbrace{p(\text{obs} \mid \boldsymbol{\theta})}_{\text{似然}} \cdot \underbrace{p(\boldsymbol{\theta})}_{\text{先验}}
$$

其中 $Z = \int p(\text{obs}|\boldsymbol{\theta})\,p(\boldsymbol{\theta})\,\mathrm{d}\boldsymbol{\theta}$ 为模型证据（归一化常数）。

### 4.2 似然函数

单工况下，观测量为比推力 $R_{ud}$ 和比油耗 $C_{ud}$，假设测量噪声独立且服从高斯分布：

$$
\ln p(\text{obs} \mid \boldsymbol{\theta}) = -\frac{1}{2}\left(\frac{R_{\text{obs}} - R_{\text{pred}}(\boldsymbol{\theta})}{\sigma_R}\right)^2 - \frac{1}{2}\left(\frac{C_{\text{obs}} - C_{\text{pred}}(\boldsymbol{\theta})}{\sigma_C}\right)^2 - \ln(2\pi\sigma_R^2) - \ln(2\pi\sigma_C^2)
$$

其中：
- $R_{\text{pred}}(\boldsymbol{\theta})$、$C_{\text{pred}}(\boldsymbol{\theta})$ 由前向热力循环模型计算
- $\sigma_R = 1\% \times R_{\text{true}}$，$\sigma_C = 1\% \times C_{\text{true}}$（相对噪声水平）

> **与多工况的区别**：多工况版本对 8 个工况的对数似然求和，信息量是单工况的 8 倍（理论上）。单工况下似然约束力弱，参数后验更宽。

### 4.3 前向模型

前向模型 `engine_forward(theta, cond)` 按照双涵道涡扇发动机热力循环公式，逐步计算：

```
T_H, M_flight, pi_k, T_g, m
        ↓
① 飞行速度 V = M × a_sound       (M=0 时 V=0)
② 进气道总压比 τ_v
③ 压气机入口总温 T_B = T_H       (M=0 时 T_B = T_H)
④ 压气机出口温度 T_k
⑤ 相对耗油量 g_T
⑥ 热恢复系数 λ_heat
⑦ 单位自由能 L_sv
⑧ 最优自由能分配系数 x_pc
⑨ 比推力 R_ud  [N·s/kg]
⑩ 比油耗 C_ud  [kg/(N·h)]
```

**地面静止工况（M=0）的特殊性**：
- 飞行速度 $V = 0$，进口动压为零
- $\tau_v = 1$（无飞行冲压效应）
- $T_B = T_H = 288\,\mathrm{K}$（压气机入口即大气温度）
- $x_{pc}$ 分子中 $mV^2/(2L_{sv}\cdot\ldots)$ 项为零，公式退化为简洁形式

---

## 五、MCMC 采样算法

### 5.1 算法类型

**随机游走 Metropolis-Hastings（RW-MH）**，在归一化参数空间 $z_j \in [0,1]$ 中执行：

$$
z_j = \frac{\theta_j - \text{lb}_j}{\text{ub}_j - \text{lb}_j}
$$

### 5.2 提议分布

$$
z_j^* = z_j^{(s)} + \epsilon, \quad \epsilon \sim \mathcal{N}(0, \sigma_{\text{prop}}^2)
$$

越界时采用**反射边界**处理（非截断），保持提议分布的对称性。

### 5.3 接受准则

$$
\alpha = \min\left(1,\; \frac{p(\boldsymbol{\theta}^* \mid \text{obs})}{p(\boldsymbol{\theta}^{(s)} \mid \text{obs})}\right)
$$

以概率 $\alpha$ 接受新点。由于均匀先验在支撑内为常数，接受比等价于似然比：

$$
\ln\alpha = \ln p(\text{obs}|\boldsymbol{\theta}^*) - \ln p(\text{obs}|\boldsymbol{\theta}^{(s)})
$$

### 5.4 自适应步长

在 burn-in 期间（第 500~3000 步），每隔 200 步根据局部接受率自动调整 $\sigma_{\text{prop}}$：

| 接受率 | 动作 |
|--------|------|
| $< 0.20$ | $\sigma_{\text{prop}} \leftarrow \sigma_{\text{prop}} / 1.10$（缩小步长） |
| $> 0.35$ | $\sigma_{\text{prop}} \leftarrow \sigma_{\text{prop}} \times 1.10$（扩大步长） |
| $[0.20, 0.35]$ | 不变（目标区间） |

### 5.5 采样配置

| 参数 | 轻量模式（默认） | 高精度模式 |
|------|----------------|-----------|
| 总采样步数 | 12,000 | 50,000 |
| Burn-in 步数 | 3,000 | 10,000 |
| 后验样本数 | 9,000 | 40,000 |
| 初始步长 $\sigma_{\text{prop}}$ | 0.015 | 0.012 |

---

## 六、程序结构

```
bayesian_inversion_single_condition.m
│
├── 主程序
│   ├── 第一部分：参数名称、上下界、真值、初值定义
│   ├── 第二部分：单工况定义（pi_k=33, T_g=1800K, M=0, T_H=288K）
│   ├── 第三部分：生成虚拟试验数据
│   ├── 第四部分：MCMC 参数配置
│   ├── 第五部分：运行 MCMC 采样
│   ├── 第六部分：输出后验统计（均值、MAP、95%CI、标准差）
│   ├── 第七部分：后验预测验证（回代检查）
│   └── 第八部分：绘图（5张图）
│
└── 局部函数
    ├── piecewise_kT(T_g)           — 燃气多变指数分段函数
    ├── piecewise_RT(T_g)           — 燃气气体常数分段函数
    ├── delta_cooling(T_g)          — 涡轮冷却引气系数
    ├── engine_forward(theta, cond) — 热力循环前向模型
    ├── generate_virtual_data(...)  — 生成带噪声虚拟观测
    ├── log_prior(theta, lb, ub)    — 对数先验（均匀分布）
    ├── log_likelihood(theta, ...)  — 对数似然（单工况高斯）
    ├── log_posterior_fn(...)       — 对数后验 = 先验 + 似然
    ├── reflect_to_unit_interval(z) — 反射边界处理
    ├── run_mcmc(...)               — MH-MCMC 主循环
    ├── plot_results(...)           — 绘制5张结果图
    └── redblue_cmap()              — 红白蓝颜色图辅助函数
```

---

## 七、输出图表说明

| 图编号 | 图名 | 内容 |
|--------|------|------|
| 图1 | MCMC 链轨迹图 | 13个参数的采样链，红线=burn-in边界，绿线=真值 |
| 图2 | 边缘后验分布 | 13个参数的后验直方图，叠加真值/均值/MAP |
| 图3 | 参数相关系数热图 | 13×13 后验相关性矩阵（红=正相关，蓝=负相关） |
| 图4 | 观测与预测对比 | $R_{ud}$、$C_{ud}$ 的观测值 vs 后验均值 vs MAP 柱状对比 |
| 图5 | 参数可辨识性指标 | 后验标准差/先验标准差（越小=可辨识性越强） |

---

## 八、与多工况版本的关键差异

| 对比项 | 单工况版本 | 多工况版本 |
|--------|-----------|-----------|
| 工况数量 | 1 | 8 |
| 观测数据量 | 2个（$R_{ud}$, $C_{ud}$） | 16个 |
| 对数似然项数 | 2项 | 16项 |
| 参数约束力 | 弱（信息少） | 强（信息多） |
| 预期可辨识性 | 部分参数后验≈先验 | 更多参数后验收缩 |
| `conds` 数据结构 | 单个结构体 `cond` | 结构体数组 `conds(1..8)` |
| `log_likelihood` 循环 | 无循环（1次前向计算） | 8次循环累加 |

---

## 九、运行方法

在 MATLAB 命令窗口直接运行：

```matlab
bayesian_inversion_single_condition
```

或在 MATLAB 编辑器中打开文件后按 **F5**（运行）。

**预期输出**：
1. 命令窗口打印参数初始化表、工况设定、虚拟数据、MCMC 进度、后验统计、预测对比
2. 自动生成5张图窗

---

## 十、扩展建议

1. **增加工况**：将 `cond` 替换为 `conds` 数组，参考 `bayesian_inversion_multicondition.m` 扩展为多工况联合反演，以提高可辨识性。
2. **修改工况参数**：直接修改第二部分的 `cond.pi_k`、`cond.T_g` 等字段即可更换工况。
3. **高精度采样**：取消注释第四部分"高精度模式"代码块，将采样步数提升至 50,000。
4. **非均匀先验**：将 `log_prior` 函数替换为正态分布或截断正态分布，以引入更多先验知识。
