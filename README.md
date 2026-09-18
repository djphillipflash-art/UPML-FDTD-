# FDTD-UPML 复现

> 基于论文 *Stability analysis of FDTD to UPML for time dependent Maxwell equations* (Fang & Ying, 2009) 的 MATLAB 数值复现

![MATLAB](https://img.shields.io/badge/MATLAB-R2018b+-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-completed-success)

---

## 项目简介

用 MATLAB 复现了论文中的二维 TM 模式 Maxwell 方程 FDTD-UPML 求解器，并通过数值实验验证了截断域上 Yee 格式的长时间稳定性。

本仓库包含：
- 完整的 FDTD-UPML 求解器（内部域 + 单轴各向异性完美匹配层）
- 论文 Figure 3 的复现结果
- 复现过程中关键问题的调试记录（Yee 网格半格对齐）

---

## 项目结构

```
FDTD-UPML-Reproduction/
├── UPML_FDTD_FANG_V2.m       # 主程序（全部逻辑在一个脚本里）
├── README.md                 # 本文件
├── results/
│   └── stability_curve.png   # 稳定性曲线输出
└── docs/
    └── debugging_notes.md    # 关键调试记录
```

---

## 快速开始

### 环境要求

- MATLAB R2018b 或更高
- 无需额外工具箱

### 运行

```matlab
% 直接打开脚本运行
>> UPML_FDTD_FANG_V2
```

运行时间约 3–5 分钟（M = 10000，桌面级 CPU）。

**快速测试**（1 秒内完成）：将脚本开头的 M 改为 100，I2 = J2 改为 20。

---

## 物理模型

二维 TM 模式 Maxwell 方程组（ε = μ = 1）：

    ∂H₁/∂t =  ∂E₃/∂x₂
    ∂H₂/∂t =  ∂E₃/∂x₁
    ∂E₃/∂t =  ∂H₂/∂x₁ + ∂H₁/∂x₂

**截断域结构**：

| 区域 | 含义 | 使用变量 |
|------|------|---------|
| Ω₁ | 内部物理域 | E₃, H₁, H₂ |
| Ω₂ | x 方向 PML | e₃, H₁, h₂ |
| Ω₃ | y 方向 PML | e₃, h₁, H₂ |
| Ω₄ | 角 PML | e₃, h₁, h₂ |

---

## 数值方法

### Yee 网格（时空交错）

```
时间:  ──●────○────●────○────●──
        E3   H1,H2  E3   H1,H2 E3
        n    n+1/2  n+1  n+3/2 ...
```

| 分量 | 空间位置 | 时间层 |
|------|---------|--------|
| E₃ | (i, j) | n |
| H₁ | (i, j+1/2) | n+1/2 |
| H₂ | (i+1/2, j) | n+1/2 |

### 蛙跳推进顺序

    E3^{n-1} ──► E3^n ──► e3^n ──► H^{n+1/2} ──► h^{n+1/2}

### UPML 辅助变量递推

论文公式 (3.3) 给出封闭求和形式，本项目采用等价递推形式（O(1) 复杂度）：

    e₃ⁿ = (β/α) · e₃ⁿ⁻¹ + (1/α) · (E₃ⁿ − E₃ⁿ⁻¹)

其中 α = 1 + στ/2，β = 1 − στ/2。

---

## 参数配置（论文例 1）

| 参数 | 值 | 说明 |
|------|-----|------|
| σ | 4.0 | PML 衰减系数 |
| τ | 0.005 | 时间步长 |
| l₁ = l₂ | 0.1 | 空间步长 |
| M | 10000 | 迭代步数 |
| I₀ = J₀ | 10 | 散射体半宽 |
| I₁ = J₁ | 50 | PML 内边界 |
| I₂ = J₂ | 55 | PML 外边界（5 层） |
| k₁ = k₂ | 1 | 初始模式阶数 |
| w | 4.5 | 初始模式系数 |

**稳定性检查**：
- CFL 条件：τ/l = 0.05 ≪ 1/√2 ≈ 0.707  ✅
- 论文条件 (5.12)：满足  ✅

---

## 调试记录（重要）

### 问题：长时间范数指数爆炸

初版实现中 l²(Ω₁) 范数在长时间步进后发散。

**根因**：Yee 网格半格错位。计算 E₃(i,j) 处 ∂H₂/∂x₁ 时，误用了 (i+1/2) 和 (i+3/2) 两点，空间错位半个格子：

```matlab
% 错误：空间错位半格
h2_right = H2(i+1 + I2 + 1, idx_j);
h2_left  = H2(i   + I2 + 1, idx_j);

% 正确：正确的中心差分
h2_right = H2(i   + I2 + 1, idx_j);   % (i+1/2, j)
h2_left  = H2(i-1 + I2 + 1, idx_j);   % (i-1/2, j)
```

**后果**：半格错位引入一阶空间误差 → 激发高频数值噪声 → Yee 格式弱稳定性下指数放大 → 范数爆炸。

**修正后**：范数曲线与论文 Figure 3 一致，二阶精度恢复。

**结论**：Yee 格式对空间错位极其敏感，是 FDTD 实现中最容易出错、也最致命的细节。

---

## 结果

### 稳定性曲线（论文 Figure 3 复现）

运行脚本后自动输出 l²(Ω₁) 范数随时间变化：

    Time (t)  ──────────────────────────────────────►
              ╱╲
    Norm     ╱  ╲
            ╱    \
           ╱.     \_____________________(stable)
          ╱

- 初始阶段短暂上升：波从内部域向 PML 传播
- 之后趋于平稳：PML 吸收效果良好
- 长时间不发散：验证论文理论结论

---

## 后续可扩展方向

- [ ] 渐变 PML：σ(s) = σ_max · (s/d)² 二次渐变，进一步减少界面反射
- [ ] TE 模式：交换 E/H 角色，验证论文对称性
- [ ] 任意散射体：三角形、圆形（论文例 3）
- [ ] 无限域对比：与无 PML 标准 FDTD 结果比较相对误差
- [ ] 向量化加速：meshgrid + 切片替代三重循环，预计提速 10–50×
- [ ] GPU 加速：gpuArray 加速差分运算

欢迎 PR！

---

## 参考文献

1. Fang N. S., Ying L. A. *Stability analysis of FDTD to UPML for time dependent Maxwell equations*. Science in China Series A, 2009, 52(4): 794–816.
2. Yee K. S. *Numerical solution of initial boundary value problems involving Maxwell's equations*. IEEE Trans. Antennas Propag., 1966, 14: 302–307.
3. Bérenger J. P. *A perfectly matched layer for the absorption of electromagnetic waves*. J. Comput. Phys., 1994, 114: 185–200.
4. Taflove A., Hagness S. C. *Computational Electrodynamics: The Finite-Difference Time-Domain Method*. 3rd ed. Artech House, 2005.

---

## License

MIT © 2024 [Your Name]

---

## 致谢

感谢论文作者 Fang N. S. 与 Ying L. A. 提供的理论框架与数值算例。

如果这个仓库对你有帮助，欢迎 star ⭐！
