# 一维 Bose-Hubbard 模型 DMRG (vMPS) 相图计算 🌐[English Version](README_English.md)

![Mott叶瓣相图](figure_result.png)

本项目使用双格点密度矩阵重整化群（DMRG）/变分矩阵乘积态（vMPS）算法，结合有限尺寸标度，高精度计算一维 Bose-Hubbard 模型单位填充下的 Mott 绝缘体相边界，重现著名的尖峰状 Mott 叶瓣结构。
---

## 背景

一维 Bose-Hubbard 模型描述光学晶格中相互作用的玻色子，哈密顿量为

$$
H = -t \sum_{\langle i,j\rangle} (a_i^\dagger a_j + \mathrm{h.c.}) + \frac{U}{2} \sum_i n_i(n_i-1),
$$

其中 $t$ 为跃迁强度， $U$ 为同格点相互作用。本代码以 $U=1$ 作为能量单位， $t$ 为可调参数。在均匀填充 $n=1$ 附近，系统存在 Mott 绝缘相，其相边界在热力学极限下呈现尖端形状，超越平均场理论的定性抛物线。平均场自洽迭代只能给出定性mott叶瓣相图，未能给出一维情况下著名的mott尖峰形状和相边界。

遵循原始论文 *The one-dimensional Bose-Hubbard Model with nearest-neighbor interaction* 提出的方法，本项目通过粒子激发与空穴激发能量的有限尺寸标度，精确提取热力学极限下的化学势 $\mu_p(t)$ 与 $\mu_h(t)$ ，从而绘制高精度 Mott 叶瓣。

本代码解决的核心问题：  
1. 构造带粒子数惩罚项的 Bose-Hubbard 哈密顿量 MPO，确保指定激发态下粒子数守恒。  
2. 实现完整的双格点 vMPS (DMRG) 算法，包含环境收缩、粗粒化、Lanczos 对角化、SVD 截断及双向扫描。  
3. 并行扫描参数空间 $(t, L, \mathrm{excite})$ ，计算不同系统尺寸的基态能量，得到有限尺寸化学势。  
4. 对每个 $t$ ，通过 $\mu$ 与 $1/L$ 的线性外推获得热力学极限结果，绘制相边界。

---

## 方法 / 算法 （和代码顺序相匹配）

### 1. MPO 构造与粒子数约束

一维 Bose‑Hubbard 模型的哈密顿量可写为

$$
H_{\mathrm{BH}} = -t \sum_{\langle i,j\rangle} (a_i^\dagger a_j + \mathrm{h.c.}) + \frac{U}{2} \sum_i n_i(n_i-1)
$$

取 $U=1$ 作为能量单位。用 MPO 表示该哈密顿量时，每个格点的局部张量是一个 $4\times4$ 的矩阵，其矩阵元为 $d\times d$ 的算符：

$$
W_{\mathrm{BH}}(t) =
\begin{pmatrix}
I & a^\dagger & -t a & \frac{1}{2}n(n-1) \\
0 & 0 & 0 & -t a \\
0 & 0 & 0 & a^\dagger \\
0 & 0 & 0 & I
\end{pmatrix}
$$

其中 $I$ 是 $d$ 维单位矩阵， $a$ 和 $a^\dagger$ 是玻色子湮灭、产生算符（截断至 $n_{\max}=4$ ）， $n$ 是粒子数对角矩阵 $\mathrm{diag}(0,1,\dots,n_{\max})$ 。左边界张量取上述矩阵的第一行，右边界张量取最后一列。

为实现指定粒子数扇区的基态，我们在哈密顿量上添加惩罚项 $\lambda(\hat{N}-N_{\mathrm{target}})^2$ ，这里 $\hat{N}=\sum_i n_i$ 为总粒子数算符， $N_{\mathrm{target}}=L+\mathrm{excite}$ 为目标粒子数， $\lambda=10$ 为惩罚强度。展开惩罚项：

$$
\lambda(\hat{N}-N_{\mathrm{target}})^2 = \lambda\hat{N}^2 - 2\lambda N_{\mathrm{target}}\hat{N} + \lambda N_{\mathrm{target}}^2
$$

分别构造这三部分的 MPO 再通过直和相加。

- **$\lambda\hat{N}^2$ 项**：其 MPO 的局部张量为 $3\times3$ 的分块矩阵

$$
W_{\lambda N^2} =
\begin{pmatrix}
I & \sqrt{\lambda}n & \lambda n^2 \\
0 & I & 2\sqrt{\lambda}n \\
0 & 0 & I
\end{pmatrix}
$$

这种构造保证了乘积展开后恰好得到 $\lambda(\sum_i n_i)^2$ 。

- **$-2\lambda N_{\mathrm{target}}\hat{N}$ 项**：先写出 $\hat{N}$ 的 MPO 张量，它是 $2\times2$ 的：

$$
W_N =
\begin{pmatrix}
I & n \\
0 & I
\end{pmatrix}
$$

于是该项的局部张量为

$$
W_{cN} =
\begin{pmatrix}
I & -2\lambda N_{\mathrm{target}} n \\
0 & I
\end{pmatrix}
$$

- **常数项 $\lambda N_{\mathrm{target}}^2$** ：对应一个辅助维数为 1 的恒等 MPO，即局部张量 $W_c = \lambda N_{\mathrm{target}}^2 I$ （作为 $1\times1$ 矩阵）。

要将它们合成总的惩罚 MPO，我们使用 **MPO 加法（直和）** 。对于两个中间格点的张量 $W_A$ （形状 $(a,b,d,d)$ ）和 $W_B$ （形状 $(c,e,d,d)$ ），它们的和为形状 $(a+c, b+e, d,d)$ 的分块对角张量：

$$
(W_A \oplus W_B)_{\alpha\beta} =
\begin{cases}
(W_A)_{\alpha\beta}, & \alpha<a,\;\beta<b\\
(W_B)_{(\alpha-a)(\beta-b)}, & \alpha\ge a,\;\beta\ge b\\
0, & \text{其它}
\end{cases}
$$

左右边界的张量加法需保证外部辅助维数始终为 1，具体为：左边界张量形状 $(1, r_A)$ 与 $(1, r_B)$ 相加得到 $(1, r_A+r_B)$ ；右边界张量形状 $(l_A, 1)$ 与 $(l_B, 1)$ 相加得到 $(l_A+l_B, 1)$ 。

按照此规则，先将 $-2\lambda N_{\mathrm{target}}\hat{N}$ 项（维数 2）与常数项（维数 1）相加，得到一个辅助维数为 $2+1=3$ 的 MPO；再将其与 $\lambda\hat{N}^2$ 项（维数 3）相加，得到辅助维数为 $3+3=6$ 的总惩罚 MPO。最后，将 Bose‑Hubbard MPO（维数 4）与总惩罚 MPO 相加，得到 **完整哈密顿量 MPO**，其辅助维数为 $4+6=10$ ，对应算符

$$
H_{\mathrm{tot}} = H_{\mathrm{BH}} + \lambda(\hat{N}-N_{\mathrm{target}})^2\; .
$$

足够大的 $\lambda$ 保证 DMRG 优化后总粒子数精确锁定在 $N_{\mathrm{target}}$ 。

---

### 2. MPS 初始化与激发态

矩阵乘积态（MPS）将多体波函数表示为每个格点上的三阶张量 $A^{s_i} _ {\alpha_i\alpha_{i+1}}$ ，其中 $s_i\in\{0,1,\dots,n_{\max}\}$ 为物理指标， $\alpha_i$ 为辅助指标（bond 维数，初始为 1）。我们希望从目标填充数对应的直积态出发，经由 DMRG 迭代收敛到真正的本征态。

- **基态（ $\mathrm{excite}=0$ ）**：每个格点占据数恰好为 1。对应的 MPS 张量为形状 $(d,1,1)$ ，只有 $s=1$ 的分量为 1：

$$
A^{s} = \delta_{s,1} \quad\Rightarrow\quad A^{1}_{0,0}=1,\; A^{s\neq 1}_{0,0}=0
$$

此即直积态 $|\psi_0\rangle = |1\rangle^{\otimes L}$ 。

- **粒子激发（ $\mathrm{excite}=+1$ ）**：总粒子数 $N_{\mathrm{target}} = L+1$ 。在系统中心格点 $i = \lfloor L/2 \rfloor$ 放置一个额外粒子，使其占据数变为 2，其余格点仍为 1。中心格点张量为

$$
A^{s} = \delta_{s,2} \quad\Rightarrow\quad A^{2}_{0,0}=1
$$

其他格点同基态。该状态对应一条额外的粒子激发。

- **空穴激发（ $\mathrm{excite}=-1$ ）**：总粒子数 $N_{\mathrm{target}} = L-1$ 。中心格点占据数设为 0，即

$$
A^{s} = \delta_{s,0} \quad\Rightarrow\quad A^{0}_{0,0}=1
$$

其余格点仍为 $|1\rangle$ ，对应一个空穴激发。

所有这些初始 MPS 的 bond 维数均为 1。在后续的双格点 DMRG 优化中，bond 维数会根据截断容差 `target_trunc=1e-8` 和上限 $D=625$ 动态增长，以描述量子涨落与纠缠。由于哈密顿量带有粒子数惩罚项，即使从直积态出发， DMRG 也会被引导到具有正确粒子数的基态，保证了三种激发扇区的计算一致性。

---

### 3. 双格点 DMRG (vMPS) 流程（更详细的入门讲解见Schollwöck的著名综述）

现代 DMRG 可理解为一种张量网络算法：将多体哈密顿量和态矢分解为局域在格点上的小张量（分别称为MPO和MPS），通过张量收缩计算期望值（这画出来像网络）。于是求解全局基态等价于依次访问每个格点的MPS小张量，将其附近的"环境"（即除它以外的所有MPO和MPS）收缩为“有效哈密顿量”，并求其最小本征矢作为该格点的优化张量。

双格点算法允许动态调整相邻态矢小张量之间收缩的求和维数。它每步同时优化相邻两个格点合并后的大张量，优化后对大张量做 SVD 分解并丢弃奇异值较小的部分，这使得要求和的张量值变少了。这一步与量子纠缠直接相关——被舍弃的正是对当前态贡献极小的纠缠自由度。具体流程包括：
- **环境张量** ：从边界真空开始，依次收缩 MPS 与 MPO 构建左、右环境，用于后续有效哈密顿量的快速构建。
- **粗粒化** ：将相邻两格点的 MPS 合并为双格点张量，相邻 MPO 合并为双格点有效 MPO。
- **有效哈密顿量对角化** ：构造 `LinearOperator`，实现双格点有效哈密顿量对波函数的乘法，利用 `scipy.sparse.linalg.eigsh` (Lanczos) 求解最小本征值与对应的本征态，得到当前迭代的最佳近似基态。
- **细粒化与截断** ：对双格点波函数进行 SVD 分解，根据奇异值平方累计和与截断容差 `target_trunc=1e-8` 以及最大保留 bond 维数 `D=625` 自适应截断。奇异值按扫描方向分配给左或右张量，维持 MPS 的正则形式。
- **双向扫描与环境更新** ：从左至右、从右至左交替优化所有相邻格点对，每次优化后更新相应的环境张量。重复扫描直至基态能量变化小于 `energy_tol=1e-7` 或达到最大扫描次数。

---

### 4. 化学势与有限尺寸外推

代码实际扫描的参数为 $t = 0.01,0.02,\dots,0.30$ （步长 0.01）、 $L=32,64,128$ 以及 $\mathrm{excite}=0,+1，-1$ 。对每个 $(t, L)$ 计算三种激发态的基态能量

- 对每个 $(t, L)$ 计算三种激发态的基态能量，得到有限尺寸化学势 $\mu_p(L)=E(L+1)-E(L)$ , $\mu_h(L)=E(L)-E(L-1)$ 。
- 固定 $t$ ，将 $\mu_p(L)$ 和 $\mu_h(L)$ 对 $1/L$ 做线性拟合，外推到热力学极限 ( $1/L \to 0\Leftrightarrow L\to\infty$ )，截距作为 $\mu_p(t)$ 和 $\mu_h(t)$ 。
- 所有参数组合通过 `multiprocessing.Pool` 并行计算，极大缩短运行时间。中间结果保存为 `pickle` 文件以支持断点续算。

### 5. 主要可调参数

| 参数 | 取值 | 说明 |
|------|------|------|
| 最大 bond 维数 $D$ | 625 | 截断保留态数上限 |
| 单格点最大占据数 $n_{\max}$ | 4 | 局部希尔伯特空间截断 |
| 惩罚系数 $\lambda$ | 10 | 粒子数约束强度 |
| 能量收敛判据 | $10^{-7}$ | DMRG 扫描停止条件 |
| SVD 截断目标 | $10^{-8}$ | 奇异值累计舍弃阈值 |

---

## 代码结构

- `DMRG_BH_ycr.py` ：主程序脚本，包含
  - MPO 生成（Bose-Hubbard、粒子数惩罚、加法）
  - MPS 初始化函数
  - 双格点 DMRG 算法核心（环境收缩、粗粒化、Lanczos 对角化、SVD 截断、扫描）
  - 化学势计算、有限尺寸外推
  - 多进程并行调度与主函数 `main()`
- `energy_results_lobe1&2/KT.pkl` ：计算过程中保存的基态能量字典（自动生成）
- `DMRG_BH_ycr_Plot.py` ：生成精美图像
- `figure_result.png` ：输出的 Mott 叶瓣相图（高分辨率）
- `tensor_basic_Note.ipynb` ：练习笔记，如何用numpy实现张量（高阶nd数组）的收缩、重塑等操作，这在张量网络算法中极其重要

---

## 依赖环境

运行代码需要以下 Python 库：

- `numpy`
- `scipy`
- `matplotlib`
- `multiprocessing` (内置)
- `pickle` (内置)
- `time` (内置)

建议使用 Python 3.7 及以上版本。

---

## 快速开始

1. **克隆仓库**（或下载文件）：
   ```bash
   git clone https://github.com/chaoranyang/QuantumManyBody_FromZreo.git
   cd QuantumManyBody_FromZreo/BH_DMRG
2. **安装依赖**：
   ```bash
   pip install numpy matplotlib SciencePlots
---

## 结果输出

- **控制台输出**：显示每组 $(t,L,\mathrm{excite})$ 的初始粒子数、迭代轮次、基态能量、最终粒子数以及真实哈密顿量期望值，并统计总耗时。
- **相图图片**： `figure_result.png` 展示 $\mu_p(t)$ （粒子激发）与 $\mu_h(t)$ （空穴激发）随跃迁强度 $t$ 的变化。两条边界曲线围成尖锐的 Mott 叶瓣，与前期文献数值结果一致，成功捕捉平均场理论无法给出的尖峰形状。

---

## 参考资料

- Till D. Kühner, Steven R. White, and H. Monien. One-dimensional Bose-Hubbard model with nearest-neighbor interaction. *Phys. Rev. B*, 61:12474–12489, May 2000.
- Ulrich Schollwöck. The density-matrix renormalization group in the age of matrix product states. *Annals of Physics*, 326(1):96–192, January 2011.
