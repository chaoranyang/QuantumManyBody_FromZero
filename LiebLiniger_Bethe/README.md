# Bethe ansatz 求解一维 Lieb-Liniger 模型状态方程

本项目通过数值求解 Bethe ansatz 积分方程，计算一维玻色气体（Lieb‑Liniger 模型）的基态能量、化学势等热力学量，并与弱耦合与强耦合极限下的渐近结果进行对比验证。

---

## 物理背景

均匀一维玻色气体由 Lieb‑Liniger 模型描述，其哈密顿量为

$$
H = -\frac{\hbar^2}{2m}\sum_{i=1}^N \frac{\partial^2}{\partial x_i^2} + g\sum_{i<j}\delta(x_i-x_j),
$$

其中 $g$ 为相互作用强度，$m$ 为粒子质量。无量纲相互作用参数定义为

$$
\gamma = \frac{mg}{\hbar^2 n},
$$

这里 $n=N/L$ 为粒子数密度。

在热力学极限下，基态可由 Bethe ansatz 化为关于准动量分布 $G(q)$ 的线性积分方程：

$$
\alpha = \gamma \int_{-1}^{1} dq \, G(q),
$$

$$
G(q) = \frac{1}{2\pi} + \int_{-1}^{1} \frac{dq'}{2\pi} \, G(q') \frac{2\alpha}{(q'-q)^2 + \alpha^2},
$$

其中 $\alpha(\gamma)$ 是一个由自洽条件确定的参数。基态能量密度 $e(\gamma)$ 由下式给出：

$$
e(\gamma) = \left( \frac{\gamma}{\alpha(\gamma)} \right)^3 \int_{-1}^{1} dq \, G(q;\gamma) \, q^2.
$$

由此可得化学势

$$
\mu = \frac{\hbar^2}{2m} n^2 \bigl[ 3e(\gamma) - \gamma e'(\gamma) \bigr].
$$

本代码解决的核心问题：  
1. 对给定的 $\gamma$ 序列自洽求解 $\alpha$ 和 $G(q)$；  
2. 计算 $e(\gamma)$ 和 $\mu(\gamma)$ 并绘图，与已知的强耦合（$\gamma \to \infty$）和弱耦合（$\gamma \to 0$）渐近表达式进行比较。

---

## 数值方法

### 1. 自洽迭代求解 Bethe 方程

- **离散化**：将 $q$ 在 $[-1,1]$ 上均匀划分为 200 个格点，网格间距 $\Delta q$。
- **初值**：  
  - 第一个 $\gamma$ 值（最大 $\gamma$）采用均匀初值 $G(q)=1/(2\pi)$ 和较大的 $\alpha=1000$。  
  - 后续 $\gamma$ 值使用前一个 $\gamma$ 的收敛结果作为初始猜测，加速迭代。
- **迭代更新**：  
  - 利用当前 $G$ 和 $\alpha$ 计算被积函数 $f(q,q') = \frac{G(q')}{\pi} \frac{\alpha}{(q'-q)^2+\alpha^2}$。  
  - 用梯形法则更新 $G(q) = \frac{1}{2\pi} + \int dq' f(q,q')$。  
  - 用新的 $G(q)$ 积分更新 $\alpha = \gamma \int dq \, G(q)$。  
- **收敛判据**：当 $G$ 和 $\alpha$ 的两次迭代值差的平方和（对 $G$）及绝对差（对 $\alpha$）均小于 $10^{-10}$ 时停止迭代，或达到最大迭代次数（默认 1000）。

### 2. 基态能量与化学势的计算

- 基态能量 $e(\gamma)$ 通过对 $G(q) q^2$ 进行梯形积分得到。
- 化学势中的导数 $e'(\gamma) = de/d\gamma$ 采用**非均匀网格中心差分**：
  - 对于内部点：$f'(x_i) \approx \frac{h_1^2 f_{i+1} + (h_2^2-h_1^2) f_i - h_2^2 f_{i-1}}{h_1 h_2 (h_1+h_2)}$，其中 $h_1=x_i-x_{i-1}$，$h_2=x_{i+1}-x_i$。
  - 边界点使用前向/后向差分。
- 最终无量纲化学势 $\tilde\mu = 3e - \gamma e'$，图中坐标定义为
  $$
  x = \frac{\tilde\mu}{2\gamma^2}, \qquad y = \frac{1}{\gamma},
  $$
  并分别对应物理量 $(\hbar^2/m)\mu/g^2$ 和 $(\hbar^2/m)n/g$。

### 3. 渐近曲线

在双对数图上绘制四条理论渐近线：

- $\gamma \to \infty$（强耦合，费米化极限）:
  $$ n = \sqrt{\frac{2m\mu}{\pi^2\hbar^2}} \quad\Rightarrow\quad y = \sqrt{\frac{2x}{\pi^2}} $$

- $\gamma \gg 1$（强耦合领头阶+次领头阶）:
  $$ n = \sqrt{\frac{2m\mu}{\pi^2\hbar^2}} + \frac{8\mu}{3\pi g} - \frac{2\sqrt{2}\,\hbar\mu^{1.5}}{\pi^2 g^2\sqrt{m}} $$

- $\gamma \to 0$（弱耦合，平均场极限）:
  $$ n = \frac{\mu}{g} \quad\Rightarrow\quad y = x $$

- $\gamma \ll 1$（弱耦合领头阶+次领头阶）:
  $$ n = \frac{\mu}{g} + \frac{1}{\pi}\sqrt{\frac{m\mu}{\hbar^2}} $$

---

## 代码结构

- `BD_mission_ycr.ipynb`：主要的 Jupyter Notebook，包含：
  - 导入库（`numpy`, `matplotlib`）
  - 迭代函数 `bethe(γ, max_iter, ϵ)`：输入 $\gamma$ 数组，返回 $e(\gamma)$ 数组
  - 数值微分函数 `μ_dimless(γ, e)`：返回无量纲化学势
  - 数据生成：调用 `bethe` 得到 $e$，计算化学势，构建绘图数据 `x` 和 `y`
  - 渐近线数据生成（标记为 `yinfy`, `ybb1`, `y0`, `yss1`）
  - 绘图：使用 `scienceplots` 包的 `science` 风格，输出对数坐标对比图
- `README_bethe.txt`：简要任务描述

---

## 依赖环境

运行代码需要以下 Python 库：

- `numpy`
- `matplotlib`
- `scienceplots` （仅用于绘图风格，可通过 `pip install SciencePlots` 安装）

建议使用 Python 3.7 及以上版本。

---

## 快速开始

1. **克隆仓库**（或下载文件）：
   ```bash
   git clone https://github.com/your-username/your-repo-name.git
   cd your-repo-name
2. **安装依赖**：
   ```bash
   pip install numpy matplotlib SciencePlots
## 结果输出

- **e(γ) 曲线**: 线形图展示基态能量随相互作用强度的变化。

- **状态方程对比图（双对数坐标）**:
  - 黑色实线为 Bethe ansatz 数值结果；
  - 红色线为强耦合渐近（虚线：$\gamma \to \infty$；实线：$\gamma \gg 1$）；
  - 蓝色线为弱耦合渐近（虚线：$\gamma \to 0$；实线：$\gamma \ll 1$）。

数值结果与渐近表达式在小 $\gamma$ 和大 $\gamma$ 极限下高度一致，验证了算法的正确性。

## 参考资料

- **姚和朋, 郭彦良. 《冷原子物理与低维量子气体》. 北京: 科学出版社, 2024.**
- **杨文力，等. 《可积模型方法及其应用 》. 北京: 科学出版社, 2019.**
