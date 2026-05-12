# Kagome 晶格 Fermi-Hubbard 模型自旋波基态 🌐[English Version](README_English.md)

![参数随 $U$ 变化图](figure_result.png)

本项目通过平均场近似与自旋波假设，解析导出 Kagome 晶格 Fermi-Hubbard 模型的能带和梯度，并在 $2/3$ 填充下使用梯度下降法迭代求解基态的磁化强度 $m_A,m_B,m_C$ 和自旋波矢 $q_{Ax},q_{Bx},q_{Cx},q_{Ay},q_{By},q_{Cy}$ ，得到序参量随相互作用 $U$ 的变化曲线。

---

## 物理背景

Kagome 晶格的原胞含有三个不等价原子 A、B、C。在平均场水平上，Fermi‑Hubbard 模型可通过引入磁化序 $m_i$ 和自旋波矢 $Q_i=(q_{ix},q_{iy})$ 进行拟设。无相互作用的三条能带 $\varepsilon_i(\mathbf{k})$ 由紧束缚给出，计入自旋波后形成六个准粒子能带：

$$
E_i^\pm(\mathbf{k}) = \frac{1}{2}\big[\varepsilon_i(\mathbf{k}) + \varepsilon_i(\mathbf{k}+Q_i)\big] \pm \sqrt{\frac{1}{4}\big[\varepsilon_i(\mathbf{k}) - \varepsilon_i(\mathbf{k}+Q_i)\big]^2 + (U m_i)^2},
$$

其中 $i=A,B,C$ ， $U$ 为相互作用强度。系统的总能量为所有占据态能量之和：

$$
E_{\text{total}} = \sum_{\mathrm{occ}} E_n(\mathbf{k}) + 2U\sum_{i}m_i^2 .
$$

在 $2/3$ 填充下， $6N_{\text{cell}}$ 个单粒子态中只有 $4N_{\text{cell}}$ 个被占据（ $N_{\text{cell}}$ 为原胞数）。利用手动推导的解析梯度，以梯度下降法同时优化总能量对 $m_i$ 和 $Q_i$ 的依赖，可以得到基态序参量。

本代码解决的核心问题：  
1. 解析推导六个能带对 $m_i,q_{ix},q_{iy}$ 的偏导数，并封装为数组运算；  
2. 在 $2/3$ 占据约束下，利用梯度下降寻找使总能量最小的参数组合，并研究 $U$ 驱动的相变行为。

---

## 数值方法

### 1. 解析导数和离散化

- **能带表达式**：  
  无相互作用能带 $\varepsilon_i(\mathbf{k})$（ $i=A,B,C$ ）在 $L\times L$ 的动量网格上计算。平移后的能带 $\varepsilon_i(\mathbf{k}+Q_i)$ 及其对 $Q_i$ 的导数 $\partial\varepsilon_i/\partial q_{ix},\partial\varepsilon_i/\partial q_{iy}$ 也同时给出。

- **六带能谱与梯度**：  
  $E_i^\pm(\mathbf{k})$ 及其对 $m_i,q_{ix},q_{iy}$ 的偏导数被解析地写为 `dEdm`、`dEdqx`、`dEdqy` 三个数组，形状均为 $(6,L,L)$ ，方便向量化运算。

### 2. 梯度下降求解基态

- **初值**：  
  对每个 $U$ ，用 $m_i = 0.1U, q_{ix}=q_{iy}=0.8\pi U/3$ 加上小随机噪声初始化，利用前一个 $U$ 的结果连续扫描时不需要特殊衔接。

- **迭代更新**：  
  1. 计算当前参数下的六条能带 $E$ 和梯度。  
  2. 将所有 $6N_{\text{cell}}$ 个能量值展平，按升序排列，仅取前 $4N_{\text{cell}}$ 个最低能级作为占据态。  
  3. 对占据态上的梯度求和得到总能量对 $m_i$ 和 $Q_i$ 的梯度，并额外计入 $2U m_i$ 的贡献。  
  4. 按固定学习率 $\alpha=0.1$ 更新参数：  

$$
m_i \leftarrow m_i - \alpha \frac{\partial E_{\text{total}}}{\partial m_i}, \quad q_{i\nu} \leftarrow q_{i\nu} - \alpha \frac{\partial E_{\text{total}}}{\partial q_{i\nu}}.
$$
  
- **收敛判据**：  
  当三个分量的最大更新步长均小于 $10^{-4}$ 时停止迭代。最大迭代次数未硬性限制，实际通常几十到几百步收敛。

### 3. 扫描和结果处理

对一系列 $U$ 值（从 $0.2$ 到 $6.4$ ，步长 $0.1$ ）重复上述优化过程，记录收敛后的 $m_i,q_{ix},q_{iy}$ 。最后绘制序参量随 $U$ 的变化，并标识已知的相变点 $U_{c1}\approx2.6$ 和 $U_{c2}\approx5.35$。

---

## 代码结构

- 主程序文件（单脚本或 Jupyter Notebook）：  
  - 参数设置：晶格大小 $L$ ， $U$ 列表，学习率 $\alpha$ ；  
  - 动量空间基矢 $\mathbf{b}_1,\mathbf{b}_2$ 及网格生成；  
  - 无相互作用能带 $\varepsilon$ 的计算（三带）；  
  - 六带能谱 $E$ 和梯度 `dEdm`,`dEdqx`,`dEdqy` 的解析表达式；  
  - 梯度下降主循环：能量排序、占据选择、梯度累积、参数更新；  
  - 收纳结果数组 `m_g`, `qx_g`, `qy_g` ；  
  - 绘图和 `np.savez` 保存数据。

---

## 依赖环境

运行代码需要以下 Python 库：

- `numpy`
- `matplotlib`

建议使用 Python 3.7 及以上版本。

---

## 快速开始

1. **克隆仓库**（或下载文件）：
   ```bash
   git clone https://github.com/chaoranyang/QuantumManyBody_FromZreo.git
   cd QuantumManyBody_FromZreo/kagomeFH_Grad
   
2. **安装依赖**：  
   ```bash
   pip install numpy matplotlib
---

## 结果输出

- **序参量随 $U$ 的变化图**：  
  三幅子图分别展示 $m_i$, $q_{ix}$, $q_{iy}$ 随相互作用强度 $U$ 的演化。
  
  图一中可看到磁化强度在临界 $U_{c1}\approx2.6$ 处开始自发产生；  
  图二显示 $q_{ix}$ 在 $U_{c1}$ 处跳变到约 $2.199$ 的平台，并在 $U_{c2}\approx5.35$ 再次变化；  
  图三给出 $q_{iy}$ 的行为，与系统对称性一致。

  数值结果清晰复现了 Kagome Hubbard 模型中自旋密度波态的相变特征。
---

## 参考资料
- **Patrik Fazekas**, *Lecture Notes on Electron Correlation and Magnetism*, in *Series in Modern Condensed Matter Physics*, Vol. 5, World Scientific, 1999, pp. 1-796.
