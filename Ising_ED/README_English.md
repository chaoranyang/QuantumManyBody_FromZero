# Exact Diagonalization of the Transverse Field Ising Model and Verification of Physical Properties

This project numerically solves the one-dimensional transverse field Ising model using exact diagonalization. It computes the energy spectrum, ground-state magnetization, and its variance for different system sizes. Through finite-size scaling analysis, the critical exponents of the quantum phase transition are obtained and compared with the theoretical universality class.

---

## Physical Background

The Hamiltonian of the one-dimensional transverse field Ising model is

$$
H = -J \sum_{i=1}^N \sigma_i^z \sigma_{i+1}^z - h \sum_{i=1}^N \sigma_i^x
$$

where $J$ is the nearest-neighbor spin exchange interaction (set to $J=1$), $h$ is the transverse magnetic field, and $\sigma_i^{z,x}$ are Pauli matrices. The model undergoes a continuous quantum phase transition from the ferromagnetic phase to the paramagnetic phase at $h = 1$, belonging to the universality class of the two-dimensional classical Ising model with critical exponents $\beta = 1/8$ and $\nu = 1$.

Core problems addressed by this code:  
1. Using binary representation of spin configurations, construct the Hamiltonian matrix for arbitrary system size $N$, obtain eigenvalues and eigenvectors via exact diagonalization, and plot the energy spectrum as a function of the transverse field $h$.  
2. Compute the ground-state average absolute magnetization $\langle|m_z|\rangle$ and its variance $\sigma_{m_z}$ for different $N$, and observe the quantum phase transition features and finite-size effects near $h=1$.  
3. Perform a scaling transformation on the magnetization data to achieve data collapse, and use high-precision grid search to fit the scaling function, finding optimal critical exponents $\beta$ and $\nu$ for comparison with the Ising universality class.

---

## Numerical Methods

### 1. Exact Diagonalization for Energy Spectrum

![Result Figure](figure_result_1.png)

- **Spin configuration representation**: Represent the $2^N$ configurations of $N$ spins as $N$-bit binary numbers, where $0$ corresponds to spin down and $1$ to spin up.
- **Hamiltonian construction**:
  - Diagonal elements: Compute the contribution of $-J \sum\sigma_i^z \sigma_{i+1}^z$ for each configuration under periodic boundary conditions.
  - Off-diagonal elements: The transverse field term $-h \sum \sigma_i^x$ flips a single spin, adding matrix element $-h$ between the two corresponding states.
- **Diagonalization**: Use `numpy.linalg.eigvalsh` to compute eigenvalues of the Hermitian matrix, taking the lowest few energy levels (default 50).
- **Parameter sweep**: Vary the transverse field $h$ from $[0, 2.0]$ with step size $0.05$, recording the energy spectrum at each $h$.

### 2. Magnetization and Variance Calculation

![Result Figure](figure_result_2.png)

- **Ground-state wavefunction**: For each $h$, use `numpy.linalg.eigh` to obtain both eigenvalues and eigenvectors, selecting the ground state (lowest energy eigenstate).
- **Expectation value of magnetization**:
  
$$
\langle|m_z|\rangle = \sum_{\text{states}} |c_i|^2 \cdot \frac{|N_{\uparrow} - N_{\downarrow}|}{N}
$$
  
  where $|c_i|^2$ is the probability of each configuration in the ground state, and $N_{\uparrow}, N_{\downarrow}$ are the numbers of up and down spins in that configuration.
- **Variance of magnetization**:

$$
\sigma_{m_z}^2 = \langle m_z^2 \rangle - \langle|m_z|\rangle^2
$$

- **System sizes**: Compute for $N = 4, 6, 8, 10, 12$ respectively.

### 3. Finite-Size Scaling and Data Collapse

![Result Figure](figure_result_3.png)

- **Scaling hypothesis**: Near the critical point, the magnetization satisfies the scaling relation

$$
|m(h, L)| = L^{-\beta/\nu}  \mathcal{M}\left( (h/h_c - 1) L^{1/\nu} \right)
$$

  where $L$ is the system size ($L=N$), $h_c=1$ is the critical field, and $\mathcal{M}$ is a universal scaling function.
- **Data collapse**: Define scaling variables

$$
X = (h/h_c - 1) L^{1/\nu}, \quad Y = |m| L^{\beta/\nu}
$$

  If the parameters $\beta,\nu$ are appropriate, data for different $L$ will collapse onto a single curve in the $(X,Y)$ plane.
- **High-precision grid search**:
  - Generate a $101\times101$ grid over $\beta \in [0.08, 0.18]$, $\nu \in [0.8, 1.2]$.
  - For each $(\beta,\nu)$ pair, compute the collapsed $X, Y$, select data with $|X|<2$ (critical region), expand $\mathcal{M}$ as a quadratic polynomial $Y = C_0 + C_1 X + C_2 X^2$ and fit, then calculate the coefficient of determination $R^2$.
  - The $(\beta,\nu)$ giving the maximum $R^2$ is taken as the optimal critical exponents.

---

## Code Structure

- `Ising_ED.ipynb`: Main Jupyter Notebook containing three parts:
  - **Part 1**: Energy spectrum calculation.
    - `build_hamiltonian(N, h, J=1.0)` constructs the Hamiltonian,
    - `cal_energy_level(N, h_values, num_states=50)` computes energy levels,
    - `plot_energy_levels(...)` plots the energy spectrum.
  - **Part 2**: Magnetization and variance calculation.
    - `cal_mag(N, state_bin)` computes magnetization for a single state,
    - `cal_mag_data(N_list, h_values)` computes magnetization and variance for all system sizes,
    - `plot_mag_and_var(...)` plots the results.
  - **Part 3**: Data collapse and optimal parameter search.
    - `grid_search_optimal_params_high_precision(mag_data, h_values, h_c=1.0)` performs high-precision grid search,
    - `plot_all_data_with_fit(...)` plots the collapsed data points and the fitted curve.

---

## Dependencies

Running the code requires the following Python libraries:

- `numpy`
- `matplotlib`

Python 3.7 or higher is recommended.

---

## Quick Start
  1. **Clone the repository** (or download the files):
     ```bash
     git clone https://github.com/chaoranyang/QuantumManyBody_FromZero.git
     cd QuantumManyBody_FromZero/Ising_ED
 2. **Install dependencies:**
    ```bash
    pip install numpy matplotlib
---

## Results Output

- **Energy spectrum plot**: Shows the lowest 50 energy levels (relative to ground state energy) as a function of transverse field $h$ for $N=12$.
- **Magnetization and variance plot**: Left panel shows magnetization for different $N$, right panel shows the corresponding variance. Phase transition features and finite-size effects are observed near $h=1$.
- **Data collapse plot**: Using the optimal $\beta,\nu$, magnetization data from different sizes are transformed onto the same scaling curve, with quadratic fit and $R^2$ value displayed.
- **Parameter search results**: Outputs the optimal $\beta$, $\nu$, and fit coefficients $C_0, C_1, C_2$, along with relative errors compared to the theoretical values ($\beta=0.125$, $\nu=1.000$).

## References

- **A. W. Sandvik, A. Avella, and F. Mancini, "Computational Studies of Quantum Spin Systems," in *AIP Conference Proceedings*, 2010, pp. 135–338. DOI: [10.1063/1.3518900](http://dx.doi.org/10.1063/1.3518900)**
