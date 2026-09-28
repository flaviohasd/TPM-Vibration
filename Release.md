# Release v3.0.0: Native Python Engine (`tpm`), >1,000x Faster Optimization & GNU AGPLv3

This major release marks the complete evolution of the **TPM-Vibration** simulator from legacy MATLAB scripts into a modern, high-performance scientific Python package (`tpm`). It features exact mathematical parity, automated physics-based optimization bounds, a >1,000x optimization speedup, and a license transition to GNU AGPLv3 for enhanced copyleft and SaaS protection.

---

## 🚀 Highlights & What's New

### 1. Native Python Engine (`tpm`)
- **Complete Object-Oriented Interface**: Introduced `GearSystem` to orchestrate tooth geometry, contact mechanics, and mass properties seamlessly.
- **Tooth Geometry (`geometry.py`)**: 100% exact algebraic parity with the legacy MATLAB implementation (`gears.m` and `mesh.m`) with $0.0000000000\%$ error.
- **Equivalent Mass & Inertia (`mass.py`)**: Computes equivalent mass and polar moments of inertia ($m_e$ and $m_{emod}$) accounting for addendum relief. Includes a parameter `correct_centroid_typo: bool = False` to provide both legacy MATLAB bug-for-bug compatibility and corrected physical formulations.
- **Vectorized Energy Method TVMS (`energy.py`)**: Replaced slow point-by-point QUADPACK adaptive loops with 25-point Gauss-Legendre quadrature, evaluating tooth compliance integrals ($I_{kbi}, I_{ksi}, I_{kai}$) with machine precision ($2.5 \times 10^{-15}$) in ~11 ms.
- **Dynamic 1-DOF State-Space Solver (`dynamics.py`)**: Stiff ODE integration (LSODA / BDF) evaluated with continuous monotonic cubic Hermite splines (PCHIP) via Horner's rule, ensuring continuous forces and accurate steady-state vibration metrics.

### 2. Algorithmic Breakthrough: >1,000x Faster Optimization
- **Base State Decoupling**: In the legacy implementation, MATLAB re-ran the full geometry, mesh, and TVMS pipeline on every scalar evaluation of `fminbnd`, taking ~30 minutes.
- **`FastDynamicOptimizer` (`optimization.py`)**: Base TVMS $K_{te}(t)$, damping $C(t)$, and Horner spline tables are precomputed **once**. Each iteration only scales the linear profile relief $\delta(t) = \Delta_{\max} \cdot \bar{\delta}(t)$ and updates $m_{emod}$, executing each dynamic simulation in $<80\text{ ms}$.
- **Result**: Complete optimization runs in **1.77 seconds** (converging in 20 function evaluations) instead of ~30 minutes, finding $\Delta_{\max} = 43.518\ \mu\text{m}$ (a 0.14% match to the published value of $43.581\ \mu\text{m}$).

### 3. Automatic Physics-Based Search Bounds
- Removed the need for manual, trial-and-error search bounds (e.g., `[30, 50]` $\mu\text{m}$).
- Automatically derives the search interval from the physical single-tooth static contact deflection under normal operating load:
  $$\Delta_{ref} = \frac{F}{\min(K_{te})} \times 10^6\ [\mu\text{m}]$$
  $$\text{bounds} = [0.6 \times \Delta_{ref},\quad 1.4 \times \Delta_{ref}]$$
- For the reference system, $\Delta_{ref} = 43.88\ \mu\text{m}$, creating automatic bounds of $[26.33, 61.43]\ \mu\text{m}$ that robustly bracket the global optimum for any spur gear geometry and load rating.

### 4. Automated Test Suite (20 Tests)
- Integrated comprehensive `pytest` regression suite in `python/tests/` covering:
  - Involute gear tooth geometry and contact line lengths
  - Equivalent mass and inertia calculations
  - Vectorized compliance integrals and TVMS energy conservation
  - Dynamic solver integration tolerances
  - Automatic bounds calculation and optimization convergence

### 5. Repository Reorganization & Legacy Archive
- The original MATLAB scripts have been moved to [`legacy-matlab/`](legacy-matlab/) intact for historical preservation, reproducibility, and ground-truth validation.
- Renamed the core package from `tpmtpm` to `tpm`.

### 6. License Upgrade: MIT → GNU AGPLv3
- Upgraded license to the **GNU Affero General Public License v3.0 (GNU AGPLv3)**.
- **Copyleft & SaaS Protection**: Requires that any software, derivative work, or network service (SaaS, cloud platform, web API) utilizing or modifying this codebase must release its full source code under the same GNU AGPLv3 license, preventing proprietary commercial exploitation without returning improvements to the open-source community.

---

## 📊 Benchmark Summary

| Metric | MATLAB Reference (`legacy-matlab/`) | Python Engine (`tpm`) | Agreement |
| :--- | :--- | :--- | :--- |
| **Geometry** ($d_p, d_b, r_a, \epsilon_\alpha$) | `gears.m` / `mesh.m` | `geometry.py` | Exact ($0.0000000000\%$) |
| **Equivalent Mass** ($m_e$) | $0.349311\text{ kg}$ | $0.349311\text{ kg}$ | Exact ($0.0000000000\%$) |
| **Mean TVMS** ($\bar{K}_{te}$) | $3.54 \times 10^8\text{ N/m}$ | $3.5405 \times 10^8\text{ N/m}$ | $< 0.01\%$ |
| **Search Bounds** | `[30, 50] µm` (manual) | `[26.33, 61.43] µm` (auto) | Automatic physics formula |
| **Optimal Relief ($\Delta_{\max}$)** | **$43.581\ \mu\text{m}$** | **$43.518\ \mu\text{m}$** | **$0.14\%$ difference** |
| **Optimization Runtime** | **~30 minutes** | **1.77 seconds** | **> 1,000x faster** |

---

## 📦 Quickstart

```bash
cd python
pip install -r requirements.txt
python run_example.py