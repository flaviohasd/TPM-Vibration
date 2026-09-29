# TPM-Vibration (Python Engine)

High-performance scientific Python engine for cylindrical spur gear tooth contact modeling, Time-Varying Mesh Stiffness (TVMS) via the Energy Method, friction-based energy dissipation mesh damping, 1-DOF dynamic mesh vibration simulation, and optimal Tooth Profile Modification (TPM) tip relief determination.

---

## Key Features

- **100% Mathematical Parity**: Tooth geometry (`gears.m`, `mesh.m`), equivalent mass and inertia (`mass.m`), TVMS energy integrals (`energy.m`, `integr.m`), and dynamic equations of motion (`dynamic.m`, `dynamictpm.m`, `model.m`) match the original MATLAB reference code to double-precision tolerances ($< 0.01\%$).
- **Energy-Based TVMS**: Evaluates tooth bending, shear, axial compression, Cai-Sainsot fillet-foundation, and Hertzian contact compliances in series along the line of action ($1 < \epsilon_\alpha < 2$).
- **Vectorized Gauss-Legendre Quadrature ($N=25$)**: Replaces slow point-by-point QUADPACK adaptive loops, matching machine precision ($2.5 \times 10^{-15}$) in milliseconds.
- **Frictional Energy Dissipation Mesh Damping**: Formulates instantaneous mesh damping $C(t) = 2 \zeta(t) \sqrt{m_e K_{te}(t)}$ based on frictional dissipation work ($W_a$) and contact strain energy ($U$), yielding the time-varying damping ratio $\zeta = \frac{W_a}{4\pi U}$.
- **Precomputed Base State Optimization**: Decouples invariant quantities ($K_{te}$, $K_{12}$, $C$, $\mu$) from the profile relief scaling $\delta(t) = \Delta_{\max} \cdot \bar{\delta}(t)$ and equivalent mass $m_{emod}(\Delta_{\max})$, evaluating ODE states with vectorized cubic Hermite splines (PCHIP) via Horner's rule.
- **Automatic Physical Bounds**: Automatically bounds the search interval from single-tooth static contact deflection under load:
  $$\Delta_{ref} = \frac{F}{\min(K_{te})} \times 10^6\ [\mu\text{m}]$$
  $$\text{bounds} = [0.6 \times \Delta_{ref},\quad 1.4 \times \Delta_{ref}]$$
- **Automated TPM Optimization**: Minimizes steady-state RMS acceleration via **Brent's bounded minimization method with Golden Section Search** (1:1 mathematical equivalent to MATLAB's `fminbnd`).
- **Parametric Resolution Presets**: Configurable meshing sweep density (`fast` = 50, `standard` = 100, `fine` = 200, `ultra` / `publication` = 400 points/tooth, up to 10,800 steps/rev) or custom `points_per_tooth`.
- **> 100x to > 1,000x Speedup**: Reduces optimization run time from **~30 minutes** in MATLAB to **~1.8–15 seconds** in Python.

---

## Package Architecture

The `tpm` package is organized into modular scientific components:

```
tpm/
├── geometry.py      # Involute tooth geometry & line of action contact relations
├── mass.py          # Equivalent mass (me) and tip-relief modified mass (memod)
├── energy.py        # TVMS compliances (Gauss-Legendre), friction & energy dissipation damping
├── dynamics.py      # 1-DOF dynamic state-space solver (LSODA / BDF)
├── optimization.py  # FastDynamicOptimizer & physics-based auto-bounds algorithm
├── plots.py         # Comparative engineering visualization routines
├── faithful.py      # High-level GearSystem orchestrator
└── friction.py      # Friction coefficient modeling
```

---

## Installation

```bash
cd python
pip install -r requirements.txt
```

### Dependencies
- Python >= 3.10
- NumPy >= 2.0
- SciPy >= 1.13
- Matplotlib >= 3.8
- pytest >= 8.0

---

## Quickstart

Run the full end-to-end example (geometry -> base TVMS -> auto-bounded optimization -> validation plots):

```bash
cd python
python run_example.py
```

### Python API Usage

```python
from tpm import GearSystem, optimize_deltamax, generate_validation_plots

# Define gear pair operating conditions
system = GearSystem(
    z1=27,                    # Pinion teeth
    z2=35,                    # Gear teeth
    module=3.0,               # Module (mm)
    alpha_deg=20.0,           # Pressure angle (deg)
    b_mm=25.0,                # Face width (mm)
    elasticity_modulus=206e9, # Modulus of elasticity (Pa)
    poisson_ratio=0.3,        # Poisson's ratio
    density=7850.0,           # Density (kg/m³)
    power_w=80e3,             # Pinion power (W)
    rpm=2000.0,               # Pinion speed (rpm)
    lubricant_density=870.0,  # Lubricant density (kg/m³)
    temperature_c=60.0,       # Operating temperature (°C)
)

# Run optimization with automatic physics-based bounds (resolution='standard' | 'ultra')
result = optimize_deltamax(system, bounds=None, opt_tol=1e-4, resolution="standard")
print(f"Optimal relief: {result['deltamax_um']:.3f} µm")
print(f"Minimum RMS acceleration: {result['objective']:.4f} m/s²")

# Generate comparative validation plots
generate_validation_plots(system, deltamax_um=result["deltamax_um"], output_dir="./plots", resolution="standard")
```

---

## Running Tests

Run the full automated test suite (unit tests, mathematical parity, dynamic response, optimization, resolution presets):

```bash
cd python
pytest tests/ -v
```
