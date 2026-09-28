# TPM-Vibration (Python Engine)

High-performance scientific Python engine for cylindrical spur gear modeling, time-varying mesh stiffness (TVMS), nonlinear elastohydrodynamic (EHL) damping, dynamic mesh response, and tooth profile modification (TPM) tip relief optimization.

---

## Key Features

- **100% Mathematical Parity**: Tooth geometry (`gears.m`, `mesh.m`), equivalent mass (`mass.m`), TVMS energy integrals (`energy.m`, `integr.m`), and dynamic equations of motion (`dynamic.m`, `dynamictpm.m`, `model.m`) match MATLAB to double-precision tolerances ($< 0.01\%$).
- **Vectorized Gauss-Legendre Quadrature ($N=25$)**: Replaces slow point-by-point QUADPACK adaptive loops, matching machine precision ($2.5 \times 10^{-15}$) in milliseconds.
- **Precomputed Base State Optimization**: Decouples invariant quantities ($K_{te}$, $K_{12}$, $C$, $\mu$) from the profile relief scaling $\delta(t) = \Delta_{\max} \cdot \bar{\delta}(t)$ and equivalent mass $m_{emod}(\Delta_{\max})$.
- **Automatic Physical Bounds**: Computes the search interval automatically from single-tooth static contact deflection under load:
  $$\Delta_{ref} = \frac{F}{\min(K_{te})} \times 10^6\ [\mu\text{m}]$$
  $$\text{bounds} = [0.6 \times \Delta_{ref},\quad 1.4 \times \Delta_{ref}]$$
- **> 1,000x Speedup**: Reduces optimization run time from **~30 minutes** in MATLAB to **~1.8 seconds** in Python.

---

## Installation

```bash
cd python
pip install -r requirements.txt
```

### Dependencies
- Python 3.10+
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
    z1=27,                  # Pinion teeth
    z2=35,                  # Gear teeth
    module=3.0,             # Module (mm)
    alpha_deg=20.0,         # Pressure angle (deg)
    b_mm=25.0,              # Face width (mm)
    elasticity_modulus=206e9, # Modulus of elasticity (Pa)
    poisson_ratio=0.3,      # Poisson ratio
    density=7850.0,         # Density (kg/m³)
    power_w=80e3,           # Pinion power (W)
    rpm=2000.0,             # Pinion speed (rpm)
    lubricant_density=870.0,# Lubricant density (kg/m³)
    temperature_c=60.0,     # Operating temperature (°C)
)

# Run optimization with automatic physics-based bounds
result = optimize_deltamax(system, bounds=None, opt_tol=1e-4)
print(f"Optimal relief: {result['deltamax_um']:.3f} µm")
print(f"Minimum RMS acceleration: {result['objective']:.4f} m/s²")

# Generate comparative validation plots
generate_validation_plots(system, deltamax_um=result["deltamax_um"], output_dir="./plots")
```

---

## Running Tests

Run the full test suite (unit tests, mathematical parity, dynamic response, optimization):

```bash
cd python
python -m pytest tests/ -v
```
