# TPM-Vibration: Tooth Profile Modification & Gear Mesh Vibration Optimization

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](python/)
[![Python Tests](https://img.shields.io/badge/pytest-20%20passed-brightgreen.svg)](#running-the-tests)
[![License: GNU AGPLv3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE.md)
[![DOI](https://img.shields.io/badge/DOI-10.1007%2Fs40430--023--04574--3-informational.svg)](https://doi.org/10.1007/s40430-023-04574-3)

High-performance scientific Python package (`tpm`) for cylindrical spur gear tooth contact modeling, Time-Varying Mesh Stiffness (TVMS) calculation via the Energy Method, 1-DOF dynamic mesh vibration simulation, and optimal Tooth Profile Modification (TPM) tip relief determination.

---

## Academic & Scientific Context

This project was developed by **Flávio Dias** at the **Federal University of Amazonas (UFAM, Brazil)** as part of his Bachelor's thesis and subsequent scientific research on gear transmission dynamics. The primary goal is to minimize the vibration of cylindrical spur gear meshes in steady-state operation ($1 < \epsilon_\alpha < 2$) by determining the optimal **Tooth Profile Modification (TPM)** tip relief amplitude ($\Delta_{\max}$) for specified operating conditions.

### Publications

- **(EN) Journal Article (Springer Nature)**:  
  Dias, F.H.A.S., Silva, G.C. & Chui, D.S. *Vibration attenuation on spur gears through optimal profile modification based on an alternative dynamic model*. **Journal of the Brazilian Society of Mechanical Sciences and Engineering** 46, 10 (2024).  
  [https://doi.org/10.1007/s40430-023-04574-3](https://doi.org/10.1007/s40430-023-04574-3) | [ReadCube Free Access](https://rdcu.be/ds3ID)

- **(PT) Bachelor's Thesis (UFAM)**:  
  Dias, F.H.A.S. *Minimizing gear mesh vibration by the tooth profile modification of cylindrical spur gears*. 2022. 113f. Final Graduation Project (Mechanical Engineering) — Federal University of Amazonas, Manaus-AM, 2022.  
  [https://riu.ufam.edu.br/bitstream/prefix/6549/3/TCC_FlavioDias.pdf](https://riu.ufam.edu.br/bitstream/prefix/6549/3/TCC_FlavioDias.pdf)

---

## Features

- **Involute Tooth Geometry & Contact Relations**: Calculates standard spur gear geometry ($d_p, d_b, d_a, d_f$), tooth thicknesses, contact ratio ($1 < \epsilon_\alpha < 2$), and division of single/double contact zones.
- **Equivalent Mass & Inertia**: Formulates equivalent system mass ($m_e$) and explicitly accounts for the mass and inertia removed by tip relief modifications ($m_{emod}$).
- **Energy-Based TVMS**: Evaluates bending, shear, axial compression, fillet-foundation, and Hertzian contact compliances using vectorized Gauss-Legendre quadrature ($N=25$).
- **Time-Varying Friction & Damping**: Formulates contact friction and instantaneous damping ratios along the line of action based on Luo & Li's formulation and Thomson's dissipated energy approach.
- **Dynamic 1-DOF State-Space Solver**: High-precision numerical integration (LSODA / BDF) with continuous monotonic cubic Hermite splines (PCHIP) evaluated via Horner's rule.
- **Automated TPM Optimization**: Minimizes steady-state RMS acceleration via Brent's bounded scalar minimization.
- **Automatic Physics-Based Search Bounds**: Automatically bounds the search interval from single-tooth static contact deflection under load:
  $$\Delta_{ref} = \frac{F}{\min(K_{te})} \times 10^6\ [\mu\text{m}]$$
  $$\text{bounds} = [0.6 \times \Delta_{ref},\quad 1.4 \times \Delta_{ref}]$$
- **Comprehensive Validation Plots**: Generates comparative engineering figures for contact velocity, TVMS interpolation, time-domain mesh dynamics, normalized coordinates, and frequency spectra.

---

## Project Structure

```
TPM-Vibration/
├── README.md               # Main project documentation
├── LICENSE.md              # GNU Affero General Public License v3.0
├── plots/                  # Output directory for validation figures
├── python/                 # Python package and test suite
│   ├── README.md           # Python engine documentation
│   ├── requirements.txt    # Package dependencies
│   ├── pytest.ini          # Pytest configuration
│   ├── run_example.py      # End-to-end example script
│   ├── tpm/                # Core Python package
│   │   ├── __init__.py     # Public API exports
│   │   ├── geometry.py     # Involute tooth geometry & mesh relations
│   │   ├── mass.py         # Equivalent mass & polar moments of inertia
│   │   ├── energy.py       # TVMS (Energy Method), friction & damping
│   │   ├── dynamics.py     # 1-DOF state-space dynamic solver
│   │   ├── optimization.py # FastDynamicOptimizer & auto-bounds algorithm
│   │   ├── plots.py        # Engineering visualization routines
│   │   └── faithful.py     # GearSystem orchestrator
│   └── tests/              # Automated test suite (20 tests)
└── legacy-matlab/          # Original legacy MATLAB source code (reference archive)
```

---

## Quickstart

### 1. Installation

```bash
cd python
pip install -r requirements.txt
```

#### Dependencies
- Python >= 3.10
- NumPy >= 2.0
- SciPy >= 1.13
- Matplotlib >= 3.8
- pytest >= 8.0

### 2. Run the End-to-End Example

Execute the automated pipeline (geometry -> TVMS -> auto-bounded optimization -> validation plots):

```bash
cd python
python run_example.py
```

Console output:
```text
Reference deflection delta_ref_um: 43.881
Automatic bounds: (26.329, 61.433)
Optimized deltamax_um: 43.519
Objective (RMS a): 90.807 m/s²
Success: True
Function evaluations: 20
Plots written to: plots
- contact_relations: plots/contact_relations.png
- interpolation_results: plots/interpolation_results.png
- meshing_dynamics: plots/meshing_dynamics.png
- normalized_comparisons: plots/normalized_comparisons.png
- frequency_analysis: plots/frequency_analysis.png
```

### 3. Programmatic Usage

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
    rpm=2000.0,               # Pinion rotation (rpm)
    lubricant_density=870.0,  # Lubricant density (kg/m³)
    temperature_c=60.0,       # Operating temperature (°C)
)

# Run optimization with automatic physics-based bounds
result = optimize_deltamax(system, bounds=None, opt_tol=1e-4)
print(f"Optimal relief: {result['deltamax_um']:.3f} µm")

# Render comparison plots
generate_validation_plots(system, deltamax_um=result["deltamax_um"], output_dir="./plots")
```

---

## Running the Tests

To run the automated test suite (geometry, mass properties, compliance integrals, TVMS physics, dynamic integration, and optimization):

```bash
cd python
python -m pytest tests/ -v
```

---

## Output Plots

The pipeline exports five validation figures into `plots/`:

| Figure | Description |
| :--- | :--- |
| `contact_relations.png` | Relative velocity $v_{rel}$ along the contact line across meshing angle $\theta_1$. |
| `interpolation_results.png` | Continuous PCHIP time-varying mesh stiffness $K_{te}(t)$ and $K_{tem}(t)$. |
| `meshing_dynamics.png` | Time-domain comparison of dynamic transmission error (DTE), meshing acceleration, velocity, and dynamic mesh force (DMF). |
| `normalized_comparisons.png` | Dynamic variables mapped along the normalized contact coordinate $\Gamma \in [-0.5, 0.5]$. |
| `frequency_analysis.png` | Frequency spectrum of baseline vs. modified gear transmission. |

<p align="center">
  <img src="https://user-images.githubusercontent.com/44821460/230794069-064d7b4e-b72a-4928-a2f7-47bb70ca2b29.png" alt="Vibration over a tooth cycle" width="750" />
</p>

---

## Legacy MATLAB Reference (`legacy-matlab/`)

The original MATLAB source code is archived in [`legacy-matlab/`](legacy-matlab/) for historical reference, verification, and research reproducibility.

- **Optimization**: Run `optimization.m` in MATLAB. It invokes `fminbnd(@simulation, 30, 50, options)` and triggers `simulation2.m`.
- **Single Simulation**: Open `simulation.m`, uncomment `deltamax` (line 7), uncomment the `plots(...)` call (line 95), and run.

---

## Citation

If you use this model or software in your research or engineering work, please cite:

```bibtex
@article{Dias2024Vibration,
  author    = {Dias, Fl{\'a}vio H. A. S. and Silva, Geraldo C. and Chui, David S.},
  title     = {Vibration attenuation on spur gears through optimal profile modification based on an alternative dynamic model},
  journal   = {Journal of the Brazilian Society of Mechanical Sciences and Engineering},
  volume    = {46},
  number    = {1},
  pages     = {10},
  year      = {2024},
  doi       = {10.1007/s40430-023-04574-3}
}
```

```bibtex
@bachelorsthesis{Dias2022TCC,
  author    = {Dias, Fl{\'a}vio H. A. S.},
  title     = {Minimizing gear mesh vibration by the tooth profile modification of cylindrical spur gears},
  school    = {Federal University of Amazonas (UFAM)},
  address   = {Manaus, Brazil},
  year      = {2022},
  note      = {113 pages}
}
```

---

## License

This project is licensed under the **GNU Affero General Public License v3.0 (GNU AGPLv3)** — see the [LICENSE.md](LICENSE.md) file for details.

### Copyleft & SaaS Protection
Under the GNU AGPLv3, any software, library, derivative work, or network/cloud service (SaaS, Web APIs, simulation platforms) that incorporates or modifies this codebase must make its complete source code publicly available under the same GNU AGPLv3 license. Proprietary closed-source commercialization without returning source code contributions to the community is strictly prohibited.
