# Surface Polariton Solver (Python)

A Python-based numerical solver for calculating surface polariton dispersion relations using the Burke, Stegeman, and Tamir (1986) model.

## 📖 Overview

This script computes the complex propagation constant $\beta$ for surface polaritons in a three-layer waveguide system ($\varepsilon_1$ / Metal / $\varepsilon_3$). It uses **Muller's Method** for robust root-finding in the complex plane and employs **$h$-continuation** to ensure accurate mode tracking.

## 🛠️ Prerequisites

Ensure you have Python 3.8+ installed along with the following libraries:
- `numpy`
- `matplotlib`

You can install them via pip:
```bash
pip install numpy matplotlib
```

## 🚀 Usage

Run the solver directly from the terminal. It will prompt you for physical parameters:

```bash
python surface_polariton_solver.py
```

### Interactive Inputs
- **$\varepsilon_1$ & $\varepsilon_3$**: Relative permittivity of the dielectric cladding and substrate.
- **Metal Mode**: Choose between (D)rude model or (F)ixed permittivity.
- **Thickness ($h$)**: Metal film thickness in nanometers (nm).
- **$\omega_p$**: Plasma frequency in rad/s.
- **$\Gamma$**: Drude damping rate in rad/s.

## 📊 Workflow & Outputs

Upon successful completion, the solver generates the following files in the `output_py_solver/` directory:
1. **CSV Data**: Raw numerical results in physical SI units.
    - Columns: $\omega$ [rad/s], $k_0$ [rad/m], $\beta_R$ [rad/m], $\beta_I$ [rad/m], $\varepsilon_m$.
2. **PNG Plot**: A two-panel publication-ready visualization of the dispersion curves (Real and Imaginary parts).
3. **HTML Report**: An interactive summary page containing the plot and a table of the computed values.

## 🔬 Core Algorithms
- **Muller's Method**: Iterative complex root-finder for implicit dispersion relations.
- **$h$-continuation**: Robustly tracks mode branches from the thick-film limit down to the physical thickness.
- **Drude Model**: Dynamically calculates metal permittivity $\varepsilon_m(\omega) = 1 - \frac{\omega_p^2}{\omega(\omega + i\Gamma)}$.

---
*Created as part of the Surface Polariton Implementation project at UEF.*
