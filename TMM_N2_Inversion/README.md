# Fabry-Perot N2 Inversion

This repository contains Python scripts for extracting the complex refractive index ($n_2 = n_{real} + i \cdot n_{imag}$) of a dielectric slab from measured complex transmission coefficients, and for visualizing/comparing those results.

## Prerequisites

The scripts require Python 3 and the following packages:
```bash
pip install numpy scipy matplotlib
```

---

## 1. Extracting Refractive Index (`extract_n2.py`)

This script numerically inverts the exact Fabry-Perot transmission formula at normal incidence across a broad frequency spectrum. Given a time-domain or frequency-domain measured transmission $|t|^2$ and phase, it uses a bounded non-linear root solver to extract the complex refractive index $n_2(\omega)$.

### Implementation Details:
The script is designed for extreme numerical stability and to prevent the common issue of solver branch-jumping when evaluating highly transparent, highly absorbing, or thick dielectric samples. It features three main components:

1. **Numerically Stable Forward Model:**
   The exact calculation $t = \frac{t_{12}t_{23}e^{i\phi}}{1 - r_{21}r_{23}e^{2i\phi}}$ diverges computationally for strongly absorbing materials ($n_{imag} \gg 0$, $\text{Im}(\phi) > 0$). When the script detects this regime, it actively factorizes the exponent $e^{i\phi}$ out of the denominator to bound all evaluated exponentials to magnitudes $\le 1$.
2. **Robust Root Solving & Guesses:**
   The script iterates frequency by frequency, mapping the experimental data against the forward model using `scipy.optimize.root`. It tries both Powell's hybrid (`hybr`) and Levenberg-Marquardt (`lm`) algorithms to hit a conservative residual tolerance of $< 10^{-10}$. To ensure continuity without forcing the solver over unphysical branch cuts, it explicitly tries up to 10 initial heuristic guesses (`build_guesses`) at every frequency, seeded heavily by the *previous* point's successful convergence.
3. **Outlier Cleanup:**
   If the solver hits an isolated root branch jump (e.g. flipping instantaneously between a real-dominated mode and an imaginary-dominated loss mode), `cleanup_outliers` will detect standard mathematical discontinuity. It sweeps broad geometric interpolations from surrounding correctly converged indices to force the root solver back onto the true physical branch, and then resweeps forward recursively.

### Input Data Format (CSV - no header)
The script expects a 3-column CSV file:
`Frequency (THz), |t|^2, phase (radians)`

### Output Data Format (CSV - with header)
The script outputs a 3-column CSV file:
`Frequency (THz), n_real, n_imag`

### Usage
Run the script with the input CSV and the slab thickness (in meters).
```bash
# General usage: python extract_n2.py <input.csv> <thickness_m> [-o <output.csv>]
python extract_n2.py modelled_trnsm.csv 5e-7
```

**Custom Outputs:**
```bash
python extract_n2.py modelled_trnsm.csv 5e-7 -o my_custom_output.csv
```

*Tip: You can modify the constants directly at the top of `extract_n2.py` (such as `INPUT_FILE`, `SLAB_THICKNESS`, `N1`, `N3`) to allow running the script without any command-line arguments.*

---

## 2. Comparing Results (`compare_n2.py`)

This script generates a visually clean 2x2 grid comparing two inverted `n2` CSV files. It overlays their real and imaginary parts side-by-side, and directly computes the absolute numerical differences ($\Delta n_{real}$, $\Delta n_{imag}$).

### Usage
Pass the two CSV results (must have headers) that you want to compare as positional arguments.
```bash
# General usage: python compare_n2.py <output1.csv> <output2.csv>
python compare_n2.py output_n2.csv result_refractive_index.csv
```

A window will pop up showing the plotted overlays and error residuals automatically.
