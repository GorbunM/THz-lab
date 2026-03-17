# Fabry-Perot N2 Inversion

Extract the complex refractive index ($n_2 = n_\text{real} + i \cdot n_\text{imag}$) of a dielectric slab from measured complex transmission coefficients, and visualize/compare those results.

## Prerequisites

Python 3 with:
```bash
pip install numpy scipy matplotlib
```

---

## 1. Extracting Refractive Index (`extract_n2.py`)

Numerically inverts the exact Fabry-Perot transmission formula at normal incidence:

$$t = \frac{t_{12}\,t_{23}\,e^{-i\phi}}{1 - r_{21}\,r_{23}\,e^{-2i\phi}}, \qquad \phi = \frac{2\pi\,n_2\,d\,f}{c}$$

for a three-layer system $n_1 \,|\, \text{slab}(n_2, d) \,|\, n_3$ at normal incidence.

### Usage

```bash
python extract_n2.py <input.csv> <thickness_m> [-p arg|delay] [-o output.csv]
```

### Command-Line Arguments

| Argument | Required | Description |
|---|---|---|
| `input.csv` | Yes* | Input CSV file (3 columns: frequency, \|t\|^2, phase) |
| `thickness_m` | Yes* | Slab thickness in metres (e.g. `80e-6` for 80 um) |
| `-p`, `--phase-convention` | No | Phase convention: `arg` or `delay` (default: `delay`) |
| `-o`, `--output` | No | Output CSV path (default: `output_n2.csv`) |

*Both can be omitted if `INPUT_FILE` and `SLAB_THICKNESS` are set at the top of the script.

### The `-p` / `--phase-convention` Flag

This flag tells the script how to interpret column 3 (phase) in the input CSV. **Choosing the wrong convention produces completely wrong results** (see `BUGFIX.md` for details).

| Flag value | Column 3 contains | Typical sign | When to use |
|---|---|---|---|
| `delay` (default) | Optical delay phase $= -\arg(t)$ | Positive, increasing | THz-TDS data, commercial THz systems |
| `arg` | $\arg(t)$, the complex argument of $t$ | Negative or wrapped | Simulation outputs, analytical calculations |

**Quick rule:** if your phase column is mostly positive and increases with frequency, use `-p delay`. If it is mostly negative or wraps around zero, use `-p arg`.

### Examples

```bash
# THz-TDS data with optical-delay phase (positive, increasing):
python extract_n2.py jin_trnsm_with_phase.csv 80e-6 -p delay

# Modelled/simulated data where phase column is arg(t):
python extract_n2.py modelled_trnsm.csv 5e-7 -p arg

# Custom output path:
python extract_n2.py data.csv 50e-6 -p delay -o result.csv

# Using built-in defaults (edit the top of extract_n2.py):
python extract_n2.py
```

### Input Data Format (CSV, no header)

Three columns, comma-separated:

| Column | Description | Units |
|--------|-------------|-------|
| 1 | Frequency | THz |
| 2 | Power transmission $\|t\|^2$ | dimensionless |
| 3 | Phase | radians |

The interpretation of column 3 depends on the **phase convention** (`-p` flag) — see above.

### Output Data Format (CSV, with header)

| Column | Description | Convention |
|--------|-------------|------------|
| Frequency (THz) | Same as input | |
| n_real | Real part of refractive index | $\ge 0$ for physical materials |
| n_imag | Imaginary part (extinction) | $\ge 0$ for absorbing materials |

The output always uses the standard physics convention where $n_\text{imag} > 0$ indicates absorption, regardless of which phase convention was used for the input.

### Configurable Defaults

Edit the constants at the top of `extract_n2.py` to run without CLI arguments:

| Variable | Purpose | Example |
|---|---|---|
| `INPUT_FILE` | Default input CSV path | `"jin_trnsm_with_phase.csv"` |
| `OUTPUT_FILE` | Default output CSV path | `"output_n2.csv"` |
| `SLAB_THICKNESS` | Slab thickness in metres | `80e-6` (80 um) |
| `PHASE_CONVENTION` | Phase convention | `"delay"` or `"arg"` |
| `N1`, `N3` | Surrounding media refractive indices | `1.0` (air) |
| `RESIDUAL_TOL` | Convergence threshold | `1e-10` |

### Implementation Details

1. **Numerically Stable Forward Model:**
   The exact Fabry-Perot formula diverges computationally for strongly absorbing materials ($\text{Im}(\phi) > 0$). When this regime is detected, the exponential is factored to bound all evaluated terms to magnitudes $\le 1$.

2. **Robust Root Solving:**
   The script iterates frequency-by-frequency using `scipy.optimize.root` with both Powell's hybrid (`hybr`) and Levenberg-Marquardt (`lm`) methods. Up to 10 initial guesses per frequency point are seeded from the previous converged solution and heuristic estimates, targeting a residual $< 10^{-10}$.

3. **Branch Selection:**
   The Fabry-Perot formula satisfies $t(n_2) = t(-n_2)$ when $n_1 = n_3$, creating a mirror-image ambiguity: both $(a, b)$ and $(-a, -b)$ are valid roots. The solver tracks the $n_\text{imag} \ge 0$ branch internally because this provides stable continuity tracking across strong resonances (where $n_\text{real} \to 0$ and $n_\text{imag}$ becomes very large). At output, $|n_\text{real}|$ is applied to recover the physical convention ($n_\text{real} \ge 0$, $n_\text{imag} \ge 0$).

4. **Outlier Cleanup:**
   Isolated branch jumps (e.g., sudden flips between real-dominated and imaginary-dominated solutions) are detected by comparing each point against its neighbours. Affected points are re-solved with broader initial guesses including neighbour interpolation, then a forward re-sweep propagates the correction.

---

## 2. Comparing Results (`compare_n2.py`)

Generates a 2x2 plot comparing two refractive index CSV files: overlaid real and imaginary parts, plus their pointwise differences.

### Usage

```bash
python compare_n2.py <output1.csv> <output2.csv>
```

**Example:**
```bash
python compare_n2.py output_n2.csv jin_original_data.csv
```

Both files must have 3 columns (Frequency, n_real, n_imag). Header rows are skipped automatically.
