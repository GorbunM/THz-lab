#!/usr/bin/env python3
"""Fabry-Perot inversion: extract complex refractive index n2.

Input CSV (no header): Frequency (THz), |t|^2, phase (radians)
Output CSV (with header): Frequency (THz), n_real, n_imag

Usage: python extract_n2.py <input.csv> <thickness_m> [-p delay|arg] [-o output.csv]
"""

import sys
import cmath
import argparse
import numpy as np
from scipy.optimize import root

# -- physical constants -------------------------------------------------------
C_LIGHT = 299792458.0          # speed of light  [m/s]

# ===========================================================================
#  USER-CONFIGURABLE PARAMETERS 
# ===========================================================================

# File paths — set to strings to use them as defaults without CLI arguments;
# leave as None to require them on the command line.
INPUT_FILE  = "modelled_trnsm.csv"   # CSV: Frequency [THz], |t|^2, phase [rad]
OUTPUT_FILE = "output_n2.csv"        # CSV output: Frequency [THz], n_real, n_imag

# Slab thickness in metres.  Set to a float to skip the CLI argument.
# Leave as None to require it on the command line.
SLAB_THICKNESS = 5e-7          # e.g. 100e-6 for 100 µm, or 50e-6 for 50 µm

# Refractive indices of the surrounding media
N1 = 1.0                       # medium before slab
N3 = 1.0                       # medium after slab

# Solver convergence threshold
RESIDUAL_TOL = 1e-10

# Phase convention of the input data.
#   "arg"   – column 3 is arg(t), the complex argument of t
#   "delay" – column 3 is the optical delay phase, i.e. -arg(t),
#             which is positive and increases with frequency.
#             Common in THz time-domain spectroscopy.
PHASE_CONVENTION = "delay"

# ===========================================================================
#  FORWARD MODEL
# ===========================================================================

def t_model(n2, freq_Hz, d, n1=N1, n3=N3):
    """Exact Fabry-Perot transmission at normal incidence."""
    t12 = 2.0 * n1 / (n1 + n2)
    t23 = 2.0 * n2 / (n2 + n3)
    r21 = (n2 - n1) / (n2 + n1)
    r23 = (n2 - n3) / (n2 + n3)
    phi = 2.0 * np.pi * freq_Hz / C_LIGHT * n2 * d
    # Numerically stable Fabry-Perot: when Im(phi) > 0 (absorbing slab),
    # |exp(-i*phi)| = exp(Im(phi)) >> 1 and overflows.  Multiply through
    # by exp(+i*phi) to get an equivalent form where the exponential
    # argument has a non-positive imaginary part:
    #   t = t12*t23*exp(+i*phi) / (exp(2i*phi) - r21*r23)
    if phi.imag > 0:
        ep = cmath.exp(1j * phi)          # |ep| = exp(-Im(phi)) < 1, safe
        return t12 * t23 * ep / (ep * ep - r21 * r23)
    else:
        em = cmath.exp(-1j * phi)         # |em| <= 1, safe
        return t12 * t23 * em / (1.0 - r21 * r23 * em * em)


def residual(x, t_meas, freq_Hz, d):
    """Residual: [Re, Im] of (t_model - t_meas)."""
    diff = t_model(complex(x[0], x[1]), freq_Hz, d) - t_meas
    return np.array([diff.real, diff.imag])


# ===========================================================================
#  HELPERS
# ===========================================================================

def to_physical_branch(x):
    """Map to the branch with n_imag >= 0 using t(n2) = t(-n2).

    The Fabry-Perot transmission is invariant under n2 -> -n2 when
    n1 = n3, so two roots always exist.  Selecting n_imag >= 0
    provides stable continuity tracking across resonances.
    When n_imag ~ 0 (transparent region), fall back to n_real > 0.

    The output stage applies abs(n_real) to recover the physical
    branch (n_real >= 0, n_imag >= 0).
    """
    x = x.copy()
    if x[1] < -1e-12:
        x = -x
    elif abs(x[1]) < 1e-12 and x[0] < -1e-12:
        x = -x
    return x

def initial_guess_from_phase(t_meas, freq_Hz, d):
    """Single-pass approximation for n_real."""
    n_est = -np.angle(t_meas) * C_LIGHT / (2*np.pi * freq_Hz * d)
    return np.array([max(n_est, 0.5), 0.0])


# ===========================================================================
#  SOLVE AT ONE POINT with multiple initial guesses
# ===========================================================================

def try_solve(x0, t_meas, freq_Hz, d):
    """Try to solve from x0. Returns (x, res) or (None, inf)."""
    sol = root(residual, x0, args=(t_meas, freq_Hz, d),
               method='hybr', tol=1e-15, options={'maxfev': 10000})
    res = np.linalg.norm(sol.fun)
    if res > RESIDUAL_TOL:
        sol2 = root(residual, x0, args=(t_meas, freq_Hz, d),
                    method='lm', tol=1e-15, options={'maxiter': 10000})
        if np.linalg.norm(sol2.fun) < res:
            sol, res = sol2, np.linalg.norm(sol2.fun)
    if res < RESIDUAL_TOL:
        return to_physical_branch(sol.x), res
    return None, np.inf


def find_candidates(t_meas, freq_Hz, d, guesses):
    """Try all guesses, return deduplicated converged solutions."""
    unique = []
    for x0 in guesses:
        x, res = try_solve(x0, t_meas, freq_Hz, d)
        if x is not None:
            dup = any(np.linalg.norm(x - xu) < 1e-6 for xu, _ in unique)
            if not dup:
                unique.append((x.copy(), res))
    return unique


def build_guesses(prev, t_meas, freq_Hz, d, prev2=None):
    """Build initial guesses -- kept minimal to avoid spurious branches."""
    g = []
    if prev is not None:
        g.append(prev.copy())
        g.append(np.array([prev[1], prev[0]]))       # swap Re<->Im
        g.append(np.array([prev[0], -prev[1]]))       # negate Im
        g.append(prev * 0.5)
        g.append(prev * 2.0)
        if prev2 is not None:
            g.append(2.0 * prev - prev2)               # linear extrapolation

    g.append(initial_guess_from_phase(t_meas, freq_Hz, d))

    pe = abs(np.angle(t_meas)) * C_LIGHT / (2*np.pi * freq_Hz * d)
    for s in [0.5, 1.0, 2.0, 4.0]:
        g.append(np.array([0.0, max(pe * s, 1.5)]))
    g.append(np.array([1.5, 0.0]))
    g.append(np.array([0.0, 1.5]))
    return g


def build_guesses_broad(prev, next_pt, t_meas, freq_Hz, d):
    """Broader guesses for cleanup: includes neighbour interpolations."""
    g = build_guesses(prev, t_meas, freq_Hz, d)
    if prev is not None:
        for s in [1.2, 1.3, 1.5, 3.0, 0.25]:
            g.append(prev * s)
    if next_pt is not None:
        g.append(next_pt.copy())
        g.append(next_pt * 0.5)
        g.append(next_pt * 2.0)
    if prev is not None and next_pt is not None:
        g.append(0.5 * (prev + next_pt))
    return g


# ===========================================================================
#  FORWARD SWEEP with distance-continuity
# ===========================================================================

def forward_sweep(freq_THz, t_meas_arr, d):
    N = len(freq_THz)
    result = np.zeros((N, 2))
    conv   = np.zeros(N, dtype=bool)
    prev = None; prev2 = None

    for i in range(N):
        freq_Hz = freq_THz[i] * 1e12
        guesses = build_guesses(prev, t_meas_arr[i], freq_Hz, d, prev2)
        cands = find_candidates(t_meas_arr[i], freq_Hz, d, guesses)

        if not cands:
            result[i] = [np.nan, np.nan]
            conv[i] = False
            print(f"WARNING: no convergence at "
                  f"freq = {freq_THz[i]:.6f} THz  (point {i})",
                  file=sys.stderr)
            continue

        if prev is not None and len(cands) > 1:
            best, _ = min(cands, key=lambda c: np.linalg.norm(c[0] - prev))
        elif prev is not None:
            best, _ = cands[0]
        else:
            best, _ = min(cands, key=lambda c: c[0][0]**2 + c[0][1]**2)

        result[i] = best
        conv[i] = True
        prev2 = prev
        prev  = best.copy()

    return result, conv


# ===========================================================================
#  OUTLIER CLEANUP: fix isolated branch jumps
# ===========================================================================

def is_real_dominated(x, thresh=0.02):
    """True if |n_imag| < thresh * |n2|."""
    mag = max(abs(x[0]), abs(x[1]), 1e-30)
    return abs(x[1]) < thresh * mag

def is_imag_dominated(x, thresh=0.02):
    """True if |n_real| < thresh * |n2|."""
    mag = max(abs(x[0]), abs(x[1]), 1e-30)
    return abs(x[0]) < thresh * mag


def cleanup_outliers(freq_THz, t_meas_arr, d, result, conv):
    """
    Detect isolated branch jumps: if the solution type at point i
    disagrees with both neighbours, re-solve using broader guesses
    including neighbour interpolation.  Then re-run forward sweep
    from the first corrected point to propagate the fix.
    """
    N = len(freq_THz)
    first_fix = None

    for i in range(1, N-1):
        if not (conv[i] and conv[i-1] and conv[i+1]):
            continue

        xi = result[i]
        xp = result[i-1]
        xn = result[i+1]

        # Check: is point i an isolated outlier?
        neighbours_real = is_real_dominated(xp) and is_real_dominated(xn)
        neighbours_imag = is_imag_dominated(xp) and is_imag_dominated(xn)
        i_real = is_real_dominated(xi)
        i_imag = is_imag_dominated(xi)

        needs_fix = False
        if neighbours_real and not i_real:
            needs_fix = True
        elif neighbours_imag and not i_imag:
            needs_fix = True

        if not needs_fix:
            continue

        # Re-solve with broader guesses including neighbour info
        freq_Hz = freq_THz[i] * 1e12
        guesses = build_guesses_broad(xp, xn, t_meas_arr[i], freq_Hz, d)
        cands = find_candidates(t_meas_arr[i], freq_Hz, d, guesses)
        if not cands:
            continue

        # Among candidates matching the neighbour regime, pick closest
        # to the midpoint of neighbours
        mid = 0.5 * (xp + xn)
        if neighbours_real:
            same = [(x, r) for x, r in cands if is_real_dominated(x)]
        else:
            same = [(x, r) for x, r in cands if is_imag_dominated(x)]

        pool = same if same else cands
        best, _ = min(pool, key=lambda c: np.linalg.norm(c[0] - mid))

        result[i] = best
        if first_fix is None:
            first_fix = i
            print(f"Cleanup: fixed outlier at point {i} "
                  f"(freq = {freq_THz[i]:.4f} THz)")

    # If any point was fixed, re-run forward sweep from there
    if first_fix is not None:
        print(f"Cleanup: re-sweeping from point {first_fix}...")
        prev  = result[first_fix].copy()
        prev2 = result[first_fix - 1].copy() if first_fix > 0 else None

        for i in range(first_fix + 1, N):
            freq_Hz = freq_THz[i] * 1e12
            # No extrapolation (prev2=None) in resweep -- extrapolation
            # near singularities can find spurious close branches.
            guesses = build_guesses(prev, t_meas_arr[i], freq_Hz, d,
                                    prev2=None)
            cands = find_candidates(t_meas_arr[i], freq_Hz, d, guesses)

            if not cands:
                conv[i] = False
                result[i] = [np.nan, np.nan]
                continue

            if len(cands) > 1:
                best, _ = min(cands,
                              key=lambda c: np.linalg.norm(c[0] - prev))
            else:
                best, _ = cands[0]

            result[i] = best
            conv[i] = True
            prev2 = prev
            prev  = best.copy()

    return result, conv


# ===========================================================================
#  MAIN SOLVER
# ===========================================================================

def solve_n2(freq_THz, T_power, phase_rad, d, phase_convention="arg"):
    """Solve for complex n2 at each frequency.

    Parameters
    ----------
    phase_convention : str
        "arg"   – phase_rad contains arg(t).
        "delay" – phase_rad contains the optical delay phase = -arg(t).
                  The sign is flipped internally so the forward model sees
                  the correct complex transmission.
    """
    N = len(freq_THz)
    t_amp = np.sqrt(T_power)

    if phase_convention == "delay":
        # delay phase = -arg(t), so arg(t) = -phase_rad
        t_meas = t_amp * np.exp(-1j * phase_rad)
    else:
        t_meas = t_amp * np.exp(1j * phase_rad)

    result, conv = forward_sweep(freq_THz, t_meas, d)
    result, conv = cleanup_outliers(freq_THz, t_meas, d, result, conv)

    return result[:, 0], result[:, 1], conv


# ===========================================================================
#  ENTRY POINT
# ===========================================================================

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="extract_n2.py",
        description=(
            "Fabry-Perot Inversion\n"
            "======================\n"
            "Extract the complex refractive index n2 of a dielectric slab\n"
            "from its measured complex transmission coefficient.\n"
            "\n"
            "Physics: three-layer system  n1 | slab(n2, d) | n3  at normal incidence.\n"
            "The exact Fabry-Perot formula is inverted numerically at every frequency.\n"
            "\n"
            "INPUT FILE FORMAT  (CSV, no header, 3 columns):\n"
            "    Frequency [THz] ,  |t|^2  ,  phase [radians]\n"
            "\n"
            "PHASE CONVENTIONS (-p / --phase-convention):\n"
            "    arg   -- column 3 = arg(t), the complex argument of t.\n"
            "    delay -- column 3 = optical delay phase = -arg(t),\n"
            "             positive and increasing with frequency.\n"
            "             Common in THz time-domain spectroscopy.\n"
            "\n"
            "OUTPUT FILE FORMAT (CSV, with header, 3 columns):\n"
            "    Frequency [THz] ,  n_real  ,  n_imag\n"
            "    (n_imag > 0 for absorbing materials, standard physics convention)\n"
            "\n"
            "CONFIGURABLE DEFAULTS (edit at the top of this file):\n"
            "    INPUT_FILE        -- default input path\n"
            "    OUTPUT_FILE       -- default output path\n"
            "    SLAB_THICKNESS    -- omit the thickness CLI argument once set\n"
            "    PHASE_CONVENTION  -- default phase convention ('arg' or 'delay')\n"
            "    N1, N3            -- refractive indices of surrounding media\n"
            "    RESIDUAL_TOL      -- solver convergence threshold"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  python extract_n2.py data.csv 100e-6\n"
            "  python extract_n2.py data.csv 80e-6 -p delay\n"
            "  python extract_n2.py data.csv 0.5e-6 -p arg -o result.csv\n"
            "  python extract_n2.py                  # uses built-in defaults\n"
            "\n"
            "tip:\n"
            "  Set INPUT_FILE, OUTPUT_FILE, SLAB_THICKNESS, and PHASE_CONVENTION\n"
            "  at the top of this file to run it without any command-line arguments."
        ),
    )
    parser.add_argument(
        "input",
        metavar="input.csv",
        nargs="?",
        default=INPUT_FILE,
        help=(
            "CSV input file: columns = Frequency [THz], |t|^2, phase [rad]. "
            "Optional when INPUT_FILE is set at the top of this file."
        ),
    )
    parser.add_argument(
        "thickness",
        metavar="thickness_m",
        nargs="?",
        type=float,
        default=SLAB_THICKNESS,
        help=(
            "Slab thickness in metres (e.g. 100e-6 for 100 µm). "
            "Optional when SLAB_THICKNESS is set at the top of this file."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        metavar="output.csv",
        default=OUTPUT_FILE,
        help=f"Output CSV file (default: {OUTPUT_FILE})",
    )
    parser.add_argument(
        "-p", "--phase-convention",
        choices=["arg", "delay"],
        default=PHASE_CONVENTION,
        dest="phase_convention",
        help=(
            "Phase convention of the input data. "
            "'arg': column 3 = arg(t). "
            "'delay': column 3 = optical delay = -arg(t), positive and "
            "increasing with frequency (common in THz-TDS). "
            f"(default: {PHASE_CONVENTION})"
        ),
    )
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if args.input is None:
        parser.error(
            "input file not provided.\n"
            "  Pass it as a positional argument  or  set INPUT_FILE at the top of this file."
        )
    if args.thickness is None:
        parser.error(
            "slab thickness not provided.\n"
            "  Pass it as a positional argument (e.g. 100e-6 for 100 µm)\n"
            "  or set SLAB_THICKNESS at the top of this file to make it permanent."
        )

    input_file  = args.input
    d           = args.thickness
    output_file = args.output
    phase_conv  = args.phase_convention

    data = np.loadtxt(input_file, delimiter=',')
    freq_THz  = data[:, 0]
    T_power   = data[:, 1]
    phase_rad = data[:, 2]

    print(f"Input:  {len(freq_THz)} frequency points, "
          f"{freq_THz[0]:.4f} - {freq_THz[-1]:.4f} THz")
    print(f"Slab thickness  d = {d:.4g} m  ({d*1e6:.4g} um)")
    print(f"Phase convention: {phase_conv}")

    n_real, n_imag, converged = solve_n2(freq_THz, T_power, phase_rad, d,
                                         phase_convention=phase_conv)

    # The solver tracks the n_imag >= 0 branch internally for stable
    # continuity across resonances.  This can leave n_real < 0 when
    # the forward model's convention places absorption at n_imag < 0.
    # Recover the physical branch (n_real >= 0, n_imag >= 0) by
    # taking abs(n_real); since n2 and -n2 are equivalent roots,
    # abs(n_real) paired with the already-positive n_imag gives the
    # standard physics convention (n_imag > 0 = absorption).
    n_real = np.abs(n_real)

    header = "Frequency (THz),n_real,n_imag"
    out = np.column_stack([freq_THz, n_real, n_imag])
    np.savetxt(output_file, out, delimiter=',', header=header,
               comments='', fmt='%.15g')
    print(f"Output: {output_file}")

    n_ok = np.sum(converged)
    print(f"\nConverged: {n_ok}/{len(freq_THz)}  "
          f"(residual < {RESIDUAL_TOL:.0e})")
    if n_ok > 0:
        v = converged
        print(f"n_real  range: [{n_real[v].min():.6f},  "
              f"{n_real[v].max():.6f}]")
        print(f"n_imag  range: [{n_imag[v].min():.6g},  "
              f"{n_imag[v].max():.6g}]")
    n_fail = np.sum(~converged)
    if n_fail:
        print(f"\nWARNING: {n_fail} point(s) did not converge.")


if __name__ == "__main__":
    main()
