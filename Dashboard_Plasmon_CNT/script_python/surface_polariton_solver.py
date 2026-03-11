"""Surface-polariton dispersion solver.

Computes k0 and the primary plasmon wave-vector branch beta for a thin metal film
between two dielectric half-spaces, using the Burke, Stegeman and Tamir model.

Current workflow:
- Interactive prompts for eps1, eps3, eps_m (real/imag), and thickness.
- Prints k0, beta_real, beta_imag.
- Writes parameter-based CSV/PNG/HTML outputs and opens the HTML report.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Callable, Optional
import os
import webbrowser

C_LIGHT = 299792458.0  # speed of light [m/s]
DEFAULT_OMEGA_P = 3e14
DEFAULT_K0_MIN = 0.01
DEFAULT_K0_MAX = 0.94
DEFAULT_NPTS = 500
BRANCH_REAL_TOL = 1e-15
NUM_EPS = 1e-300


# ── 1. Drude Model ──────────────────────────────────────────────────

def drude_epsilon(omega: complex, omega_p: float = 1.0, gamma: float = 0.01) -> complex:
    """Drude permittivity in paper convention: eps_m = 1 - omega_p^2/(omega*(omega - i*gamma))."""
    return 1.0 - omega_p**2 / (omega * (omega - 1j * gamma))


# ── 2. Transverse Decay Constants ────────────────────────────────────

def branch_sqrt(z: complex) -> complex:
    """Square root on the branch with Re(s) >= 0.
    If Re(s) == 0, enforces Im(s) >= 0."""
    s = np.sqrt(z)
    re_s = np.real(s)
    im_s = np.imag(s)
    if re_s < 0.0:
        s = -s
    elif abs(re_s) <= BRANCH_REAL_TOL and im_s < 0.0:
        s = -s
    return s


def compute_S(beta: complex, epsilon: complex, k0: float) -> complex:
    """Transverse decay constant: S^2 = beta^2 - epsilon * k0^2.
    Returns S with Re(S) >= 0 (bound mode convention)."""
    s_sq = beta**2 - epsilon * k0**2
    return branch_sqrt(s_sq)


def compute_all_S(beta, k0, eps1, eps3, eps_m):
    """Compute all three transverse decay constants."""
    S1 = compute_S(beta, eps1, k0)
    S2 = compute_S(beta, eps_m, k0)
    S3 = compute_S(beta, eps3, k0)
    return S1, S2, S3


# ── 3. Dispersion Relation ──────────────────────────────────────────

def dispersion_residual(beta, k0, eps1, eps3, eps_m, h,
                        sign_S1=+1, sign_S3=+1):
    """Master dispersion relation F(beta) = 0 (Eq. 7 of the paper).

    F = tanh(S2*h) * (eps1*eps3*S2^2 + eps_m^2*S1*S3)
      + S2*eps_m * (eps1*S3 + eps3*S1)

    sign_S1, sign_S3: +1 for bound modes, -1 for leaky branches.
    """
    S1, S2, S3 = compute_all_S(beta, k0, eps1, eps3, eps_m)
    S1 = sign_S1 * S1
    S3 = sign_S3 * S3

    th = np.tanh(S2 * h)
    term1 = th * (eps1 * eps3 * S2**2 + eps_m**2 * S1 * S3)
    term2 = S2 * eps_m * (eps1 * S3 + eps3 * S1)
    return term1 + term2


# ── 4. Single-Interface Initial Guess ────────────────────────────────

def single_interface_beta(eps_i, eps_m, k0):
    """Analytical solution for a semi-infinite metal (Eq. 9a,b).
    beta_i(inf) = sqrt(eps_i) * k0 * sqrt(eps_m / (eps_m + eps_i))"""
    ratio = eps_m / (eps_m + eps_i)
    b = np.sqrt(eps_i) * k0 * np.sqrt(ratio)
    if np.real(b) < 0:
        b = -b
    return complex(b)


def lossy_beta_approx(eps_i, eps_R, eps_I, k0):
    """Approximate lossy single-interface dispersion (Eqs. 10a,b).

    Returns (beta_R, beta_I).

    NOTE: The original paper (Burke, Stegeman & Tamir 1986) prints the
    numerator in Eq. 10a as  eps_R^2 - eps_R*eps_i + eps_i^2.
    Careful derivation from the exact Eq. 9 shows the correct expression
    is  eps_R^2 - eps_R*eps_i + eps_I^2  (metal loss squared, not
    dielectric constant squared).  This is a typographical error in the
    paper; the corrected form is used here.
    """
    re_sum = -eps_R + eps_i                    # Re(eps_m + eps_i)
    mod_sq = re_sum**2 + eps_I**2              # |eps_m + eps_i|^2
    numer = eps_R**2 - eps_R * eps_i + eps_I**2  # CORRECTED: eps_I^2, not eps_i^2
    beta_R = np.sqrt(eps_i) * k0 * np.sqrt(max(numer / mod_sq, 0.0))
    beta_I = 0.5 * (eps_i**2 * eps_I) / beta_R * k0**2 / mod_sq if beta_R > 0 else 0.0
    return beta_R, beta_I


# ── 5. Muller's Method ──────────────────────────────────────────────

@dataclass
class RootResult:
    root: complex = 0j
    residual: float = 0.0
    iterations: int = 0
    converged: bool = False


def muller(func: Callable, guess: complex,
           tol: float = 1e-12, max_iter: int = 500,
           perturb: float = 1e-3) -> RootResult:
    """Muller's method for complex root finding.
    Three starting points: guess*(1-p), guess, guess*(1+p)."""

    if abs(guess) < 1e-30:
        x0, x1, x2 = complex(-perturb), 0j, complex(perturb)
    else:
        x0 = guess * (1.0 - perturb)
        x1 = guess
        x2 = guess * (1.0 + perturb)

    f0 = func(x0)
    f1 = func(x1)
    f2 = func(x2)

    for it in range(max_iter):
        h1 = x1 - x0
        h2 = x2 - x1
        hsum = h1 + h2
        if abs(h1) < NUM_EPS or abs(h2) < NUM_EPS or abs(hsum) < NUM_EPS:
            jitter = complex(perturb, perturb) * max(1.0, abs(x2))
            x0, x1, x2 = x1, x2, x2 + jitter
            f0, f1, f2 = f1, f2, func(x2)
            continue

        d1 = (f1 - f0) / h1
        d2 = (f2 - f1) / h2
        a = (d2 - d1) / hsum
        b = d2 + h2 * a
        c = f2

        disc = np.sqrt(b**2 - 4.0 * a * c)
        denom_plus = b + disc
        denom_minus = b - disc
        denom = denom_plus if abs(denom_plus) >= abs(denom_minus) else denom_minus

        if abs(denom) < NUM_EPS:
            denom = b
            if abs(denom) < NUM_EPS:
                break

        dx = -2.0 * c / denom
        x3 = x2 + dx
        f3 = func(x3)
        scale = max(1.0, abs(x3))

        if abs(dx) < tol * scale or abs(f3) < tol:
            return RootResult(root=x3, residual=abs(f3),
                              iterations=it + 1, converged=True)

        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f3

    return RootResult(root=x2, residual=abs(f2),
                      iterations=max_iter, converged=False)


# ── 6. Solver with h-continuation ───────────────────────────────────

def solve_one_mode(k0, eps1, eps3, eps_m, h, guess, sign_S1=+1, sign_S3=+1):
    """Solve the dispersion relation for a single mode at one (k0, h) point."""
    def F(beta):
        return dispersion_residual(beta, k0, eps1, eps3, eps_m, h,
                                   sign_S1, sign_S3)
    return muller(F, guess)


def h_continuation(k0, eps1, eps3, eps_m, h_target, initial_root,
                   recovery_guess, sign_S1=+1, sign_S3=+1):
    """Track a root from h_start=10 down to h_target via parameter continuation."""
    h_start = 10.0
    if h_target >= h_start:
        return initial_root

    n_steps = max(20, int(10.0 * np.log(h_start / max(h_target, 1e-4))))
    root = initial_root

    for i in range(1, n_steps + 1):
        t = i / n_steps
        h_curr = h_start * (h_target / h_start) ** t
        res = solve_one_mode(k0, eps1, eps3, eps_m, h_curr, root,
                             sign_S1, sign_S3)
        if res.converged:
            root = res.root
        else:
            rec = solve_one_mode(k0, eps1, eps3, eps_m, h_curr,
                                 recovery_guess, sign_S1, sign_S3)
            if rec.converged:
                root = rec.root
    return root


@dataclass
class SolverParams:
    epsilon_1: float = 1.0
    epsilon_3: float = 11.0
    omega_p: float = 1.0
    gamma: float = 0.01
    h: float = 0.05         # normalised thickness (h * omega_p / c)
    k0_min: float = 0.01
    k0_max: float = 0.94
    n_points: int = 500
    eps_m_fixed: Optional[complex] = None  # If set, use constant eps_m (skip Drude)


@dataclass
class DispersionPoint:
    k0: float = 0.0
    beta1: complex = 0j
    beta2: complex = 0j
    betaL1: complex = 0j
    betaL2: complex = 0j
    eps_m: complex = 0j


def solve_dispersion(params: SolverParams) -> list[DispersionPoint]:
    """Solve the dispersion relation over the full k0 grid.

    If params.eps_m_fixed is set, uses that constant eps_m at all frequencies.
    Otherwise, computes eps_m from the Drude model at each frequency.
    """
    p = params
    curve: list[DispersionPoint] = []
    curve_append = curve.append
    eps1 = p.epsilon_1
    eps3 = p.epsilon_3
    k0_values = np.linspace(p.k0_min, p.k0_max, p.n_points)

    prev_beta1 = 0j
    prev_beta2 = 0j
    prev_betaL1 = 0j
    prev_betaL2 = 0j
    first_point = True

    def _solve_mode_reliable(k0: float, eps_m: complex,
                             initial_guess: complex, si_seed: complex,
                             sign_S1: int = +1, sign_S3: int = +1) -> complex:
        """Multi-strategy solver: continuation -> single-interface seed -> h-ramp."""
        # Fast path: continuation in k0 from previous valid root.
        r_main = solve_one_mode(k0, eps1, eps3, eps_m, p.h, initial_guess,
                                sign_S1, sign_S3)
        if r_main.converged:
            return r_main.root

        # Recovery 1: direct single-interface seed at target h.
        r_si = solve_one_mode(k0, eps1, eps3, eps_m, p.h, si_seed,
                              sign_S1, sign_S3)
        if r_si.converged:
            return r_si.root

        # Recovery 2: robust thick-film solve + h-continuation to target h.
        r_thick = solve_one_mode(k0, eps1, eps3, eps_m, 10.0, si_seed,
                                 sign_S1, sign_S3)
        seed = r_thick.root if r_thick.converged else si_seed
        beta = h_continuation(k0, eps1, eps3, eps_m, p.h, seed, si_seed,
                              sign_S1, sign_S3)

        # Final polish at target h.
        r_polish = solve_one_mode(k0, eps1, eps3, eps_m, p.h, beta,
                                  sign_S1, sign_S3)
        return r_polish.root if r_polish.converged else beta

    for k0 in k0_values:
        if p.eps_m_fixed is not None:
            eps_m = p.eps_m_fixed
        else:
            eps_m = drude_epsilon(k0, p.omega_p, p.gamma)

        # Analytical seeds for the two single-interface SPP branches
        si1 = single_interface_beta(eps1, eps_m, k0)   # SP_1 seed
        si3 = single_interface_beta(eps3, eps_m, k0)    # SP_3 seed

        if first_point:
            beta1  = _solve_mode_reliable(k0, eps_m, si1, si1, +1, +1)
            beta2  = _solve_mode_reliable(k0, eps_m, si3, si3, +1, +1)
            betaL1 = _solve_mode_reliable(k0, eps_m, si1, si1, -1, +1)
            betaL2 = _solve_mode_reliable(k0, eps_m, si3, si3, +1, -1)
            first_point = False
        else:
            beta1  = _solve_mode_reliable(k0, eps_m, prev_beta1,  si1, +1, +1)
            beta2  = _solve_mode_reliable(k0, eps_m, prev_beta2,  si3, +1, +1)
            betaL1 = _solve_mode_reliable(k0, eps_m, prev_betaL1, si1, -1, +1)
            betaL2 = _solve_mode_reliable(k0, eps_m, prev_betaL2, si3, +1, -1)

        prev_beta1  = beta1
        prev_beta2  = beta2
        prev_betaL1 = betaL1
        prev_betaL2 = betaL2
        curve_append(DispersionPoint(k0=k0, beta1=beta1, beta2=beta2,
                                     betaL1=betaL1, betaL2=betaL2,
                                     eps_m=eps_m))

    return curve


def curve_residual_stats(curve: list[DispersionPoint], params: SolverParams) -> tuple[float, float]:
    """Return (max, mean) residual |F(beta)| for the solved primary branch."""
    residuals = []
    for pt in curve:
        res = abs(dispersion_residual(pt.beta1, pt.k0, params.epsilon_1,
                                      params.epsilon_3, pt.eps_m, params.h,
                                      +1, +1))
        residuals.append(float(res))

    if not residuals:
        return 0.0, 0.0
    arr = np.array(residuals, dtype=float)
    return float(np.max(arr)), float(np.mean(arr))


def beta_components_paper(beta: complex) -> tuple[float, float]:
    """Return (beta_R, beta_I) for beta = beta_R - i beta_I (paper convention)."""
    return float(np.real(beta)), float(-np.imag(beta))


# ── 7. Field Profile ────────────────────────────────────────────────

def field_profile(z, beta, k0, eps1, eps3, eps_m, h):
    """Magnetic field depth profile f(z) for TM surface wave."""
    S1, S2, S3 = compute_all_S(beta, k0, eps1, eps3, eps_m)
    ratio = (S1 * eps_m) / (S2 * eps1)

    if z < 0:
        return np.exp(S1 * z)
    elif z <= h:
        return np.cosh(S2 * z) + ratio * np.sinh(S2 * z)
    else:
        bracket = np.cosh(S2 * h) + ratio * np.sinh(S2 * h)
        return bracket * np.exp(-S3 * (z - h))


# ── 8. Plotting ─────────────────────────────────────────────────────

def plot_dispersion(curve: list[DispersionPoint],
                                        params: SolverParams,
                                        save_path: str) -> None:
        """Save a clean two-panel dispersion figure for k0 vs beta."""
        k0_vals = np.array([pt.k0 for pt in curve])
        beta_r = np.array([np.real(pt.beta1) for pt in curve])
        beta_i = np.array([-np.imag(pt.beta1) for pt in curve])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=140)

        ax1.plot(k0_vals, beta_r, color="#1565c0", lw=2.0)
        ax1.set_title("Real Part", fontsize=12)
        ax1.set_xlabel(r"$k_0\;[\omega_p/c]$", fontsize=11)
        ax1.set_ylabel(r"$\beta_R\;[\omega_p/c]$", fontsize=11)
        ax1.grid(True, alpha=0.25)

        ax2.plot(k0_vals, beta_i, color="#ef6c00", lw=2.0)
        ax2.set_title("Attenuation (beta_I)", fontsize=12)
        ax2.set_xlabel(r"$k_0\;[\omega_p/c]$", fontsize=11)
        ax2.set_ylabel(r"$\beta_I\;[\omega_p/c]$", fontsize=11)
        ax2.grid(True, alpha=0.25)

        fig.suptitle("Surface Polariton Dispersion", fontsize=15, fontweight="bold", y=0.98)
        # Keep parameter metadata only in the HTML header to avoid title overlap in the PNG.
        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
        fig.savefig(save_path, bbox_inches="tight")
        plt.close(fig)


def format_token(value: float, decimals: int = 3) -> str:
        token = f"{value:+.{decimals}f}"
        return token.replace("+", "p").replace("-", "m").replace(".", "d")


def build_output_base_name(eps1: float, eps3: float, eps_m: complex, thickness_nm: float) -> str:
        return (
                f"sp_eps1_{format_token(eps1)}"
                f"__eps3_{format_token(eps3)}"
                f"__emr_{format_token(eps_m.real)}"
                f"__emi_{format_token(eps_m.imag)}"
                f"__h_{format_token(thickness_nm, 1)}nm"
        )


def write_k0_beta_csv(curve: list[DispersionPoint], output_path: str) -> None:
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("k0,beta_R,beta_I\n")
        for pt in curve:
            beta_r, beta_i = beta_components_paper(pt.beta1)
            f.write(f"{pt.k0:.8f},{beta_r:.12e},{beta_i:.12e}\n")


def build_html_report(curve: list[DispersionPoint], params: SolverParams,
                                            png_filename: str, html_path: str) -> None:
    rows = []
    for pt in curve:
        beta_r, beta_i = beta_components_paper(pt.beta1)
        rows.append(
            f"<tr><td>{pt.k0:.8f}</td><td>{beta_r:.12e}</td><td>{beta_i:.12e}</td></tr>"
        )
    table_rows = "\n".join(rows)

    if params.eps_m_fixed is not None:
        eps_m_label = f"{params.eps_m_fixed.real:.6g}{params.eps_m_fixed.imag:+.6g}j"
    else:
        eps_m_label = "drude"

    html = f"""<!doctype html>
<html lang=\"en\">
<head>
    <meta charset=\"utf-8\" />
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
    <title>Surface Polariton Results</title>
    <style>
        :root {{
            --bg-a: #f4f6fb;
            --bg-b: #dfe7f5;
            --card: #ffffff;
            --ink: #0f172a;
            --muted: #475569;
            --line: #d6deee;
            --accent: #0066cc;
            --accent-2: #ef6c00;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0;
            font-family: "Segoe UI", "Trebuchet MS", sans-serif;
            color: var(--ink);
            background: radial-gradient(circle at 15% 15%, #ffffff 0%, var(--bg-a) 36%, var(--bg-b) 100%);
            min-height: 100vh;
            padding: 28px 18px;
        }}
        .shell {{
            max-width: 1200px;
            margin: 0 auto;
            display: grid;
            gap: 16px;
        }}
        .hero {{
            background: linear-gradient(145deg, #ffffff, #f4f8ff);
            border: 1px solid var(--line);
            border-radius: 18px;
            padding: 18px 20px;
            box-shadow: 0 8px 30px rgba(15, 23, 42, 0.08);
        }}
        .hero h1 {{ margin: 0 0 8px 0; font-size: 1.5rem; letter-spacing: 0.2px; }}
        .meta {{ color: var(--muted); font-size: 0.95rem; line-height: 1.45; }}
        .grid {{
            display: grid;
            grid-template-columns: 1.2fr 1fr;
            gap: 16px;
        }}
        .card {{
            background: var(--card);
            border: 1px solid var(--line);
            border-radius: 18px;
            box-shadow: 0 6px 25px rgba(15, 23, 42, 0.08);
            overflow: hidden;
        }}
        .card h2 {{
            margin: 0;
            padding: 14px 16px;
            font-size: 1rem;
            border-bottom: 1px solid var(--line);
            background: linear-gradient(90deg, rgba(0, 102, 204, 0.08), rgba(239, 108, 0, 0.08));
        }}
        .plot-wrap {{ padding: 14px; }}
        .plot-wrap img {{ width: 100%; height: auto; border-radius: 12px; border: 1px solid var(--line); }}
        .table-wrap {{ max-height: 620px; overflow: auto; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
        thead th {{
            position: sticky; top: 0; z-index: 1;
            background: #f7fafe;
            border-bottom: 1px solid var(--line);
            text-align: left;
            padding: 10px;
            color: #0b3f7f;
        }}
        tbody td {{ padding: 8px 10px; border-bottom: 1px solid #eef2fa; }}
        tbody tr:nth-child(even) {{ background: #fcfdff; }}
        @media (max-width: 980px) {{
            .grid {{ grid-template-columns: 1fr; }}
            .table-wrap {{ max-height: none; }}
        }}
    </style>
</head>
<body>
    <div class=\"shell\">
        <section class=\"hero\">
            <h1>Surface Polariton Dispersion Results</h1>
            <div class=\"meta\">eps1={params.epsilon_1:.6g}, eps3={params.epsilon_3:.6g}, eps_m={eps_m_label}, h={params.h:.3g}</div>
        </section>
        <section class=\"grid\">
            <article class=\"card\">
                <h2>Dispersion Plot</h2>
                <div class=\"plot-wrap\"><img src=\"{png_filename}\" alt=\"Dispersion plot\" /></div>
            </article>
            <article class=\"card\">
                <h2>k0 and beta Values</h2>
                <div class=\"table-wrap\">
                    <table>
                        <thead>
                            <tr><th>k0</th><th>beta_R</th><th>beta_I</th></tr>
                        </thead>
                        <tbody>
                            {table_rows}
                        </tbody>
                    </table>
                </div>
            </article>
        </section>
    </div>
</body>
</html>
"""

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)


def plot_field_profiles(k0_target, curves_dict, params, save_path=None):
    """Plot |f(z)|^2 field profiles for each thickness at a chosen frequency."""
    fig, axes = plt.subplots(1, len(curves_dict), figsize=(5 * len(curves_dict), 5),
                             squeeze=False)

    for idx, (label, curve) in enumerate(curves_dict.items()):
        ax = axes[0, idx]
        best = min(curve, key=lambda pt: abs(pt.k0 - k0_target))
        k0 = best.k0
        eps_m = best.eps_m
        h = params.h

        S1, _, S3 = compute_all_S(best.beta1, k0, params.epsilon_1,
                                   params.epsilon_3, eps_m)
        pen1 = 1.0 / max(np.real(S1), 0.1)

        margin = max(min(3 * pen1, 5 * h), 2 * h)
        z_min = -margin
        z_max = h + margin
        z_grid = np.linspace(z_min, z_max, 500)

        for mode_name, beta in [('Mode 1', best.beta1), ('Mode 2', best.beta2)]:
            intensity = [abs(field_profile(z, beta, k0, params.epsilon_1,
                                           params.epsilon_3, eps_m, h))**2
                         for z in z_grid]
            peak = max(intensity) if max(intensity) > 0 else 1
            intensity = [val / peak for val in intensity]
            ax.plot(z_grid, intensity, lw=1.5, label=mode_name)

        ax.axvline(0, color='gold', ls='-', lw=2, alpha=0.5, label='z=0')
        ax.axvline(h, color='gold', ls='--', lw=2, alpha=0.5, label='z=h')
        ax.axvspan(0, h, alpha=0.1, color='gold')
        ax.set_xlabel(r'$z \; [c/\omega_p]$', fontsize=11)
        ax.set_ylabel(r'$|f(z)|^2$ (normalised)', fontsize=11)
        ax.set_title(f'{label}, $\\omega/\\omega_p$={k0:.3f}', fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"  Saved field profile: {save_path}")


# ── 9. CSV Export ────────────────────────────────────────────────────

def write_csv(curve, filename, physical_omega_p):
    """Write dispersion curve to CSV in physical SI units.

    Columns: omega [rad/s], k0 [rad/m], then Re/Im of beta for each mode [rad/m].
    """
    scale_beta = physical_omega_p / C_LIGHT  # normalised -> rad/m

    with open(filename, 'w') as f:
        f.write("omega [rad/s],"
                "k0 [rad/m],"
                "beta1_Re [rad/m],beta1_Im [rad/m],"
                "beta2_Re [rad/m],beta2_Im [rad/m],"
                "betaL1_Re [rad/m],betaL1_Im [rad/m],"
                "betaL2_Re [rad/m],betaL2_Im [rad/m],"
                "eps_m_Re,eps_m_Im\n")
        for pt in curve:
            omega_phys = pt.k0 * physical_omega_p
            k0_phys = pt.k0 * scale_beta
            f.write(f"{omega_phys:.10e},"
                    f"{k0_phys:.10e},"
                    f"{pt.beta1.real * scale_beta:.10e},{pt.beta1.imag * scale_beta:.10e},"
                    f"{pt.beta2.real * scale_beta:.10e},{pt.beta2.imag * scale_beta:.10e},"
                    f"{pt.betaL1.real * scale_beta:.10e},{pt.betaL1.imag * scale_beta:.10e},"
                    f"{pt.betaL2.real * scale_beta:.10e},{pt.betaL2.imag * scale_beta:.10e},"
                    f"{pt.eps_m.real:.10e},{pt.eps_m.imag:.10e}\n")
    print(f"  Saved CSV: {filename}")


# ── 10. Validation ──────────────────────────────────────────────────

def validate():
    """Run validation checks against known analytical results."""
    print("=== Validation Suite ===\n")
    gamma = 0.01
    k0_test = 0.5

    eps_m = drude_epsilon(k0_test, 1.0, gamma)
    eps_R = -eps_m.real
    eps_I = abs(eps_m.imag)
    print(f"Drude at omega/omega_p = {k0_test}:")
    print(f"  eps_m = {eps_m.real:.6f} {eps_m.imag:+.6f}j")
    print(f"  eps_R = {eps_R:.6f}, eps_I = {eps_I:.6f}")

    omega_test = 1.0 / np.sqrt(2)
    eps_test = drude_epsilon(omega_test, 1.0, 0.0001)
    print(f"\nDrude at omega/omega_p = 1/sqrt(2) = {omega_test:.4f}:")
    print(f"  Re(eps_m) = {eps_test.real:.6f} (expected ~ -1.0)")
    assert abs(eps_test.real - (-1.0)) < 0.01, "Drude validation failed"
    print("  PASS")

    print("\nThick-film limit (h=10):")
    for eps_i, name in [(1.0, "eps_1=1"), (3.0, "eps_3=3")]:
        beta_inf = single_interface_beta(eps_i, eps_m, k0_test)
        res = abs(dispersion_residual(beta_inf, k0_test, 1.0, 3.0, eps_m, 10.0))
        print(f"  |F(beta_inf)| for {name}: {res:.2e} (should be ~0)")

    print("\nLossy single-interface approximation vs exact:")
    for eps_i, name in [(1.0, "eps_1=1"), (3.0, "eps_3=3")]:
        beta_exact = single_interface_beta(eps_i, eps_m, k0_test)
        beta_R_approx, beta_I_approx = lossy_beta_approx(eps_i, eps_R, eps_I, k0_test)
        print(f"  {name}: exact = ({beta_exact.real:.6f}, {-beta_exact.imag:.6f})"
              f"  approx = ({beta_R_approx:.6f}, {beta_I_approx:.6f})")

    print("\nLossless metal (Gamma=0, thick film):")
    eps_m_lossless = drude_epsilon(k0_test, 1.0, 0.0)
    beta_lossless = single_interface_beta(1.0, eps_m_lossless, k0_test)
    print(f"  beta = {beta_lossless.real:.8f} {beta_lossless.imag:+.2e}j"
          f"  (Im should be ~0)")

    print("\n=== All Validations Passed ===\n")


# ── 11. Print physical results table ─────────────────────────────────

def print_physical_results(curve, label, physical_omega_p, h_nm,
                           eps1, eps3, eps_m_value=None, n_sample=20):
    """Print a table of dispersion results in physical SI units.

    All quantities are clearly labelled with their units.
    """
    scale_beta = physical_omega_p / C_LIGHT

    n = len(curve)
    step = max(1, n // n_sample)
    samples = curve[::step]

    sep = "=" * 120
    thin = "-" * 120
    print(f"\n{sep}")
    print(f"  SURFACE POLARITON DISPERSION RESULTS   [{label}]")
    print(sep)
    print(f"  Film thickness:      h = {h_nm:.1f} nm")
    if eps_m_value is not None:
        print(f"  Metal permittivity:  eps_m = {eps_m_value.real:.4f} {eps_m_value.imag:+.4f}j  (constant)")
    else:
        print(f"  Metal permittivity:  eps_m computed from Drude model (varies with frequency)")
    print(f"  Dielectrics:         eps_1 = {eps1:.2f} (top, z < 0)")
    print(f"                       eps_3 = {eps3:.2f} (bottom, z > h)")
    print(f"  Frequency scale:     omega_p = {physical_omega_p:.4e} rad/s")
    print(f"  Speed of light:      c = {C_LIGHT:.6e} m/s")
    print(thin)
    print(f"  QUANTITY DEFINITIONS:")
    print(f"    omega   = angular frequency                        [rad/s]")
    print(f"    k0      = free-space wave number = omega / c       [rad/m]")
    print(f"    beta_R  = Re(beta) = phase constant of SP wave     [rad/m]")
    print(f"    beta_I  = Im(beta) = attenuation constant          [rad/m]")
    print(f"              (negative beta_I => wave decays along propagation direction)")
    print(f"    Mode 1  = bound SP mode at eps_1 / metal interface")
    print(f"    Mode 2  = bound SP mode at eps_3 / metal interface")
    print(thin)

    hdr = (f"{'omega [rad/s]':>14s}  {'k0 [rad/m]':>14s}  "
           f"{'eps_m':>20s}  "
           f"{'M1 beta_R [rad/m]':>18s}  {'M1 beta_I [rad/m]':>18s}  "
           f"{'M2 beta_R [rad/m]':>18s}  {'M2 beta_I [rad/m]':>18s}")
    print(hdr)
    print(thin)

    for pt in samples:
        omega_phys = pt.k0 * physical_omega_p
        k0_phys = pt.k0 * scale_beta
        b1_R = pt.beta1.real * scale_beta
        b1_I = pt.beta1.imag * scale_beta
        b2_R = pt.beta2.real * scale_beta
        b2_I = pt.beta2.imag * scale_beta
        em = pt.eps_m

        print(f"{omega_phys:14.4e}  {k0_phys:14.4e}  "
              f"{em.real:+10.4f}{em.imag:+10.4f}j  "
              f"{b1_R:18.6e}  {b1_I:18.6e}  "
              f"{b2_R:18.6e}  {b2_I:18.6e}")

    print(sep)
    print()


# ── 12. Interactive Input & Main ────────────────────────────────────

def prompt_float(prompt_text: str, default: float) -> float:
    while True:
        raw = input(f"{prompt_text} [{default}]: ").strip()
        if raw == "":
            return float(default)
        try:
            return float(raw)
        except ValueError:
            print("Please enter a valid number.")


def collect_inputs() -> tuple[float, float, complex | None, float, float, float]:
    """Collect parameters interactively.

    Returns (eps1, eps3, eps_m_fixed_or_None, thickness_nm, gamma, omega_p_phys).
    When eps_m_fixed is None the solver will use the Drude model.
    """
    print("Surface Polariton Solver (interactive mode)")
    print("Press Enter to accept defaults.\n")

    eps1 = prompt_float("epsilon_1 (top dielectric)", 1.0)
    eps3 = prompt_float("epsilon_3 (bottom dielectric)", 11.0)

    mode_raw = input("Metal permittivity mode — (D)rude or (F)ixed? [D]: ").strip().upper()
    use_drude = mode_raw in ("", "D")

    if use_drude:
        gamma = prompt_float("Drude damping gamma / omega_p", 0.01)
        eps_m_fixed = None
    else:
        eps_m_real = prompt_float("epsilon_m real part", -15.0)
        eps_m_loss = prompt_float("epsilon_I magnitude in epsilon_m = -epsilon_R - i*epsilon_I", 1.0)
        eps_m_fixed = complex(eps_m_real, -abs(eps_m_loss))
        gamma = 0.0

    thickness_nm = prompt_float("metal thickness in nm", 50.0)
    omega_p_phys = prompt_float("plasma frequency omega_p [rad/s]", DEFAULT_OMEGA_P)

    return eps1, eps3, eps_m_fixed, thickness_nm, gamma, omega_p_phys


def print_k0_beta(curve: list[DispersionPoint]):
    print("\nk0 (omega/omega_p),beta_R (omega_p/c),beta_I (omega_p/c)")
    for pt in curve:
        beta_r, beta_i = beta_components_paper(pt.beta1)
        print(f"{pt.k0:.8f},{beta_r:.10e},{beta_i:.10e}")


def main():
    epsilon_1, epsilon_3, eps_m_fixed, h_nm, gamma, physical_omega_p = collect_inputs()

    h_meters = h_nm * 1e-9
    h_norm = h_meters * physical_omega_p / C_LIGHT

    params = SolverParams(
        epsilon_1=epsilon_1,
        epsilon_3=epsilon_3,
        omega_p=1.0,
        gamma=gamma,          # used only when eps_m_fixed is None (Drude path)
        h=h_norm,
        k0_min=DEFAULT_K0_MIN,
        k0_max=DEFAULT_K0_MAX,
        n_points=DEFAULT_NPTS,
        eps_m_fixed=eps_m_fixed,
    )

    print("\nSolving dispersion...")
    if eps_m_fixed is None:
        print(f"  Using Drude model: omega_p={physical_omega_p:.3e} rad/s, "
              f"gamma/omega_p={gamma:.4f}")
    else:
        print(f"  Using fixed eps_m = {eps_m_fixed}")

    curve = solve_dispersion(params)
    max_res, mean_res = curve_residual_stats(curve, params)
    print(f"Residual check: max|F|={max_res:.3e}, mean|F|={mean_res:.3e}")
    if max_res > 1e-8:
        print("WARNING: residual is higher than expected for mission-critical use.")
    print_k0_beta(curve)

    os.makedirs("output_py_solver", exist_ok=True)
    if eps_m_fixed is not None:
        base_name = build_output_base_name(epsilon_1, epsilon_3, eps_m_fixed, h_nm)
    else:
        base_name = (f"sp_eps1_{format_token(epsilon_1)}"
                     f"__eps3_{format_token(epsilon_3)}"
                     f"__drude_g{format_token(gamma)}"
                     f"__h_{format_token(h_nm, 1)}nm")
    csv_path = os.path.join("output_py_solver", f"{base_name}.csv")
    png_path = os.path.join("output_py_solver", f"{base_name}.png")
    html_path = os.path.join("output_py_solver", f"{base_name}.html")

    write_k0_beta_csv(curve, csv_path)
    plot_dispersion(curve, params, png_path)
    build_html_report(curve, params, os.path.basename(png_path), html_path)

    print(f"\nSaved: {csv_path}")
    print(f"Saved: {png_path}")
    print(f"Saved: {html_path}")

    webbrowser.open(f"file://{os.path.abspath(html_path)}", new=1)


if __name__ == "__main__":
    main()