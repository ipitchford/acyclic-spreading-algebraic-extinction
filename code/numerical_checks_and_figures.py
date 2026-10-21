#!/usr/bin/env python3
"""Independent ODE comparison and reproducible scientific figures.

The analytic proofs do not depend on these finite floating-point checks.
Run with Python 3, NumPy, SciPy and Matplotlib. Outputs are written relative
to this script so that the bundle is reproducible from any working directory.
"""

from pathlib import Path
import csv
import json
import os
import platform
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "spreading_decay_mpl_cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
EVIDENCE = ROOT / "evidence"
FIGURES.mkdir(exist_ok=True)
EVIDENCE.mkdir(exist_ok=True)

PLUM = "#794B71"
TEAL = "#218C8D"
OCHRE = "#B77936"
LAVENDER = "#9B98C7"
INK = "#29262B"
GREY = "#777277"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#AAA4AA", "axes.titleweight": "bold",
    "grid.color": "#E5E0E5", "grid.linewidth": 0.6,
    "legend.frameon": False, "pdf.fonttype": 42,
})


PDFMETA = {"CreationDate": None, "ModDate": None, "Creator": None, "Producer": None}
PNGMETA = {"Software": None}


def rnd(obj, sig=10):
    """Round floats to a fixed number of significant digits for stable output."""
    if isinstance(obj, float):
        return float(f"{obj:.{sig}g}")
    if isinstance(obj, dict):
        return {k: rnd(v, sig) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [rnd(v, sig) for v in obj]
    return obj


def fmt_rows(arr):
    return [[f"{v:.12g}" for v in row] for row in arr]


def exact_four(t):
    """State order S, W, X, R, as in the paper."""
    t = np.asarray(t, dtype=float)
    d = (16.0 + 2.0 * t / 3.0) ** -0.5
    z = 1.0 + d
    return np.array([z*z/3.0, z*(3.0-2.0*z)/3.0,
                     d**3/(3.0*z), 1.0/(3.0*z)])


def original_rhs(_t, y):
    """Reconstruct the original vector field from ordered contact reactions."""
    # (donor, destination, catalyst, coefficient)
    reactions = [(0, 1, 2, 2.0), (1, 2, 2, 1.0), (2, 3, 3, 1.0)]
    result = np.zeros(4)
    for donor, destination, catalyst, coefficient in reactions:
        flux = coefficient * y[donor] * y[catalyst]
        result[donor] -= flux
        result[destination] += flux
    return result


def perturbed_curve(t, epsilon):
    """Positive scalar reduction; log distance avoids crossing the final root."""
    a0 = 25.0/48.0
    h = 4.0/5.0
    shift = np.cbrt(-epsilon/a0)
    root = h + shift
    d0 = 1.0-root
    assert 0.0 < root < 1.0 and d0 > 0.0

    def logarithmic_distance_rhs(_t, log_distance):
        d = np.exp(log_distance[0])
        return [-a0*(d*d + 3.0*shift*d + 3.0*shift*shift)]

    solution = solve_ivp(logarithmic_distance_rhs, (0.0, float(t[-1])),
                         [np.log(d0)], t_eval=t, method="DOP853",
                         rtol=2e-12, atol=2e-12)
    assert solution.success
    d = np.exp(solution.y[0])
    x = a0*d*(d*d+3.0*shift*d+3.0*shift*shift)/(root+d)
    return x, float(3.0*a0**(1.0/3.0)*abs(epsilon)**(2.0/3.0))


def main():
    comparison_times = np.r_[0.0, np.geomspace(1e-3, 2e4, 301)]
    initial = np.array([25.0/48.0, 5.0/24.0, 1.0/240.0, 4.0/15.0])
    ode = solve_ivp(original_rhs, (0.0, float(comparison_times[-1])),
                    initial, t_eval=comparison_times, method="DOP853",
                    rtol=2e-12, atol=1e-14)
    assert ode.success, ode.message
    analytic = exact_four(comparison_times)
    absolute_error = np.abs(ode.y-analytic)
    max_error = float(absolute_error.max())
    max_active_relative_error = float((absolute_error[2]/analytic[2]).max())
    conservation_error = float(np.abs(ode.y.sum(axis=0)-1.0).max())
    assert max_error < 2e-9, max_error
    assert max_active_relative_error < 1e-5, max_active_relative_error
    assert np.all(ode.y > 0.0)
    assert conservation_error < 1e-12, conservation_error

    tail_times = np.array([1e2, 1e4, 1e6, 1e8, 1e10])
    tail_x = exact_four(tail_times)[2]
    leading = np.sqrt(3.0/8.0)
    tail_ratio = tail_times**1.5*tail_x/leading
    report = {
        "status": "PASS",
        "scope": "Finite numerical sanity checks; analytic theorems require the paper proofs.",
        "original_ode_comparison": {
            "method": "DOP853", "rtol": 2e-12, "atol": 1e-14,
            "last_time": float(comparison_times[-1]),
            "number_of_output_times": len(comparison_times),
            "rhs_evaluations": ode.nfev,
            "maximum_absolute_coordinate_error": max_error,
            "maximum_relative_X_error": max_active_relative_error,
            "maximum_total_mass_error": conservation_error,
            "minimum_numerical_coordinate": float(ode.y.min()),
        },
        "tail_sanity_check": [
            {"t": float(t), "x": float(x),
             "t_to_power_3_over_2_x_div_sqrt_3_over_8": float(ratio)}
            for t, x, ratio in zip(tail_times, tail_x, tail_ratio)
        ],
    }
    with (EVIDENCE / "original_ode_comparison.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["t", "S_exact", "W_exact", "X_exact", "R_exact",
                         "S_ode", "W_ode", "X_ode", "R_ode"])
        writer.writerows(fmt_rows(np.column_stack([comparison_times, analytic.T, ode.y.T])))

    t = np.geomspace(1.0, 2e5, 600)
    critical = exact_four(t)[2]
    epsilons = [-1e-5, 1e-5]
    perturbed = {eps: perturbed_curve(t, eps) for eps in epsilons}
    report["perturbation_rates"] = {
        str(eps): rate for eps, (_, rate) in perturbed.items()}
    with (EVIDENCE / "figure_critical_data.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["t", "X_critical", "asymptotic_X",
                         "X_epsilon_minus_1e_minus_5", "X_epsilon_plus_1e_minus_5"])
        writer.writerows(fmt_rows(np.column_stack([t, critical, leading*t**-1.5,
                                          perturbed[-1e-5][0], perturbed[1e-5][0]])))

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), layout="constrained")
    left, right = axes
    left.loglog(t, critical, color=PLUM, lw=2.2, label="Explicit solution")
    left.loglog(t, leading*t**-1.5, color=GREY, ls="--", lw=1.5,
                label=r"$\sqrt{3/8}\,t^{-3/2}$")
    sample = np.flatnonzero(comparison_times >= 1.0)[::15]
    left.loglog(comparison_times[sample], ode.y[2, sample], linestyle="none",
                marker="o", ms=4.1, mec=TEAL, mfc="white", label="Original ODE")
    left.set(title="Four states: exact algebraic extinction", xlabel="Time $t$",
             ylabel="Active fraction $x(t)$", xlim=(1.0, 2e5), ylim=(1e-9, 1.0))
    left.legend(loc="upper right", fontsize=8.6)

    right.loglog(t, critical, color=PLUM, lw=2.2, label=r"$\varepsilon=0$")
    for eps, colour, style in [(-1e-5, TEAL, "-"), (1e-5, OCHRE, "-.")]:
        sign = "-" if eps < 0 else "+"
        right.loglog(t, perturbed[eps][0], color=colour, lw=1.8, ls=style,
                     label=rf"$\varepsilon={sign}10^{{-5}}$")
    right.set(title="Nearby initial states: slow exponential tails",
              xlabel="Time $t$", ylabel="Active fraction $x(t)$",
              xlim=(1.0, 2e5), ylim=(1e-12, 1e-2))
    right.legend(loc="lower left", fontsize=8.6)
    for ax in axes:
        ax.grid(True, which="major")
    fig.savefig(FIGURES / "critical_dynamics.pdf", metadata=PDFMETA)
    fig.savefig(FIGURES / "critical_dynamics.png", dpi=190, metadata=PNGMETA)
    plt.close(fig)

    tau = np.geomspace(1e-2, 1e6, 600)
    fig, ax = plt.subplots(figsize=(8.3, 4.6), layout="constrained")
    rows = [tau]
    columns = ["normalised_time_t_over_Tm"]
    scales = {}
    for m, colour in zip([2, 4, 6, 10], [PLUM, TEAL, OCHRE, LAVENDER]):
        d0 = 1.0/(2.0*m)
        d = d0*(1.0+tau)**(-1.0/m)
        normalised_x = ((1.0+d0)/(1.0+d))*(1.0+tau)**(-1.0-1.0/m)
        scale = (m+1)*(2*m)**m/m
        scales[str(m)] = scale
        rows.append(normalised_x)
        columns.append(f"X_over_X0_m_{m}")
        ax.loglog(tau, normalised_x, color=colour, lw=2.1,
                  label=rf"$m={m}$: exponent ${m+1}/{m}$; {m+2} states")
    ax.set(title="Extinction exponents approach one in acyclic chains",
           xlabel=r"Normalised time $t/T_m$",
           ylabel=r"Normalised active fraction $x(t)/x(0)$",
           xlim=(1e-2, 1e6), ylim=(1e-10, 2.0))
    ax.grid(True, which="major")
    ax.legend(loc="lower left", fontsize=9.5)
    fig.savefig(FIGURES / "rate_hierarchy.pdf", metadata=PDFMETA)
    fig.savefig(FIGURES / "rate_hierarchy.png", dpi=190, metadata=PNGMETA)
    plt.close(fig)
    report["normalisation_time_scales_Tm"] = scales
    with (EVIDENCE / "figure_hierarchy_data.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(fmt_rows(np.column_stack(rows)))

    (EVIDENCE / "numerical_verification.json").write_text(
        json.dumps(rnd(report), indent=2, sort_keys=True)+"\n")
    print("PASS numerical checks")


if __name__ == "__main__":
    main()
