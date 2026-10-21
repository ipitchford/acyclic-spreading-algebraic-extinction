#!/usr/bin/env python3
"""Replay every check in this package and compare outputs byte for byte.

The package is copied to a temporary directory, every script is rerun there,
and each regenerated output file is compared with the archived copy.  Nothing
in the package itself is modified.  Exit status 0 means every output was
reproduced exactly.

Usage:  python verify_all.py            (all checks; about 90 seconds)
        python verify_all.py --skip-simulation   (omit the stochastic simulation)
"""
import filecmp
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable

# (script, outputs it regenerates, label)
JOBS = [
    ("code/check_chain_family.py", ["code/chain_verification.json"], "exact"),
    ("code/check_extensions.py", ["code/extensions_verification.json"], "exact/high-precision"),
    ("code/check_four_state.py", ["code/four_state_verification.json"], "exact/numerical"),
    ("code/numerical_checks_and_figures.py",
     ["evidence/numerical_verification.json", "evidence/original_ode_comparison.csv",
      "evidence/figure_critical_data.csv", "evidence/figure_hierarchy_data.csv",
      "figures/critical_dynamics.pdf", "figures/critical_dynamics.png",
      "figures/rate_hierarchy.pdf", "figures/rate_hierarchy.png"], "numerical"),
    ("code/finite_population_gillespie.py", ["code/finite_population_gillespie.json"], "exploratory"),
]


def main():
    skip_sim = "--skip-simulation" in sys.argv
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "package"
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns("__pycache__", "*.aux", "*.log", "*.out"))
        for script, outputs, label in JOBS:
            if skip_sim and "gillespie" in script:
                print(f"SKIP  {script}")
                continue
            if not (ROOT / script).exists():
                print(f"MISSING {script}")
                failures += 1
                continue
            for out in outputs:
                (work / out).unlink(missing_ok=True)
            proc = subprocess.run([PY, str(work / script)], cwd=work / Path(script).parent,
                                  capture_output=True, text=True)
            if proc.returncode != 0:
                print(f"FAIL  {script} (exit {proc.returncode})\n{proc.stderr[-2000:]}")
                failures += 1
                continue
            for out in outputs:
                same = (work / out).exists() and filecmp.cmp(ROOT / out, work / out, shallow=False)
                print(f"{'OK  ' if same else 'DIFF'}  [{label}] {out}")
                failures += 0 if same else 1
    print("ALL OUTPUTS REPRODUCED BYTE FOR BYTE" if failures == 0 else f"{failures} PROBLEM(S)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
