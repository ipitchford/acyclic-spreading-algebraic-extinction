"""EXPLORATORY SIMULATION (not a proof): finite-population base chain.

Gillespie direct-method simulation of the stochastic version of the paper's
  s'=-2sx, w'=2sx-wx, x'=(w-r)x, r'=rx,  (s0,w0,x0,r0)=(25/48,5/24,1/240,4/15).
(Internal variable names A,B,X,C correspond to S,W,X,R in the paper.)
Events: A->B rate 2AX/N, B->X rate BX/N, X->C rate CX/N. Extinction time T = first time X=0.
Deterministic: fixed seeds from SeedSequence(20261021); output floats rounded to 6 sig. digits.
Usage: python finite_population_gillespie.py [--workers 4]
"""
import json
import math
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
from multiprocessing import Pool

MASTER_SEED = 20261021
BOOT_SEED = 20261022
N_BOOT = 2000
BLOCK = 1_000_000
CONFIG = [(24_000, 200), (240_000, 60), (2_400_000, 20), (24_000_000, 20)]


def r6(v):
    if v is None:
        return None
    v = float(v)
    if v == 0.0 or not math.isfinite(v):
        return v
    return float(f"{v:.6g}")


def simulate(args):
    N, seed_seq = args
    rng = np.random.Generator(np.random.PCG64(seed_seq))
    A = 25 * N // 48
    B = 5 * N // 24
    X = N // 240
    C = 4 * N // 15
    assert A + B + X + C == N
    t = 0.0
    nev = 0
    while X > 0:
        E = rng.standard_exponential(BLOCK).tolist()
        U = rng.random(BLOCK).tolist()
        for i in range(BLOCK):
            w = 2 * A + B + C  # total rate = X*w/N
            t += E[i] * N / (X * w)
            r = U[i] * w
            if r < 2 * A:
                A -= 1
                B += 1
            elif r < 2 * A + B:
                B -= 1
                X += 1
            else:
                X -= 1
                C += 1
                if X == 0:
                    nev += i + 1
                    break
        else:
            nev += BLOCK
    return (t, C / N, A / N, B / N, nev)


def t1_of_N(N):
    # solve d^3/(3(1+d)) = 1/N for d in (0, d0], d0 = 1/4
    lo, hi = 0.0, 0.25
    f = lambda d: d ** 3 / (3 * (1 + d)) - 1.0 / N
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    d = 0.5 * (lo + hi)
    return 1.5 * (d ** -2 - 16.0)


def fit_loglog(Ns, meds):
    x = np.log(np.asarray(Ns, float))
    y = np.log(np.asarray(meds, float))
    theta, logK = np.polyfit(x, y, 1)
    return theta, math.exp(logK)


def main():
    workers = 4
    if "--workers" in sys.argv:
        workers = int(sys.argv[sys.argv.index("--workers") + 1])
    workers = min(workers, 4)
    config = list(CONFIG)

    master = np.random.SeedSequence(MASTER_SEED)
    per_N_ss = master.spawn(len(config))
    tasks = []
    for k, (N, runs) in enumerate(config):
        for j, ss in enumerate(per_N_ss[k].spawn(runs)):
            tasks.append((k, N, ss))
    with Pool(workers) as pool:
        out = pool.map(simulate, [(N, ss) for (_, N, ss) in tasks], chunksize=1)

    by_k = {k: [] for k in range(len(config))}
    for (k, N, _), res in zip(tasks, out):
        by_k[k].append(res)

    summary = []
    Ts_all = []
    for k, (N, runs) in enumerate(config):
        arr = np.array([r[:4] for r in by_k[k]])
        T, Cf, Af, Bf = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
        nev = np.array([r[4] for r in by_k[k]])
        Ts_all.append(T)
        t1 = t1_of_N(N)
        medT = float(np.median(T))
        summary.append({
            "N": N,
            "runs": runs,
            "median_T": r6(medT),
            "T_q10": r6(np.quantile(T, 0.1)),
            "T_q90": r6(np.quantile(T, 0.9)),
            "T_min": r6(T.min()),
            "T_max": r6(T.max()),
            "t1_N": r6(t1),
            "ratio_t1_over_medianT": r6(t1 / medT),
            "median_T_over_N^(2/5)": r6(medT / N ** 0.4),
            "median_final_C_over_N": r6(np.median(Cf)),
            "final_C_q10": r6(np.quantile(Cf, 0.1)),
            "final_C_q90": r6(np.quantile(Cf, 0.9)),
            "median_final_A_over_N": r6(np.median(Af)),
            "final_A_q10": r6(np.quantile(Af, 0.1)),
            "final_A_q90": r6(np.quantile(Af, 0.9)),
            "median_final_B_over_N": r6(np.median(Bf)),
            "median_events": int(np.median(nev)),
        })

    Ns = [c[0] for c in config]
    meds = [float(np.median(T)) for T in Ts_all]
    theta, K = fit_loglog(Ns, meds)
    brng = np.random.Generator(np.random.PCG64(BOOT_SEED))
    thetas = []
    for _ in range(N_BOOT):
        bm = [float(np.median(brng.choice(T, size=len(T), replace=True))) for T in Ts_all]
        thetas.append(fit_loglog(Ns, bm)[0])
    lo, hi = np.quantile(thetas, [0.025, 0.975])
    t1s = [t1_of_N(N) for N in Ns]
    th1, K1 = fit_loglog(Ns, t1s)

    result = {
        "label": "EXPLORATORY SIMULATION -- Gillespie direct method; not a proof",
        "model": {
            "events": ["A->B rate 2AX/N", "B->X rate BX/N", "X->C rate CX/N"],
            "initial_fractions": "(25/48, 5/24, 1/240, 4/15)",
            "T": "first time X=0",
            "deterministic_t1": "N*x(t1)=1 on exact solution, t=(3/2)(d^-2-16)",
        },
        "seeds": {"master": MASTER_SEED, "bootstrap": BOOT_SEED,
                  "scheme": "SeedSequence(master).spawn(len(config))[k].spawn(runs)[j]; PCG64"},
        "per_N": summary,
        "fit_median_T": {
            "form": "median T ~ K N^theta (least squares on log-log over listed N)",
            "theta": r6(theta),
            "K": r6(K),
            "theta_bootstrap_95": [r6(lo), r6(hi)],
            "n_bootstrap": N_BOOT,
        },
        "fit_t1_reference": {"theta": r6(th1), "K": r6(K1),
                             "note": "t1 ~ 1.5 (N/3)^(2/3) asymptotically"},
        "diffusion_heuristic_theta": 0.4,
    }
    fn = "finite_population_gillespie.json"
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fn), "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")


if __name__ == "__main__":
    main()
