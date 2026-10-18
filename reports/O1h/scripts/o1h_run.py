#!/usr/bin/env python
"""O1h: frozen O1g Test-R top-axis diagnostic."""
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ortho_group

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "O1h"
O1G_SCRIPT = Path("/Users/terooatt/Documents/scopalign-manus-O1g/reports/O1g/scripts/o1g_run.py")
CASE, N, R, D = "R_n2759_r0", 2759, 0, 30
KS = (1, 2, 3, 5, 8, 10, 15, 20, 30)
O1G_ERR, O1G_L = 0.6037296396289642, 0.14564707632646573
KNOWN = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def import_o1g():
    if not O1G_SCRIPT.is_file():
        raise RuntimeError(f"Missing authorized O1g script: {O1G_SCRIPT}")
    spec = importlib.util.spec_from_file_location("o1g_o1h_reference", O1G_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not import authorized O1g script")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


O1G = import_o1g()


def case_data():
    z = np.load(O1G.LAT1).astype(np.float32)
    if z.shape != (5518, D):
        raise RuntimeError(f"Expected O1g latent shape (5518, {D}); found {z.shape}")
    perm = np.random.default_rng(20261030 + R).permutation(len(z))
    qstar = np.asarray(ortho_group.rvs(D, random_state=20261031 + R), dtype=np.float64)
    x = z[np.sort(perm[:N])]
    y = (z[np.sort(perm[N:2 * N])].astype(np.float64) @ qstar.T).astype(np.float32)
    return x, y, qstar


def fit_path(side, seed):
    local = O1G.local_fit(CASE, side, seed)
    if local.is_file():
        return local
    old = O1G.old_fit(CASE, side, seed)
    if seed in O1G.FIRST5 and old.is_file():
        return old
    raise RuntimeError(f"Missing cached fit: {CASE}_{side}_{seed}")


def cached_models():
    mx, my, rows = [], [], []
    for seed in O1G.SEEDS20:
        for side, bucket in (("X", mx), ("Y", my)):
            path = fit_path(side, seed)
            bucket.append(O1G.build_model(O1G.load_saved(path)))
            rows.append({"case": CASE, "side": side, "seed": int(seed), "path": str(path), "sha256": sha(path)})
    pd.DataFrame(sorted(rows, key=lambda v: (v["side"], v["seed"]))).to_csv(OUT / "fits_used.csv", index=False)
    return mx, my


def load_qhats():
    path = OUT / "Q_hats.npz"
    if not path.is_file():
        raise RuntimeError("Missing Q_hats.npz; run --stage estimate first")
    with np.load(path, allow_pickle=False) as a:
        q = np.asarray(a["Q"], dtype=np.float64)
        loss = np.asarray(a["L"], dtype=np.float64)
        qstar = np.asarray(a["Q_star"], dtype=np.float64)
        best = int(a["best"])
        variance = np.asarray(a["axis_variance"], dtype=np.float64)
        order = np.asarray(a["axis_order"], dtype=np.int64)
        sorted_axes = bool(a["axes_sorted"])
    if q.shape != (20, D, D) or loss.shape != (20,) or qstar.shape != (D, D) or not 0 <= best < 20:
        raise RuntimeError("Invalid Q_hats.npz")
    return q, loss, qstar, best, variance, order, sorted_axes


def err(q, qstar, k, order):
    return float(np.linalg.norm((q - qstar)[:, order[:k]], ord="fro") / np.sqrt(2 * k))


def disp(q, qstar, x):
    x = np.asarray(x, dtype=np.float64)
    numerator = np.mean(np.sum((x @ (q - qstar).T) ** 2, axis=1))
    denominator = 2.0 * np.mean(np.sum(x ** 2, axis=1))
    return float(np.sqrt(numerator / denominator))


def summary():
    x, _, constructed = case_data()
    q, loss, qstar, best, variance, order, sorted_axes = load_qhats()
    if not np.array_equal(qstar, constructed):
        raise RuntimeError("Saved Q_star differs from O1g construction")
    qhat = q[best]
    e30 = err(qhat, qstar, 30, order)
    return {"x": x, "Q": q, "L": loss, "Q_star": qstar, "best": best, "Q_hat": qhat, "axis_variance": variance, "axis_order": order, "axes_sorted": sorted_axes, "err_30": e30, "L_hat": float(loss[best]), "repro_err_abs": float(abs(e30 - O1G_ERR)), "repro_L_abs": float(abs(float(loss[best]) - O1G_L))}


def stopped(s):
    return s["repro_err_abs"] > 1e-6 or s["repro_L_abs"] > 1e-6


def estimate():
    OUT.mkdir(parents=True, exist_ok=True)
    required = (OUT / "Q_hats.npz", OUT / "starts.csv", OUT / "fits_used.csv")
    if not all(p.is_file() for p in required):
        if O1G.start_workers() > 3:
            raise RuntimeError(f"O1g start workers exceed limit: {O1G.start_workers()}")
        x, y, qstar = case_data()
        mx, my = cached_models()
        best, rows = O1G.estimate(O1G.Problem(mx, my, x, y, R))
        rows = sorted(rows, key=lambda v: int(v["start"]))
        q = np.stack([np.asarray(v["Q"], dtype=np.float64) for v in rows])
        loss = np.asarray([float(v["L_hat"]) for v in rows], dtype=np.float64)
        variance = np.var(np.asarray(x, dtype=np.float64), axis=0)
        axes_sorted = bool(np.all(variance[:-1] >= variance[1:]))
        order = np.arange(D, dtype=np.int64) if axes_sorted else np.argsort(-variance, kind="stable").astype(np.int64)
        np.savez(OUT / "Q_hats.npz", Q=q, L=loss, Q_star=qstar, best=np.int64(best["start"]), axis_variance=variance, axis_order=order, axes_sorted=np.bool_(axes_sorted))
        pd.DataFrame({"start": [int(v["start"]) for v in rows], "L": loss, "err_30": [err(v["Q"], qstar, 30, order) for v in rows]}).sort_values("start").to_csv(OUT / "starts.csv", index=False)
    s = summary()
    print(json.dumps({"best_start": s["best"], "L_hat": s["L_hat"], "err_30": s["err_30"], "repro_err_abs": s["repro_err_abs"], "repro_L_abs": s["repro_L_abs"], "axes_sorted": s["axes_sorted"], "start_workers": O1G.start_workers()}, sort_keys=True), flush=True)
    return s


def measure():
    s = summary()
    if stopped(s):
        raise RuntimeError("Pre-registered validity stop: the O1g case did not reproduce")
    random_q = [np.asarray(ortho_group.rvs(D, random_state=20264000 + i), dtype=np.float64) for i in range(200)]
    start_errors = np.asarray([[err(q, s["Q_star"], k, s["axis_order"]) for k in KS] for q in s["Q"]])
    identity = np.eye(D, dtype=np.float64)
    rows = []
    for j, k in enumerate(KS):
        random_errors = [err(q, s["Q_star"], k, s["axis_order"]) for q in random_q]
        rows.append({"k": k, "err_hat": err(s["Q_hat"], s["Q_star"], k, s["axis_order"]), "err_identity": err(identity, s["Q_star"], k, s["axis_order"]), "err_random_median": float(np.median(random_errors)), "n_starts_le_010": int(np.sum(start_errors[:, j] <= .10))})
    topk = pd.DataFrame(rows)
    topk.to_csv(OUT / "topk.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for values in start_errors:
        ax.plot(KS, values, color="0.65", linewidth=.8, alpha=.7)
    ax.plot(KS, topk.err_hat, color="C0", marker="o", linewidth=2.8, label="Q-hat")
    ax.plot(KS, topk.err_identity, color="C1", marker="s", linestyle="--", label="identity")
    ax.plot(KS, topk.err_random_median, color="C2", marker="^", linestyle="--", label="random median")
    ax.axhline(.10, color="black", linestyle=":", linewidth=1, label="0.10")
    ax.set(xlabel="top PCA axes (k)", ylabel="top-k rotation error", xticks=KS, ylim=(0, None))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_topk.png", dpi=180)
    plt.close(fig)
    m = {"err_top3": float(topk.loc[topk.k == 3, "err_hat"].iloc[0]), "err_top5": float(topk.loc[topk.k == 5, "err_hat"].iloc[0]), "err_top10": float(topk.loc[topk.k == 10, "err_hat"].iloc[0]), "err_top5_random": float(topk.loc[topk.k == 5, "err_random_median"].iloc[0]), "n_starts_top5_ok": int(topk.loc[topk.k == 5, "n_starts_le_010"].iloc[0]), "disp_hat": disp(s["Q_hat"], s["Q_star"], s["x"]), "disp_random": float(np.median([disp(q, s["Q_star"], s["x"]) for q in random_q]))}
    print(json.dumps(m, sort_keys=True), flush=True)
    return m


def command():
    return 'cd "$HOME/Documents/scopalign-manus-O1h" && .venv/bin/python reports/O1h/scripts/o1h_run.py'


def entry(path):
    return {"path": path, "sha256": sha(OUT / path)}


def report(s, measured=None):
    with (OUT / "requirements.lock").open("w") as f:
        subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=f, check=True)
    invalid = stopped(s)
    gated = {"repro_err_abs": s["repro_err_abs"], "repro_L_abs": s["repro_L_abs"], "axes_sorted": s["axes_sorted"]}
    if not invalid:
        if measured is None:
            topk = pd.read_csv(OUT / "topk.csv").set_index("k")
            random_q = [np.asarray(ortho_group.rvs(D, random_state=20264000 + i), dtype=np.float64) for i in range(200)]
            measured = {"err_top3": float(topk.loc[3, "err_hat"]), "err_top5": float(topk.loc[5, "err_hat"]), "err_top10": float(topk.loc[10, "err_hat"]), "err_top5_random": float(topk.loc[5, "err_random_median"]), "n_starts_top5_ok": int(topk.loc[5, "n_starts_le_010"]), "disp_hat": disp(s["Q_hat"], s["Q_star"], s["x"]), "disp_random": float(np.median([disp(q, s["Q_star"], s["x"]) for q in random_q]))}
        gated.update(measured)
        gated["reading"] = "top_recovered" if measured["err_top5"] <= .10 else ("not_recovered" if measured["err_top5"] > .25 else "partial")
        paths = ["requirements.lock", "starts.csv", "topk.csv", "fits_used.csv", "fig_topk.png", "scripts/o1h_run.py"]
    else:
        paths = ["requirements.lock", "starts.csv", "fits_used.csv", "scripts/o1h_run.py"]
    base = command()
    regen = {p: base + " --stage report" for p in paths}
    regen.update({"repro_err_abs": base + " --stage estimate", "repro_L_abs": base + " --stage estimate", "axes_sorted": base + " --stage estimate"})
    if not invalid:
        regen.update({k: base + (" --stage report" if k == "reading" else " --stage measure") for k in ("err_top3", "err_top5", "err_top10", "err_top5_random", "n_starts_top5_ok", "disp_hat", "disp_random", "reading")})
    rep = {"task_id": "O1h", "executor": "manus", "branch": "task/O1h", "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(), "seed": 20261014, "started": datetime.fromtimestamp((OUT / "starts.csv").stat().st_mtime, timezone.utc).isoformat(), "finished": datetime.now(timezone.utc).isoformat(), "status": "stopped" if invalid else "complete", "gated": gated, "deliverables": [entry(p) for p in paths], "artifacts": [entry("Q_hats.npz")], "regeneration": regen, "deviations": [], "deviations_requested": [], "observations": [f"Authorized O1g import source SHA-256: {sha(O1G_SCRIPT)}.", "fits_used.csv records the absolute path and SHA-256 of every cached drift fit loaded for this one case."], "known_failures": KNOWN, "lockfile_sha256": sha(OUT / "requirements.lock")}
    if invalid:
        rep["stop_rule"] = "Validity: repro_err_abs > 1e-6 or repro_L_abs > 1e-6 → REVISE (the case did not reproduce)."
    (OUT / "report.json").write_text(json.dumps(rep, indent=2, sort_keys=True) + "\n")
    print(json.dumps(gated, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("estimate", "measure", "report", "all"))
    stage = parser.parse_args().stage
    if stage in ("estimate", "all"):
        s = estimate()
        if stage == "all" and stopped(s):
            report(s)
            return
    if stage in ("measure", "all"):
        m = measure()
    else:
        m = None
    if stage in ("report", "all"):
        report(summary(), m)


if __name__ == "__main__":
    main()
