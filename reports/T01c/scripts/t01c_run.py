#!/usr/bin/env python
"""Pre-registered T01c basin-operator kill test on T00c FA latents."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import scanpy as sc
import torch
from scipy.stats import spearmanr
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP
from scjdo2.opspace import classify_basin_latent

ROOT = Path(__file__).resolve().parents[3]
REPORT = ROOT / "reports" / "T01c"
DATA = ROOT / "data" / "pbmc_v2_v3"
LATENTS = ROOT / "reports" / "T00c"
SEED, D, K = 20261014, 30, 40
CLASSES = ("", "holder", "fork", "drift")
KNOWN_FAILURES = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def config():
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=SEED, batch_size=256,
                             loss="esm", lr=5e-4, beta=1.0, esm_laplacian="hutchinson", device="cpu")


def latent_path(dataset):
    return LATENTS / f"latent_fa_{dataset}.npy"


def drift_path(dataset):
    return REPORT / f"drift_{dataset}.pt"


def ops_path(dataset):
    return REPORT / f"ops_{dataset}.npz"


def load_inputs(dataset):
    z = np.asarray(np.load(latent_path(dataset)), dtype=np.float32)
    a = sc.read_h5ad(DATA / f"pbmc_{dataset}.h5ad")
    labels = a.obs.loc[a.obs["qc_pass"].to_numpy(dtype=bool), "cell_type"].astype(str).to_numpy()
    if z.shape != (len(labels), D):
        raise RuntimeError(f"{dataset}: latent shape {z.shape} does not match QC cells ({len(labels)}, {D})")
    return z, labels


def save_drift(dataset, drift):
    torch.save({"state_dict": drift._model.state_dict(), "dim": D, "config": asdict(drift.config)}, drift_path(dataset))


def load_drift(dataset):
    saved = torch.load(drift_path(dataset), map_location="cpu", weights_only=False)
    if saved["dim"] != D:
        raise RuntimeError("saved drift dimension mismatch")
    cfg = EnergyDriftConfig(**saved["config"])
    drift = EnergyDriftAdapter(D, config=cfg)
    model = _EnergyMLP(D, cfg.hidden, cfg.depth).to("cpu")
    model.load_state_dict(saved["state_dict"])
    model.eval()
    drift._model, drift._device, drift._fitted = model, "cpu", True
    return drift


def stage_drift():
    for dataset in ("v2", "v3"):
        z, _ = load_inputs(dataset)
        save_drift(dataset, EnergyDriftAdapter(D, config=config()).fit(z, verbose=True))


def read_ops(dataset):
    with np.load(ops_path(dataset), allow_pickle=False) as data:
        return {name: data[name] for name in data.files}


def write_ops(dataset, arrays):
    np.savez_compressed(ops_path(dataset), **arrays)


def stage_operators():
    gates = {}
    for dataset in ("v2", "v3"):
        z, labels = load_inputs(dataset)
        drift = load_drift(dataset)
        jacobian = drift.jacobian_autodiff(z).astype(np.float32)
        symmetry = drift.symmetry_check(z, n_cells=100)
        write_ops(dataset, {
            "J": jacobian,
            "eig": np.linalg.eigvalsh(jacobian).astype(np.float32),
            "cell_type": np.asarray(labels, dtype=str),
            "symmetry_median": np.asarray(symmetry["median_asym"], dtype=np.float64),
            "symmetry_max": np.asarray(symmetry["max_asym"], dtype=np.float64),
            "symmetry_mean": np.asarray(symmetry["mean_asym"], dtype=np.float64),
        })
        gates[f"symmetry_{dataset}"] = float(symmetry["median_asym"])
    print(json.dumps(gates, sort_keys=True))


def majority(values):
    count = pd.Series(values).value_counts()
    if count.empty:
        return ""
    winners = count[count == count.max()].index.astype(str).tolist()
    if len(winners) != 1:
        raise RuntimeError("Unspecified tie in majority op_behaviour_basin")
    return winners[0]


def stage_basin():
    ops = {}
    for dataset in ("v2", "v3"):
        arrays = read_ops(dataset)
        z, labels = load_inputs(dataset)
        tol = 0.1 * float(np.median(np.abs(arrays["eig"])))
        a = ad.AnnData(X=np.zeros((len(labels), 1), dtype=np.float32))
        a.obs["cell_type"] = pd.Categorical(labels)
        a.obsm["X_latent"] = z
        a.uns["scjdo2"] = {"J_latent": arrays["J"]}
        classify_basin_latent(a, tol=tol, use_rep="X_latent", rule="lead_eig")
        arrays["basin_index"] = a.obs["op_basin_index"].to_numpy(dtype=np.float32)
        arrays["basin_class"] = a.obs["op_behaviour_basin"].astype(str).to_numpy(dtype=str)
        arrays["tol"] = np.asarray(tol, dtype=np.float64)
        write_ops(dataset, arrays)
        ops[dataset] = arrays
    rows = []
    for cell_type in sorted(set(ops["v2"]["cell_type"].astype(str)) | set(ops["v3"]["cell_type"].astype(str))):
        row = {"cell_type": cell_type}
        for dataset in ("v2", "v3"):
            mask = ops[dataset]["cell_type"].astype(str) == cell_type
            row[f"n_{dataset}"] = int(mask.sum())
            row[f"median_bi_{dataset}"] = float(np.median(ops[dataset]["basin_index"][mask])) if mask.any() else np.nan
            row[f"class_{dataset}"] = majority(ops[dataset]["basin_class"][mask]) if mask.any() else ""
        row["agree"] = bool(row["class_v2"] == row["class_v3"])
        rows.append(row)
    pd.DataFrame(rows, columns=["cell_type", "n_v2", "n_v3", "median_bi_v2", "median_bi_v3", "class_v2", "class_v3", "agree"]).to_csv(REPORT / "basin_table.csv", index=False)


def tested_types(ops):
    common = set(ops["v2"]["cell_type"].astype(str)) & set(ops["v3"]["cell_type"].astype(str))
    return sorted([x for x in common if (ops["v2"]["cell_type"].astype(str) == x).sum() >= 50 and (ops["v3"]["cell_type"].astype(str) == x).sum() >= 50])


def w2(a, b):
    return float(np.sqrt(ot.wasserstein_1d(a, b, p=2)))


def w1(a, b):
    return float(ot.wasserstein_1d(a, b, p=1))


def tv(a, b):
    return float(0.5 * sum(abs(np.mean(a == c) - np.mean(b == c)) for c in CLASSES))


def index_map(ops, types):
    return {d: {t: np.flatnonzero(ops[d]["cell_type"].astype(str) == t) for t in types} for d in ("v2", "v3")}


def pooled_draw(ops, idx, types, sizes, rng):
    v2 = {t: rng.choice(idx["v2"][t], sizes[t], replace=False) for t in types}
    v3 = {t: rng.choice(idx["v3"][t], sizes[t], replace=False) for t in types}
    sp, bi, cl = [], [], []
    for t in types:
        sp.append(w2(ops["v2"]["eig"][v2[t]].ravel(), ops["v3"]["eig"][v3[t]].ravel()))
        bi.append(w1(ops["v2"]["basin_index"][v2[t]], ops["v3"]["basin_index"][v3[t]]))
        cl.append(tv(ops["v2"]["basin_class"][v2[t]], ops["v3"]["basin_class"][v3[t]]))
    return float(np.mean(sp)), float(np.mean(bi)), float(np.mean(cl))


def null_draw(ops, idx, types, sizes, dataset, rng):
    sp, bi, cl = [], [], []
    for t in types:
        chosen = rng.choice(idx[dataset][t], 2 * sizes[t], replace=False)
        left, right = chosen[:sizes[t]], chosen[sizes[t]:]
        sp.append(w2(ops[dataset]["eig"][left].ravel(), ops[dataset]["eig"][right].ravel()))
        bi.append(w1(ops[dataset]["basin_index"][left], ops[dataset]["basin_index"][right]))
        cl.append(tv(ops[dataset]["basin_class"][left], ops[dataset]["basin_class"][right]))
    return float(np.mean(sp)), float(np.mean(bi)), float(np.mean(cl))


def separation_draw(ops, idx, types, sizes, dataset, rng):
    sp, bi = [], []
    for first, second in itertools.combinations(types, 2):
        size = min(sizes[first], sizes[second])
        a, b = rng.choice(idx[dataset][first], size, replace=False), rng.choice(idx[dataset][second], size, replace=False)
        sp.append(w2(ops[dataset]["eig"][a].ravel(), ops[dataset]["eig"][b].ravel()))
        bi.append(w1(ops[dataset]["basin_index"][a], ops[dataset]["basin_index"][b]))
    return float(np.mean(sp)), float(np.mean(bi))


def quotient(a, b):
    return 1.0 if b == 0 and a == 0 else (float("inf") if b == 0 else a / b)


def draw_figures(ops, types, records, table, rho):
    fig, axes = plt.subplots(1, len(types), figsize=(3.1 * len(types), 3.3), squeeze=False)
    axes = axes[0]
    for axis, cell_type in zip(axes, types):
        for dataset, color in (("v2", "C0"), ("v3", "C1")):
            mask = ops[dataset]["cell_type"].astype(str) == cell_type
            axis.hist(ops[dataset]["eig"][mask].ravel(), density=True, histtype="step", color=color, label=dataset)
        axis.set_title(cell_type)
        axis.set_xlabel("eigenvalue")
    axes[0].set_ylabel("density")
    axes[-1].legend(title="dataset")
    inset = axes[-1].inset_axes([0.46, 0.48, 0.50, 0.45])
    null = [r["spectral_w2"] for r in records if r["kind"].startswith("null_")]
    pooled = [r["spectral_w2"] for r in records if r["kind"] == "pooled"]
    inset.hist(null, density=True, color="0.7")
    inset.axvline(np.median(pooled), color="C3")
    inset.set_title("null / pooled", fontsize=7)
    fig.tight_layout()
    fig.savefig(REPORT / "fig_spectra.png", dpi=180)
    plt.close(fig)
    frame = table.set_index("cell_type").loc[types]
    x, y = frame["median_bi_v2"].to_numpy(), frame["median_bi_v3"].to_numpy()
    fig, axis = plt.subplots(figsize=(4.4, 4.0))
    axis.scatter(x, y)
    for label, x_i, y_i in zip(types, x, y):
        axis.annotate(label, (x_i, y_i), fontsize=8)
    low, high = min(x.min(), y.min()), max(x.max(), y.max())
    axis.plot([low, high], [low, high], color="0.3", linestyle="--")
    axis.set_xlabel("median basin_index v2")
    axis.set_ylabel("median basin_index v3")
    axis.set_title(f"Spearman rho = {rho:.3f}")
    fig.tight_layout()
    fig.savefig(REPORT / "fig_basin.png", dpi=180)
    plt.close(fig)


def stage_compare(write=True):
    ops = {d: read_ops(d) for d in ("v2", "v3")}
    types = tested_types(ops)
    idx = index_map(ops, types)
    sizes = {t: min(len(idx["v2"][t]), len(idx["v3"][t])) // 2 for t in types}
    rng, records = np.random.default_rng(SEED), []
    for draw in range(200):
        sp, bi, cl = pooled_draw(ops, idx, types, sizes, rng)
        records.append(dict(kind="pooled", draw=draw, spectral_w2=sp, basin_w1=bi, class_tv=cl))
    for dataset in ("v2", "v3"):
        for draw in range(100):
            sp, bi, cl = null_draw(ops, idx, types, sizes, dataset, rng)
            records.append(dict(kind=f"null_{dataset}", draw=draw, spectral_w2=sp, basin_w1=bi, class_tv=cl))
    for dataset in ("v2", "v3"):
        for draw in range(100):
            sp, bi = separation_draw(ops, idx, types, sizes, dataset, rng)
            records.append(dict(kind=f"sep_{dataset}", draw=draw, spectral_w2=sp, basin_w1=bi, class_tv=np.nan))
    draws = pd.DataFrame(records, columns=["kind", "draw", "spectral_w2", "basin_w1", "class_tv"])
    if write:
        draws.to_csv(REPORT / "draws.csv", index=False)
    median = lambda mask, key: float(np.median(draws.loc[mask, key]))
    sp_pool, sp_null = median(draws.kind == "pooled", "spectral_w2"), median(draws.kind.str.startswith("null_"), "spectral_w2")
    bi_pool, bi_null = median(draws.kind == "pooled", "basin_w1"), median(draws.kind.str.startswith("null_"), "basin_w1")
    cl_pool, cl_null = median(draws.kind == "pooled", "class_tv"), median(draws.kind.str.startswith("null_"), "class_tv")
    table = pd.read_csv(REPORT / "basin_table.csv")
    lookup = table.set_index("cell_type").loc[types]
    rho = float(spearmanr(lookup["median_bi_v2"], lookup["median_bi_v3"]).statistic)
    gate = {
        "n_cell_types_tested": int(len(types)),
        "symmetry_v2": float(ops["v2"]["symmetry_median"]), "symmetry_v3": float(ops["v3"]["symmetry_median"]),
        "spectral_w2_pooled": sp_pool, "spectral_w2_null": sp_null, "spectral_ratio": quotient(sp_pool, sp_null),
        "basin_w1_pooled": bi_pool, "basin_w1_null": bi_null, "basin_ratio": quotient(bi_pool, bi_null),
        "class_tv_pooled": cl_pool, "class_tv_null": cl_null, "class_tv_ratio": quotient(cl_pool, cl_null),
        "spectral_sep_v2": quotient(median(draws.kind == "sep_v2", "spectral_w2"), median(draws.kind == "null_v2", "spectral_w2")),
        "spectral_sep_v3": quotient(median(draws.kind == "sep_v3", "spectral_w2"), median(draws.kind == "null_v3", "spectral_w2")),
        "basin_sep_v2": quotient(median(draws.kind == "sep_v2", "basin_w1"), median(draws.kind == "null_v2", "basin_w1")),
        "basin_sep_v3": quotient(median(draws.kind == "sep_v3", "basin_w1"), median(draws.kind == "null_v3", "basin_w1")),
        "basin_index_rho": rho, "basin_class_agreement": float(lookup["agree"].mean()),
    }
    if write:
        draw_figures(ops, types, records, table, rho)
    print(json.dumps(gate, sort_keys=True))
    return gate


def stage_report():
    gate = stage_compare(write=True)
    files = ["requirements.lock", "ops_v2.npz", "ops_v3.npz", "basin_table.csv", "draws.csv", "fig_spectra.png", "fig_basin.png", "scripts/t01c_run.py"]
    base = 'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T01c/scripts/t01c_run.py'
    regen = {name: base + " --stage all" for name in files}
    regen.update({name: base + " --stage compare" for name in gate})
    report = {
        "task_id": "T01c", "executor": "manus", "branch": "task/T01c",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": SEED, "started": "2026-10-14T22:30:11-06:00", "finished": datetime.now(timezone.utc).isoformat(), "status": "complete",
        "gated": gate, "deliverables": [{"path": f, "sha256": sha256(REPORT / f)} for f in files], "regeneration": regen,
        "deviations": [], "deviations_requested": [],
        "observations": [f"Input latent_fa_v2.npy SHA-256: {sha256(latent_path('v2'))}.", f"Input latent_fa_v3.npy SHA-256: {sha256(latent_path('v3'))}."],
        "known_failures": KNOWN_FAILURES, "lockfile_sha256": sha256(REPORT / "requirements.lock"),
    }
    (REPORT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("drift", "operators", "basin", "compare", "report", "all"))
    stage = parser.parse_args().stage
    if stage in ("drift", "all"):
        stage_drift()
    if stage in ("operators", "all"):
        stage_operators()
    if stage in ("basin", "all"):
        stage_basin()
    if stage in ("compare", "all"):
        stage_compare()
    if stage in ("report", "all"):
        stage_report()


if __name__ == "__main__":
    main()
