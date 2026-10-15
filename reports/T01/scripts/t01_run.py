#!/usr/bin/env python
"""Pre-registered T01 basin-operator kill test."""
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
REPORT = ROOT / "reports" / "T01"
DATA = ROOT / "data" / "pbmc_v2_v3"
T00B = ROOT / "reports" / "T00b"
SEED = 20261014
D = 30
N_DRAWS = 200
N_NULL = 100
N_SEP = 100
CLASSES = ("", "holder", "fork", "drift")
KNOWN_FAILURES = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]


def cfg() -> EnergyDriftConfig:
    return EnergyDriftConfig(
        hidden=96, depth=3, n_epochs=1500, seed=SEED, batch_size=256,
        loss="esm", lr=5e-4, beta=1.0, esm_laplacian="hutchinson", device="cpu",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def latent_path(dataset: str) -> Path:
    return T00B / f"latent_fa_{dataset}.npy"


def drift_path(dataset: str) -> Path:
    return REPORT / f"drift_{dataset}.pt"


def ops_path(dataset: str) -> Path:
    return REPORT / f"ops_{dataset}.npz"


def load_inputs(dataset: str) -> tuple[np.ndarray, np.ndarray]:
    latent = np.asarray(np.load(latent_path(dataset)), dtype=np.float32)
    adata = sc.read_h5ad(DATA / f"pbmc_{dataset}.h5ad")
    qc = adata.obs["qc_pass"].to_numpy(dtype=bool)
    labels = adata.obs.loc[qc, "cell_type"].astype(str).to_numpy()
    if latent.shape != (len(labels), D):
        raise RuntimeError(
            f"{dataset}: latent shape {latent.shape} does not match QC rows ({len(labels)}, {D})"
        )
    return latent, labels


def save_drift(dataset: str, drift: EnergyDriftAdapter) -> None:
    torch.save({
        "state_dict": drift._model.state_dict(),
        "dim": D,
        "config": asdict(drift.config),
    }, drift_path(dataset))


def load_drift(dataset: str) -> EnergyDriftAdapter:
    payload = torch.load(drift_path(dataset), map_location="cpu", weights_only=False)
    if int(payload["dim"]) != D:
        raise RuntimeError(f"{dataset}: saved drift dimension is not {D}")
    saved = payload["config"]
    config = EnergyDriftConfig(**saved)
    drift = EnergyDriftAdapter(dim=D, config=config)
    model = _EnergyMLP(D, config.hidden, config.depth).to("cpu")
    model.load_state_dict(payload["state_dict"])
    model.eval()
    drift._model = model
    drift._device = "cpu"
    drift._fitted = True
    return drift


def stage_drift() -> None:
    for dataset in ("v2", "v3"):
        z, _ = load_inputs(dataset)
        drift = EnergyDriftAdapter(dim=D, config=cfg()).fit(z, verbose=True)
        save_drift(dataset, drift)


def load_ops(dataset: str) -> dict[str, np.ndarray]:
    with np.load(ops_path(dataset), allow_pickle=False) as loaded:
        return {key: loaded[key] for key in loaded.files}


def save_ops(dataset: str, arrays: dict[str, np.ndarray]) -> None:
    np.savez_compressed(ops_path(dataset), **arrays)


def stage_operators() -> dict[str, float]:
    result = {}
    for dataset in ("v2", "v3"):
        z, labels = load_inputs(dataset)
        drift = load_drift(dataset)
        jacobian = drift.jacobian_autodiff(z)
        eigenvalues = np.linalg.eigvalsh(jacobian).astype(np.float32)
        symmetry = drift.symmetry_check(z, n_cells=100)
        save_ops(dataset, {
            "J": jacobian.astype(np.float32),
            "eig": eigenvalues,
            "cell_type": np.asarray(labels, dtype=str),
            "symmetry_median": np.asarray(symmetry["median_asym"], dtype=np.float64),
            "symmetry_max": np.asarray(symmetry["max_asym"], dtype=np.float64),
            "symmetry_mean": np.asarray(symmetry["mean_asym"], dtype=np.float64),
        })
        result[f"symmetry_{dataset}"] = float(symmetry["median_asym"])
    print(json.dumps(result, sort_keys=True))
    return result


def majority(values: np.ndarray) -> str:
    counts = pd.Series(values).value_counts()
    if counts.empty:
        return ""
    winners = counts[counts == counts.max()].index.astype(str).tolist()
    if len(winners) != 1:
        raise RuntimeError("Unspecified tie in majority op_behaviour_basin")
    return winners[0]


def stage_basin() -> None:
    results = {}
    for dataset in ("v2", "v3"):
        arrays = load_ops(dataset)
        z, labels = load_inputs(dataset)
        eigenvalues = np.asarray(arrays["eig"], dtype=np.float64)
        tol = 0.1 * float(np.median(np.abs(eigenvalues)))
        operator_adata = ad.AnnData(X=np.zeros((len(labels), 1), dtype=np.float32))
        operator_adata.obs["cell_type"] = pd.Categorical(labels)
        operator_adata.obsm["X_latent"] = z
        operator_adata.uns["scjdo2"] = {"J_latent": arrays["J"]}
        classify_basin_latent(operator_adata, tol=tol, use_rep="X_latent", rule="lead_eig")
        arrays["basin_index"] = operator_adata.obs["op_basin_index"].to_numpy(dtype=np.float32)
        arrays["basin_class"] = operator_adata.obs["op_behaviour_basin"].astype(str).to_numpy(dtype=str)
        arrays["tol"] = np.asarray(tol, dtype=np.float64)
        save_ops(dataset, arrays)
        results[dataset] = arrays

    types = sorted(set(results["v2"]["cell_type"].astype(str)) | set(results["v3"]["cell_type"].astype(str)))
    rows = []
    for cell_type in types:
        row = {"cell_type": cell_type}
        for dataset in ("v2", "v3"):
            mask = results[dataset]["cell_type"].astype(str) == cell_type
            row[f"n_{dataset}"] = int(mask.sum())
            row[f"median_bi_{dataset}"] = float(np.median(results[dataset]["basin_index"][mask])) if mask.any() else np.nan
            row[f"class_{dataset}"] = majority(results[dataset]["basin_class"][mask]) if mask.any() else ""
        row["agree"] = bool(row["class_v2"] == row["class_v3"])
        rows.append(row)
    pd.DataFrame(rows, columns=[
        "cell_type", "n_v2", "n_v3", "median_bi_v2", "median_bi_v3", "class_v2", "class_v3", "agree"
    ]).to_csv(REPORT / "basin_table.csv", index=False)


def tested_types(ops: dict[str, dict[str, np.ndarray]]) -> list[str]:
    common = set(ops["v2"]["cell_type"].astype(str)) & set(ops["v3"]["cell_type"].astype(str))
    return sorted([
        cell_type for cell_type in common
        if (ops["v2"]["cell_type"].astype(str) == cell_type).sum() >= 50
        and (ops["v3"]["cell_type"].astype(str) == cell_type).sum() >= 50
    ])


def w2(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.sqrt(ot.wasserstein_1d(left, right, p=2)))


def w1(left: np.ndarray, right: np.ndarray) -> float:
    return float(ot.wasserstein_1d(left, right, p=1))


def class_tv(left: np.ndarray, right: np.ndarray) -> float:
    return 0.5 * sum(abs(np.mean(left == category) - np.mean(right == category)) for category in CLASSES)


def type_indices(ops: dict[str, dict[str, np.ndarray]], kinds: list[str]) -> dict[str, dict[str, np.ndarray]]:
    return {
        dataset: {
            cell_type: np.flatnonzero(ops[dataset]["cell_type"].astype(str) == cell_type)
            for cell_type in kinds
        }
        for dataset in ("v2", "v3")
    }


def evaluate_same(ops, idx_left, idx_right, kinds: list[str]) -> tuple[float, float, float]:
    spectral, basin, tv = [], [], []
    for cell_type in kinds:
        l, r = idx_left[cell_type], idx_right[cell_type]
        spectral.append(w2(ops["v2"]["eig"][l].ravel(), ops["v3"]["eig"][r].ravel()))
        basin.append(w1(ops["v2"]["basin_index"][l], ops["v3"]["basin_index"][r]))
        tv.append(class_tv(ops["v2"]["basin_class"][l], ops["v3"]["basin_class"][r]))
    return float(np.mean(spectral)), float(np.mean(basin)), float(np.mean(tv))


def same_draw(ops, indices, kinds, mt, rng):
    v2 = {t: rng.choice(indices["v2"][t], size=mt[t], replace=False) for t in kinds}
    v3 = {t: rng.choice(indices["v3"][t], size=mt[t], replace=False) for t in kinds}
    return evaluate_same(ops, v2, v3, kinds)


def null_draw(ops, indices, kinds, mt, dataset, rng):
    spectral, basin, tv = [], [], []
    for cell_type in kinds:
        chosen = rng.choice(indices[dataset][cell_type], size=2 * mt[cell_type], replace=False)
        left, right = chosen[:mt[cell_type]], chosen[mt[cell_type]:]
        spectral.append(w2(ops[dataset]["eig"][left].ravel(), ops[dataset]["eig"][right].ravel()))
        basin.append(w1(ops[dataset]["basin_index"][left], ops[dataset]["basin_index"][right]))
        tv.append(class_tv(ops[dataset]["basin_class"][left], ops[dataset]["basin_class"][right]))
    return float(np.mean(spectral)), float(np.mean(basin)), float(np.mean(tv))


def separation_draw(ops, indices, kinds, mt, dataset, rng):
    spectral, basin = [], []
    for first, second in itertools.combinations(kinds, 2):
        size = min(mt[first], mt[second])
        a = rng.choice(indices[dataset][first], size=size, replace=False)
        b = rng.choice(indices[dataset][second], size=size, replace=False)
        spectral.append(w2(ops[dataset]["eig"][a].ravel(), ops[dataset]["eig"][b].ravel()))
        basin.append(w1(ops[dataset]["basin_index"][a], ops[dataset]["basin_index"][b]))
    return float(np.mean(spectral)), float(np.mean(basin))


def ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return 1.0 if numerator == 0 else float("inf")
    return numerator / denominator


def write_figures(ops, kinds, rows, basin_frame, rho) -> None:
    fig, axes = plt.subplots(1, len(kinds), figsize=(3.1 * len(kinds), 3.3), squeeze=False)
    axes = axes[0]
    for axis, cell_type in zip(axes, kinds):
        for dataset, color in (("v2", "C0"), ("v3", "C1")):
            mask = ops[dataset]["cell_type"].astype(str) == cell_type
            axis.hist(ops[dataset]["eig"][mask].ravel(), density=True, histtype="step", color=color, label=dataset)
        axis.set_title(cell_type)
        axis.set_xlabel("eigenvalue")
    axes[0].set_ylabel("density")
    axes[-1].legend(title="dataset")
    pooled = [row["spectral_w2"] for row in rows if row["kind"] == "pooled"]
    null = [row["spectral_w2"] for row in rows if row["kind"].startswith("null_")]
    inset = axes[-1].inset_axes([0.46, 0.48, 0.50, 0.45])
    inset.hist(null, density=True, color="0.7")
    inset.axvline(np.median(pooled), color="C3")
    inset.set_title("null / pooled", fontsize=7)
    fig.tight_layout()
    fig.savefig(REPORT / "fig_spectra.png", dpi=180)
    plt.close(fig)

    subset = basin_frame.set_index("cell_type").loc[kinds]
    x, y = subset["median_bi_v2"].to_numpy(), subset["median_bi_v3"].to_numpy()
    fig, axis = plt.subplots(figsize=(4.4, 4.0))
    axis.scatter(x, y)
    for cell_type, x_i, y_i in zip(kinds, x, y):
        axis.annotate(cell_type, (x_i, y_i), fontsize=8)
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    axis.plot([lo, hi], [lo, hi], color="0.3", linestyle="--")
    axis.set_xlabel("median basin_index v2")
    axis.set_ylabel("median basin_index v3")
    axis.set_title(f"Spearman rho = {rho:.3f}")
    fig.tight_layout()
    fig.savefig(REPORT / "fig_basin.png", dpi=180)
    plt.close(fig)


def stage_compare(write_outputs: bool = True) -> dict[str, float]:
    ops = {dataset: load_ops(dataset) for dataset in ("v2", "v3")}
    kinds = tested_types(ops)
    indices = type_indices(ops, kinds)
    mt = {t: int(min(len(indices["v2"][t]), len(indices["v3"][t])) // 2) for t in kinds}
    rng = np.random.default_rng(SEED)
    rows = []
    for draw in range(N_DRAWS):
        spectral, basin, tv = same_draw(ops, indices, kinds, mt, rng)
        rows.append(dict(kind="pooled", draw=draw, spectral_w2=spectral, basin_w1=basin, class_tv=tv))
    for dataset in ("v2", "v3"):
        for draw in range(N_NULL):
            spectral, basin, tv = null_draw(ops, indices, kinds, mt, dataset, rng)
            rows.append(dict(kind=f"null_{dataset}", draw=draw, spectral_w2=spectral, basin_w1=basin, class_tv=tv))
    for dataset in ("v2", "v3"):
        for draw in range(N_SEP):
            spectral, basin = separation_draw(ops, indices, kinds, mt, dataset, rng)
            rows.append(dict(kind=f"sep_{dataset}", draw=draw, spectral_w2=spectral, basin_w1=basin, class_tv=np.nan))
    draws = pd.DataFrame(rows, columns=["kind", "draw", "spectral_w2", "basin_w1", "class_tv"])
    if write_outputs:
        draws.to_csv(REPORT / "draws.csv", index=False)
    def median(kind, column):
        return float(np.median(draws.loc[draws.kind == kind, column]))
    spectral_pooled, spectral_null = median("pooled", "spectral_w2"), float(np.median(draws.loc[draws.kind.str.startswith("null_"), "spectral_w2"]))
    basin_pooled, basin_null = median("pooled", "basin_w1"), float(np.median(draws.loc[draws.kind.str.startswith("null_"), "basin_w1"]))
    class_pooled, class_null = median("pooled", "class_tv"), float(np.median(draws.loc[draws.kind.str.startswith("null_"), "class_tv"]))
    basin_frame = pd.read_csv(REPORT / "basin_table.csv")
    med_v2 = basin_frame.set_index("cell_type").loc[kinds, "median_bi_v2"].to_numpy()
    med_v3 = basin_frame.set_index("cell_type").loc[kinds, "median_bi_v3"].to_numpy()
    rho = float(spearmanr(med_v2, med_v3).statistic)
    agreement = float(basin_frame.set_index("cell_type").loc[kinds, "agree"].mean())
    gated = {
        "n_cell_types_tested": int(len(kinds)),
        "symmetry_v2": float(ops["v2"]["symmetry_median"]),
        "symmetry_v3": float(ops["v3"]["symmetry_median"]),
        "spectral_w2_pooled": spectral_pooled,
        "spectral_w2_null": spectral_null,
        "spectral_ratio": ratio(spectral_pooled, spectral_null),
        "basin_w1_pooled": basin_pooled,
        "basin_w1_null": basin_null,
        "basin_ratio": ratio(basin_pooled, basin_null),
        "class_tv_pooled": class_pooled,
        "class_tv_null": class_null,
        "class_tv_ratio": ratio(class_pooled, class_null),
        "spectral_sep_v2": ratio(median("sep_v2", "spectral_w2"), median("null_v2", "spectral_w2")),
        "spectral_sep_v3": ratio(median("sep_v3", "spectral_w2"), median("null_v3", "spectral_w2")),
        "basin_sep_v2": ratio(median("sep_v2", "basin_w1"), median("null_v2", "basin_w1")),
        "basin_sep_v3": ratio(median("sep_v3", "basin_w1"), median("null_v3", "basin_w1")),
        "basin_index_rho": rho,
        "basin_class_agreement": agreement,
    }
    if write_outputs:
        write_figures(ops, kinds, rows, basin_frame, rho)
    print(json.dumps(gated, sort_keys=True))
    return gated


def stage_report() -> None:
    gated = stage_compare(write_outputs=True)
    files = [
        "requirements.lock", "ops_v2.npz", "ops_v3.npz", "basin_table.csv", "draws.csv",
        "fig_spectra.png", "fig_basin.png", "scripts/t01_run.py",
    ]
    regeneration = {path: f'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T01/scripts/t01_run.py --stage all' for path in files}
    for name in gated:
        regeneration[name] = 'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T01/scripts/t01_run.py --stage compare'
    report = {
        "task_id": "T01", "executor": "manus", "branch": "task/T01",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": SEED, "started": "2026-10-14T21:31:10-06:00",
        "finished": datetime.now(timezone.utc).isoformat(), "status": "complete",
        "gated": gated,
        "deliverables": [{"path": path, "sha256": sha256(REPORT / path)} for path in files],
        "regeneration": regeneration,
        "deviations": [], "deviations_requested": [],
        "observations": [
            f"Input latent_fa_v2.npy SHA-256: {sha256(latent_path('v2'))}.",
            f"Input latent_fa_v3.npy SHA-256: {sha256(latent_path('v3'))}.",
        ],
        "known_failures": KNOWN_FAILURES,
        "lockfile_sha256": sha256(REPORT / "requirements.lock"),
    }
    (REPORT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("drift", "operators", "basin", "compare", "report", "all"))
    args = parser.parse_args()
    if args.stage in ("drift", "all"):
        stage_drift()
    if args.stage in ("operators", "all"):
        stage_operators()
    if args.stage in ("basin", "all"):
        stage_basin()
    if args.stage in ("compare", "all"):
        stage_compare()
    if args.stage in ("report", "all"):
        stage_report()


if __name__ == "__main__":
    main()
