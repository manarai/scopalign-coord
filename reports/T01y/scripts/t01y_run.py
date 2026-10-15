#!/usr/bin/env python
"""Pre-registered T01y scale exploration from committed T01c operators."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
REPORT = ROOT / "reports" / "T01y"
INPUT = ROOT / "reports" / "T01c"
SEED = 20261014
CLASSES = ("", "holder", "fork", "drift")
REFERENCE = {
    "class_tv_ratio": 6.508187871784667,
    "spectral_sep_v2": 2.958350209622332,
    "spectral_sep_v3": 13.765335171400633,
    "basin_sep_v2": 1.2213431209981245,
    "basin_sep_v3": 3.5299139192227402,
}
KNOWN_FAILURES = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]


def digest(path):
    result = hashlib.sha256()
    with open(path, "rb") as handle:
        for part in iter(lambda: handle.read(1 << 20), b""):
            result.update(part)
    return result.hexdigest()


def input_path(dataset):
    return INPUT / f"ops_{dataset}.npz"


def load_data():
    data = {}
    for dataset, n_cells in (("v2", 992), ("v3", 1113)):
        with np.load(input_path(dataset), allow_pickle=False) as arrays:
            eig = np.asarray(arrays["eig"], dtype=np.float64)
            labels = arrays["cell_type"].astype(str)
        if eig.shape != (n_cells, 30) or labels.shape != (n_cells,) or not np.isrealobj(eig):
            raise RuntimeError(f"{dataset}: unexpected T01c input arrays")
        scale = float(np.median(np.abs(eig)))
        if scale == 0:
            raise RuntimeError(f"{dataset}: scale is zero")
        scaled_eig = eig / scale
        lead = scaled_eig.max(axis=1)
        basin = -lead
        classes = np.where(lead > 0.1, "fork", np.where(lead < -0.1, "holder", "drift"))
        data[dataset] = {"eig": scaled_eig, "labels": labels, "scale": scale, "basin": basin, "classes": classes}
    return data


def w2(left, right):
    return float(np.sqrt(ot.wasserstein_1d(left, right, p=2)))


def w1(left, right):
    return float(ot.wasserstein_1d(left, right, p=1))


def tv(left, right):
    return float(0.5 * sum(abs(np.mean(left == item) - np.mean(right == item)) for item in CLASSES))


def ratio(numerator, denominator):
    if denominator == 0:
        return 1.0 if numerator == 0 else float("inf")
    return numerator / denominator


def selected(data):
    types = sorted([cell_type for cell_type in set(data["v2"]["labels"]) & set(data["v3"]["labels"])
                    if (data["v2"]["labels"] == cell_type).sum() >= 50 and (data["v3"]["labels"] == cell_type).sum() >= 50])
    indexes = {dataset: {cell_type: np.flatnonzero(data[dataset]["labels"] == cell_type) for cell_type in types}
               for dataset in ("v2", "v3")}
    sizes = {cell_type: min(len(indexes["v2"][cell_type]), len(indexes["v3"][cell_type])) // 2 for cell_type in types}
    return types, indexes, sizes


def pair_metrics(data, v2_indices, v3_indices, types):
    spectral, basin, classes = [], [], []
    for cell_type in types:
        spectral.append(w2(data["v2"]["eig"][v2_indices[cell_type]].ravel(), data["v3"]["eig"][v3_indices[cell_type]].ravel()))
        basin.append(w1(data["v2"]["basin"][v2_indices[cell_type]], data["v3"]["basin"][v3_indices[cell_type]]))
        classes.append(tv(data["v2"]["classes"][v2_indices[cell_type]], data["v3"]["classes"][v3_indices[cell_type]]))
    return float(np.mean(spectral)), float(np.mean(basin)), float(np.mean(classes))


def within_metrics(data, left, right, types, dataset):
    spectral, basin, classes = [], [], []
    for cell_type in types:
        spectral.append(w2(data[dataset]["eig"][left[cell_type]].ravel(), data[dataset]["eig"][right[cell_type]].ravel()))
        basin.append(w1(data[dataset]["basin"][left[cell_type]], data[dataset]["basin"][right[cell_type]]))
        classes.append(tv(data[dataset]["classes"][left[cell_type]], data[dataset]["classes"][right[cell_type]]))
    return float(np.mean(spectral)), float(np.mean(basin)), float(np.mean(classes))


def separation_metrics(data, types, sizes, indexes, dataset, rng):
    spectral, basin = [], []
    for first, second in itertools.combinations(types, 2):
        draw_size = min(sizes[first], sizes[second])
        left = rng.choice(indexes[dataset][first], draw_size, replace=False)
        right = rng.choice(indexes[dataset][second], draw_size, replace=False)
        spectral.append(w2(data[dataset]["eig"][left].ravel(), data[dataset]["eig"][right].ravel()))
        basin.append(w1(data[dataset]["basin"][left], data[dataset]["basin"][right]))
    return float(np.mean(spectral)), float(np.mean(basin))


def make_draws(data, types, indexes, sizes):
    rng = np.random.default_rng(SEED)
    rows = []
    for draw in range(200):
        v2 = {cell_type: rng.choice(indexes["v2"][cell_type], sizes[cell_type], replace=False) for cell_type in types}
        v3 = {cell_type: rng.choice(indexes["v3"][cell_type], sizes[cell_type], replace=False) for cell_type in types}
        spectral, basin, classes = pair_metrics(data, v2, v3, types)
        rows.append(dict(kind="pooled", draw=draw, spectral_w2=spectral, basin_w1=basin, class_tv=classes))
    for dataset in ("v2", "v3"):
        for draw in range(100):
            left, right = {}, {}
            for cell_type in types:
                values = rng.choice(indexes[dataset][cell_type], 2 * sizes[cell_type], replace=False)
                left[cell_type], right[cell_type] = values[:sizes[cell_type]], values[sizes[cell_type]:]
            spectral, basin, classes = within_metrics(data, left, right, types, dataset)
            rows.append(dict(kind=f"null_{dataset}", draw=draw, spectral_w2=spectral, basin_w1=basin, class_tv=classes))
    for dataset in ("v2", "v3"):
        for draw in range(100):
            spectral, basin = separation_metrics(data, types, sizes, indexes, dataset, rng)
            rows.append(dict(kind=f"sep_{dataset}", draw=draw, spectral_w2=spectral, basin_w1=basin, class_tv=np.nan))
    return pd.DataFrame(rows, columns=["kind", "draw", "spectral_w2", "basin_w1", "class_tv"])


def median(frame, mask, field):
    return float(np.median(frame.loc[mask, field]))


def render(data, types, draws):
    fig, axes = plt.subplots(1, len(types), figsize=(3.1 * len(types), 3.3), squeeze=False)
    axes = axes[0]
    for axis, cell_type in zip(axes, types):
        for dataset, color in (("v2", "C0"), ("v3", "C1")):
            mask = data[dataset]["labels"] == cell_type
            axis.hist(data[dataset]["eig"][mask].ravel(), density=True, histtype="step", color=color, label=dataset)
        axis.set_title(cell_type)
        axis.set_xlabel("scaled eigenvalue")
    axes[0].set_ylabel("density")
    axes[-1].legend(title="dataset")
    inset = axes[-1].inset_axes([0.46, 0.48, 0.50, 0.45])
    inset.hist(draws.loc[draws.kind.str.startswith("null_"), "spectral_w2"], density=True, color="0.7")
    inset.axvline(median(draws, draws.kind == "pooled", "spectral_w2"), color="C3")
    inset.set_title("null / pooled", fontsize=7)
    fig.tight_layout()
    fig.savefig(REPORT / "fig_spectra_scaled.png", dpi=180)
    plt.close(fig)


def stage_compare(write=True):
    data = load_data()
    types, indexes, sizes = selected(data)
    draws = make_draws(data, types, indexes, sizes)
    if write:
        draws.to_csv(REPORT / "draws_scaled.csv", index=False)
        render(data, types, draws)
    spectral_pooled = median(draws, draws.kind == "pooled", "spectral_w2")
    spectral_null = median(draws, draws.kind.str.startswith("null_"), "spectral_w2")
    basin_pooled = median(draws, draws.kind == "pooled", "basin_w1")
    basin_null = median(draws, draws.kind.str.startswith("null_"), "basin_w1")
    class_pooled = median(draws, draws.kind == "pooled", "class_tv")
    class_null = median(draws, draws.kind.str.startswith("null_"), "class_tv")
    sep = {
        "spectral_sep_v2": ratio(median(draws, draws.kind == "sep_v2", "spectral_w2"), median(draws, draws.kind == "null_v2", "spectral_w2")),
        "spectral_sep_v3": ratio(median(draws, draws.kind == "sep_v3", "spectral_w2"), median(draws, draws.kind == "null_v3", "spectral_w2")),
        "basin_sep_v2": ratio(median(draws, draws.kind == "sep_v2", "basin_w1"), median(draws, draws.kind == "null_v2", "basin_w1")),
        "basin_sep_v3": ratio(median(draws, draws.kind == "sep_v3", "basin_w1"), median(draws, draws.kind == "null_v3", "basin_w1")),
    }
    gated = {
        "scale_v2": data["v2"]["scale"], "scale_v3": data["v3"]["scale"], "scale_ratio": ratio(data["v2"]["scale"], data["v3"]["scale"]),
        "spectral_w2_pooled_s": spectral_pooled, "spectral_w2_null_s": spectral_null, "spectral_ratio_s": ratio(spectral_pooled, spectral_null),
        "basin_w1_pooled_s": basin_pooled, "basin_w1_null_s": basin_null, "basin_ratio_s": ratio(basin_pooled, basin_null),
        "class_tv_ratio_s": ratio(class_pooled, class_null),
        "max_sep_change": float(max(abs(sep[key] - REFERENCE[key]) / abs(REFERENCE[key]) for key in sep)),
    }
    print(json.dumps(gated, sort_keys=True))
    return gated


def stage_report():
    gated = stage_compare(write=True)
    files = ["requirements.lock", "draws_scaled.csv", "fig_spectra_scaled.png", "scripts/t01y_run.py"]
    command = 'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T01y/scripts/t01y_run.py'
    report = {
        "task_id": "T01y", "executor": "manus", "branch": "task/T01y",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": SEED, "started": "2026-10-15T04:49:54-06:00", "finished": datetime.now(timezone.utc).isoformat(), "status": "complete",
        "gated": gated, "deliverables": [{"path": path, "sha256": digest(REPORT / path)} for path in files],
        "regeneration": {**{path: command + " --stage all" for path in files}, **{key: command + " --stage compare" for key in gated}},
        "deviations": [], "deviations_requested": [],
        "observations": [f"Input ops_v2.npz SHA-256: {digest(input_path('v2'))}.", f"Input ops_v3.npz SHA-256: {digest(input_path('v3'))}."],
        "known_failures": KNOWN_FAILURES, "lockfile_sha256": digest(REPORT / "requirements.lock"),
    }
    (REPORT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("compare", "report", "all"))
    stage = parser.parse_args().stage
    if stage in ("compare", "all"):
        stage_compare()
    if stage in ("report", "all"):
        stage_report()


if __name__ == "__main__":
    main()
