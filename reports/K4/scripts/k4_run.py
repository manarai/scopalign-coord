#!/usr/bin/env python3
"""K4: label-blind Leiden anchors from cached centroid-operator spectra.

This runner only consumes the K3 and T04 caches named by tasks/K4.md.  It does
not fit drifts or calculate new per-cell operators.  All state it creates is
kept below reports/K4, and per-copy clustering checkpoints are deliberately
kept below its ignored cache/ directory for restart recovery.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scjdo2.compare import centroid_operator

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "K4"
CACHE = OUT / "cache"
K3 = ROOT / "reports" / "K3"
K3_CACHE = K3 / "cache"
T04 = ROOT / "reports" / "T04"
T04_CACHE = T04 / "cache"
K0 = ROOT / "reports" / "K0"
DATA = ROOT / "data" / "pbmc_10k_v3"
RAW = DATA / "pbmc_10k_v3.h5ad"
HALVES = DATA / "k0_halves.csv"
GENES_PATH = K0 / "genes.txt"

BASE_SEED = 20261014
NULL_SEED = 20261400
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
D = 30
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 1.0
MIN_CLUSTER_CELLS = 30
N_CENTROID = 30
N_PERMUTATIONS = 1000
PAIR_ORDER = ("P1", "P2", "P3")
COPY_ORDER = {"P1": ("H1", "H2"), "P2": ("A", "B"), "P3": ("A", "B")}
KNOWN_FAILURES = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
GATED_MATCH_KEYS = (
    "n_pairs_P1",
    "n_pairs_P2",
    "acc_op_P1",
    "acc_op_P2",
    "acc_op_P3",
    "acc_expr_P1",
    "acc_expr_P2",
    "acc_expr_P3",
    "chance_P2",
    "p_op_P2",
)
ALL_GATED_KEYS = (*GATED_MATCH_KEYS, "symmetry_max", "k4_outcome")


def utcnow() -> str:
    """Return an RFC 3339-compatible UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    """Calculate a content hash without loading a full input into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def is_below(path: Path, parent: Path) -> bool:
    """Return whether an absolute resolved path is contained by a parent."""
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        return False
    return True


def output_path(path: Path) -> Path:
    """Reject any write destination outside the one authorized report tree."""
    if not is_below(path, OUT):
        raise RuntimeError(f"K4 may only write below {OUT}: {path}")
    return path


def atomic_text(path: Path, value: str) -> None:
    """Atomically replace a text result after creating its K4-only parent."""
    target = output_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json(path: Path, value: Any) -> None:
    """Atomically serialize a JSON artifact with stable key ordering."""
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_csv(path: Path, table: pd.DataFrame, columns: list[str]) -> None:
    """Atomically write an exact-column CSV result."""
    target = output_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        table.to_csv(temporary, index=False, columns=columns, float_format="%.17g")
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    """Atomically write a compressed NPZ result below reports/K4."""
    target = output_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_figure(path: Path, figure: plt.Figure) -> None:
    """Atomically publish a PNG figure only after matplotlib has finished it."""
    target = output_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        figure.savefig(temporary, format="png", dpi=180)
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        plt.close(figure)
        if temporary.exists():
            temporary.unlink()


def ensure_layout() -> None:
    """Create only K4 paths and the required exclusions for transient state."""
    for directory in (OUT, CACHE, CACHE / "clusters"):
        output_path(directory).mkdir(parents=True, exist_ok=True)
    atomic_text(OUT / ".gitignore", "cache/\nrun.log\nrun.pid\n__pycache__/\n")


def require_input(path: Path) -> Path:
    """Ensure an authorized input exists without opening it for writing."""
    allowed_roots = (K3_CACHE, T04_CACHE, RAW.parent, K0)
    if not any(is_below(path, root) for root in allowed_roots):
        raise RuntimeError(f"Unexpected K4 input outside authorized locations: {path}")
    if not path.is_file():
        raise RuntimeError(f"Required K4 input is missing: {path}")
    return path


def p1_drift_path(copy_name: str, seed: int) -> Path:
    return K3_CACHE / "drifts" / f"P1_{copy_name.lower()}_{seed}.pt"


def p1_operator_path(copy_name: str, seed: int) -> Path:
    return K3_CACHE / "operators" / f"P1_{copy_name.lower()}_{seed}.npy"


def t04_drift_path(scenario: str, copy_name: str, seed: int) -> Path:
    return T04_CACHE / "drifts" / f"{scenario}_{copy_name}_{seed}.pt"


def t04_operator_path(scenario: str, copy_name: str, seed: int) -> Path:
    return T04_CACHE / "operators" / f"{scenario}_{copy_name}_{seed}.npy"


def p1_latents_path() -> Path:
    return K3_CACHE / "p1_latents.npz"


def t04_embedding_path(scenario: str) -> Path:
    return T04_CACHE / "embeddings" / f"{scenario}.npz"


def t04_data_path(scenario: str) -> Path:
    return T04_CACHE / "data" / f"{scenario}.npz"


def all_input_paths() -> list[Path]:
    """Return every cache and source-data file K4 consumes, in fixed order."""
    paths = [RAW, HALVES, GENES_PATH, p1_latents_path()]
    for copy_name in ("H1", "H2"):
        for seed in SEEDS:
            paths.append(p1_drift_path(copy_name, seed))
            paths.append(p1_operator_path(copy_name, seed))
    for pair, scenario in (("P2", "S0"), ("P3", "S1")):
        del pair
        paths.append(t04_embedding_path(scenario))
        paths.append(t04_data_path(scenario))
        for copy_name in ("A", "B"):
            for seed in SEEDS:
                paths.append(t04_drift_path(scenario, copy_name, seed))
                paths.append(t04_operator_path(scenario, copy_name, seed))
    return paths


def input_manifest_path() -> Path:
    return OUT / "input_manifest.json"


def write_input_manifest() -> None:
    """Record SHA-256s for every input cache and source-data file K4 uses."""
    entries = []
    for path in all_input_paths():
        source = require_input(path)
        entries.append(
            {
                "path": str(source.relative_to(ROOT)),
                "sha256": sha256(source),
                "bytes": int(source.stat().st_size),
            }
        )
    atomic_json(
        input_manifest_path(),
        {
            "inputs": entries,
            "read_only_directories": [
                str(K3_CACHE.relative_to(ROOT)),
                str(T04_CACHE.relative_to(ROOT)),
            ],
        },
    )


def manifest_hash() -> str:
    path = input_manifest_path()
    if not path.exists():
        raise RuntimeError("Run --stage check before a K4 computational stage")
    return sha256(path)


def gates_path() -> Path:
    return OUT / "gates.json"


def get_gates() -> dict[str, Any]:
    path = gates_path()
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def update_gates(**values: Any) -> dict[str, Any]:
    """Merge checkpointed values without discarding an earlier validity stop."""
    gates = get_gates()
    gates.update(values)
    atomic_json(gates_path(), gates)
    return gates


def stop_at(rule: str) -> None:
    """Persist the exact pre-registered stop rule once it has fired."""
    update_gates(stopped=True, stop_rule=rule)


def is_stopped() -> bool:
    return bool(get_gates().get("stopped", False))


def require_active() -> None:
    """Prevent later computation after K4's mandatory validity stop."""
    if is_stopped():
        raise RuntimeError("K4 has a pre-registered validity stop; only --stage report may run")


def requirements_path() -> Path:
    return OUT / "requirements.lock"


def write_requirements() -> None:
    """Capture the runner environment once, preserving it across resume calls."""
    if not requirements_path().exists():
        frozen = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
        atomic_text(requirements_path(), frozen)


def run_meta_path() -> Path:
    return OUT / "run_meta.json"


def make_meta() -> dict[str, Any]:
    """Create immutable descriptive run metadata on the first K4 stage."""
    path = run_meta_path()
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    value = {
        "started": utcnow(),
        "seed": BASE_SEED,
        "null_seed": NULL_SEED,
        "drift_seeds": list(SEEDS),
        "n_neighbors": N_NEIGHBORS,
        "leiden_resolution": LEIDEN_RESOLUTION,
        "leiden_flavor": "igraph",
        "leiden_iterations": 2,
        "minimum_cluster_cells": MIN_CLUSTER_CELLS,
        "operator_centroid_cells": N_CENTROID,
        "permutations": N_PERMUTATIONS,
        "permutation_pair_order": list(PAIR_ORDER),
    }
    atomic_json(path, value)
    return value


def load_npz(path: Path) -> dict[str, np.ndarray]:
    """Read a cache archive as ordinary arrays, never as writable mappings."""
    source = require_input(path)
    with np.load(source, allow_pickle=False) as arrays:
        return {key: np.asarray(arrays[key]) for key in arrays.files}


def require_keys(values: dict[str, np.ndarray], keys: tuple[str, ...], path: Path) -> None:
    missing = [key for key in keys if key not in values]
    if missing:
        raise RuntimeError(f"Missing keys in {path}: {missing}")


def validate_operator(path: Path, n_cells: int) -> None:
    """Read and validate a saved operator array without calculating one."""
    source = require_input(path)
    matrix = np.load(source, mmap_mode="r", allow_pickle=False)
    try:
        if matrix.shape != (n_cells, D, D) or matrix.dtype != np.float64:
            raise RuntimeError(f"Unexpected cached operator shape or dtype: {source}")
    finally:
        del matrix


def saved_symmetry(path: Path) -> tuple[float, float]:
    """Read only saved checkpoint symmetry, never instantiate or fit a model."""
    source = require_input(path)
    saved = torch.load(source, map_location="cpu", weights_only=False)
    if not isinstance(saved, dict) or "symmetry" not in saved:
        raise RuntimeError(f"Checkpoint lacks saved symmetry information: {source}")
    symmetry = saved["symmetry"]
    if not isinstance(symmetry, dict) or "median_asym" not in symmetry:
        raise RuntimeError(f"Checkpoint lacks saved median_asym: {source}")
    median = float(symmetry["median_asym"])
    maximum = float(symmetry.get("max_asym", float("nan")))
    if not np.isfinite(median):
        raise RuntimeError(f"Checkpoint saved a non-finite median_asym: {source}")
    return median, maximum


def load_p1_latents() -> dict[str, np.ndarray]:
    """Load K3's complete P1 latent archive and validate its aligned arrays."""
    path = p1_latents_path()
    values = load_npz(path)
    require_keys(values, ("z_h1", "z_h2", "w_h1", "w_h2", "types_h1", "types_h2"), path)
    for copy_name in ("h1", "h2"):
        latent = values[f"z_{copy_name}"]
        loading = values[f"w_{copy_name}"]
        labels = values[f"types_{copy_name}"]
        if latent.ndim != 2 or latent.shape[1] != D:
            raise RuntimeError(f"K3 P1 latent has unexpected shape for {copy_name}")
        if loading.shape != (D, 2000) or len(labels) != len(latent):
            raise RuntimeError(f"K3 P1 cached arrays do not align for {copy_name}")
    return values


def load_t04_case(scenario: str) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Load one T04 latent/count cache pair and validate A/B row alignment."""
    embedding_path = t04_embedding_path(scenario)
    data_path = t04_data_path(scenario)
    embedding = load_npz(embedding_path)
    data = load_npz(data_path)
    require_keys(embedding, ("z_a", "z_b", "w_a", "w_b", "types"), embedding_path)
    require_keys(data, ("A", "B", "types_a", "types_b"), data_path)
    n_a = len(embedding["z_a"])
    n_b = len(embedding["z_b"])
    if embedding["z_a"].shape != (n_a, D) or embedding["z_b"].shape != (n_b, D):
        raise RuntimeError(f"T04 {scenario} latents must each be 30-dimensional")
    if embedding["w_a"].shape != (D, 2000) or embedding["w_b"].shape != (D, 2000):
        raise RuntimeError(f"T04 {scenario} PCA loading shapes are invalid")
    if len(embedding["types"]) != n_a + n_b:
        raise RuntimeError(f"T04 {scenario} combined labels do not align with latents")
    if data["A"].shape != (n_a, 2000) or data["B"].shape != (n_b, 2000):
        raise RuntimeError(f"T04 {scenario} cached count dimensions do not align")
    if not np.array_equal(data["types_a"].astype(str), embedding["types"][:n_a].astype(str)):
        raise RuntimeError(f"T04 {scenario} A labels differ between data and embedding caches")
    if not np.array_equal(data["types_b"].astype(str), embedding["types"][n_a:].astype(str)):
        raise RuntimeError(f"T04 {scenario} B labels differ between data and embedding caches")
    return embedding, data


def validate_all_inputs() -> list[dict[str, Any]]:
    """Open all named latent, fit, and operator caches without deriving anything."""
    p1 = load_p1_latents()
    cases = {"S0": load_t04_case("S0"), "S1": load_t04_case("S1")}
    rows: list[dict[str, Any]] = []
    for copy_name in ("H1", "H2"):
        n_cells = len(p1[f"z_{copy_name.lower()}"])
        for seed in SEEDS:
            drift = p1_drift_path(copy_name, seed)
            operator = p1_operator_path(copy_name, seed)
            median, maximum = saved_symmetry(drift)
            validate_operator(operator, n_cells)
            rows.append(
                {
                    "pair": "P1",
                    "copy": copy_name,
                    "seed": int(seed),
                    "drift_path": str(drift.relative_to(ROOT)),
                    "operator_path": str(operator.relative_to(ROOT)),
                    "median_asym": median,
                    "saved_max_asym": maximum,
                }
            )
    for pair, scenario in (("P2", "S0"), ("P3", "S1")):
        embedding, data = cases[scenario]
        del data
        for copy_name in ("A", "B"):
            n_cells = len(embedding[f"z_{copy_name.lower()}"])
            for seed in SEEDS:
                drift = t04_drift_path(scenario, copy_name, seed)
                operator = t04_operator_path(scenario, copy_name, seed)
                median, maximum = saved_symmetry(drift)
                validate_operator(operator, n_cells)
                rows.append(
                    {
                        "pair": pair,
                        "copy": copy_name,
                        "seed": int(seed),
                        "drift_path": str(drift.relative_to(ROOT)),
                        "operator_path": str(operator.relative_to(ROOT)),
                        "median_asym": median,
                        "saved_max_asym": maximum,
                    }
                )
    return rows


def check_stage() -> None:
    """Hash read-only inputs and gate saved checkpoint symmetry without fitting."""
    ensure_layout()
    make_meta()
    write_requirements()
    write_input_manifest()
    rows = validate_all_inputs()
    columns = ["pair", "copy", "seed", "drift_path", "operator_path", "median_asym", "saved_max_asym"]
    atomic_csv(OUT / "symmetry_stats.csv", pd.DataFrame(rows), columns)
    symmetry_max = float(max(float(row["median_asym"]) for row in rows))
    update_gates(symmetry_max=symmetry_max)
    if symmetry_max > 1e-4:
        stop_at("symmetry_max > 1e-4 → REVISE.")
    print(json.dumps({"symmetry_max": symmetry_max}, sort_keys=True), flush=True)


def pair_latents_and_labels(pair: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
    """Return the two specified latent copies and their coarse labels."""
    if pair == "P1":
        values = load_p1_latents()
        return (
            np.asarray(values["z_h1"], dtype=np.float64),
            np.asarray(values["types_h1"].astype(str)),
            np.asarray(values["z_h2"], dtype=np.float64),
            np.asarray(values["types_h2"].astype(str)),
            sha256(p1_latents_path()),
        )
    scenario = "S0" if pair == "P2" else "S1"
    embedding, data = load_t04_case(scenario)
    return (
        np.asarray(embedding["z_a"], dtype=np.float64),
        np.asarray(data["types_a"].astype(str)),
        np.asarray(embedding["z_b"], dtype=np.float64),
        np.asarray(data["types_b"].astype(str)),
        sha256(t04_embedding_path(scenario)),
    )


def cluster_cache_path(pair: str, copy_name: str) -> Path:
    return CACHE / "clusters" / f"{pair}_{copy_name}.npz"


def cluster_sort_key(item: tuple[str, int]) -> tuple[int, str]:
    """Order raw Leiden IDs by decreasing size and alphabetically on ties."""
    return (-int(item[1]), str(item[0]))


def majority_label(labels: np.ndarray) -> tuple[str, float]:
    """Return the alphabetically first majority label and its cluster purity."""
    values, counts = np.unique(np.asarray(labels, dtype=str), return_counts=True)
    winner_count = int(counts.max())
    winners = values[counts == winner_count]
    winner = sorted(winners.astype(str).tolist())[0]
    return winner, float(winner_count / len(labels))


def validate_cluster_cache(
    path: Path, pair: str, copy_name: str, n_cells: int, source_hash: str
) -> bool:
    """Check whether a saved per-copy clustering checkpoint can be resumed."""
    if not path.exists():
        return False
    try:
        with np.load(path, allow_pickle=False) as values:
            required = {
                "version", "pair", "copy", "source_sha256", "n_cells", "cell_cluster",
                "cluster_id", "cluster_n_cells", "cluster_label", "cluster_label_purity",
                "dropped_cluster", "dropped_n_cells",
            }
            if not required.issubset(set(values.files)):
                return False
            if str(values["version"].item()) != "K4_cluster_v1":
                return False
            if str(values["pair"].item()) != pair or str(values["copy"].item()) != copy_name:
                return False
            if str(values["source_sha256"].item()) != source_hash or int(values["n_cells"].item()) != n_cells:
                return False
            if len(values["cell_cluster"]) != n_cells:
                return False
            n_clusters = len(values["cluster_id"])
            if not all(
                len(values[key]) == n_clusters
                for key in ("cluster_n_cells", "cluster_label", "cluster_label_purity")
            ):
                return False
    except (OSError, ValueError, KeyError):
        return False
    return True


def leiden_cluster(
    pair: str, copy_name: str, latent: np.ndarray, labels: np.ndarray, source_hash: str
) -> None:
    """Run the task-specified latent AnnData neighbor graph and Leiden clustering."""
    target = cluster_cache_path(pair, copy_name)
    if validate_cluster_cache(target, pair, copy_name, len(latent), source_hash):
        return
    work = ad.AnnData(X=np.asarray(latent, dtype=np.float64).copy())
    sc.pp.neighbors(work, n_neighbors=15, use_rep="X", random_state=20261014)
    sc.tl.leiden(
        work,
        resolution=1.0,
        random_state=20261014,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    raw = work.obs["leiden"].astype(str).to_numpy()
    raw_ids, raw_counts = np.unique(raw, return_counts=True)
    ranked = sorted(zip(raw_ids.astype(str).tolist(), raw_counts.astype(int).tolist()), key=cluster_sort_key)
    keep = [(raw_id, n_cells) for raw_id, n_cells in ranked if n_cells >= MIN_CLUSTER_CELLS]
    dropped = [(raw_id, n_cells) for raw_id, n_cells in ranked if n_cells < MIN_CLUSTER_CELLS]
    assignments = np.full(len(raw), "", dtype="<U16")
    cluster_ids: list[str] = []
    cluster_sizes: list[int] = []
    cluster_labels: list[str] = []
    cluster_purities: list[float] = []
    for rank, (raw_id, n_cells) in enumerate(keep, start=1):
        cluster_id = f"C{rank:02d}"
        member = raw == raw_id
        label, purity = majority_label(labels[member])
        assignments[member] = cluster_id
        cluster_ids.append(cluster_id)
        cluster_sizes.append(int(n_cells))
        cluster_labels.append(label)
        cluster_purities.append(purity)
    atomic_npz(
        target,
        version=np.asarray("K4_cluster_v1"),
        pair=np.asarray(pair),
        copy=np.asarray(copy_name),
        source_sha256=np.asarray(source_hash),
        n_cells=np.asarray(len(latent), dtype=np.int64),
        cell_cluster=assignments,
        cluster_id=np.asarray(cluster_ids, dtype=str),
        cluster_n_cells=np.asarray(cluster_sizes, dtype=np.int64),
        cluster_label=np.asarray(cluster_labels, dtype=str),
        cluster_label_purity=np.asarray(cluster_purities, dtype=np.float64),
        dropped_cluster=np.asarray([item[0] for item in dropped], dtype=str),
        dropped_n_cells=np.asarray([item[1] for item in dropped], dtype=np.int64),
    )


def load_cluster_cache(pair: str, copy_name: str) -> dict[str, np.ndarray]:
    """Read K4's own completed clustering checkpoint as immutable arrays."""
    path = cluster_cache_path(pair, copy_name)
    if not path.exists():
        raise RuntimeError(f"Cluster cache is missing: run --stage cluster ({path})")
    with np.load(path, allow_pickle=False) as values:
        return {key: np.asarray(values[key]) for key in values.files}


def cluster_rows_and_drops() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build the committed cluster CSV and small dropped-count artifact from cache."""
    rows: list[dict[str, Any]] = []
    dropped: dict[str, Any] = {"minimum_cluster_cells": MIN_CLUSTER_CELLS, "pairs": {}}
    for pair in PAIR_ORDER:
        dropped["pairs"][pair] = {}
        for copy_name in COPY_ORDER[pair]:
            values = load_cluster_cache(pair, copy_name)
            ids = values["cluster_id"].astype(str)
            sizes = values["cluster_n_cells"].astype(int)
            labels = values["cluster_label"].astype(str)
            purities = values["cluster_label_purity"].astype(float)
            for cluster_id, n_cells, label, purity in zip(ids, sizes, labels, purities):
                rows.append(
                    {
                        "pair": pair,
                        "copy": copy_name,
                        "cluster": cluster_id,
                        "n_cells": int(n_cells),
                        "label": label,
                        "label_purity": float(purity),
                    }
                )
            dropped_sizes = values["dropped_n_cells"].astype(int)
            dropped["pairs"][pair][copy_name] = {
                "n_input_cells": int(values["n_cells"].item()),
                "n_input_clusters": int(len(ids) + len(dropped_sizes)),
                "n_kept_clusters": int(len(ids)),
                "n_kept_cells": int(sizes.sum()),
                "n_dropped_clusters": int(len(dropped_sizes)),
                "n_dropped_cells": int(dropped_sizes.sum()),
            }
    return rows, dropped


def cluster_stage() -> None:
    """Cluster each copy separately with the exact specified Leiden invocation."""
    ensure_layout()
    require_active()
    manifest_hash()
    for pair in PAIR_ORDER:
        z1, labels1, z2, labels2, source_hash = pair_latents_and_labels(pair)
        copy1, copy2 = COPY_ORDER[pair]
        leiden_cluster(pair, copy1, z1, labels1, source_hash)
        leiden_cluster(pair, copy2, z2, labels2, source_hash)
    rows, dropped = cluster_rows_and_drops()
    columns = ["pair", "copy", "cluster", "n_cells", "label", "label_purity"]
    atomic_csv(OUT / "clusters.csv", pd.DataFrame(rows), columns)
    atomic_json(OUT / "cluster_drops.json", dropped)
    print(json.dumps({"clusters": int(len(rows))}, sort_keys=True), flush=True)


def copy_cluster_metadata(pair: str, copy_name: str) -> list[dict[str, Any]]:
    """Return kept clusters in their mandated descending-size rank order."""
    values = load_cluster_cache(pair, copy_name)
    result = []
    for cluster_id, n_cells, label, purity in zip(
        values["cluster_id"].astype(str),
        values["cluster_n_cells"].astype(int),
        values["cluster_label"].astype(str),
        values["cluster_label_purity"].astype(float),
    ):
        result.append(
            {
                "cluster": str(cluster_id),
                "n_cells": int(n_cells),
                "label": str(label),
                "label_purity": float(purity),
            }
        )
    return result


def copy_cluster_indices(pair: str, copy_name: str, clusters: list[dict[str, Any]]) -> list[np.ndarray]:
    """Map kept rank-order cluster IDs to their original latent row indices."""
    values = load_cluster_cache(pair, copy_name)
    assignments = values["cell_cluster"].astype(str)
    indices = []
    for cluster in clusters:
        member = np.flatnonzero(assignments == cluster["cluster"])
        if len(member) != int(cluster["n_cells"]):
            raise RuntimeError(f"Cluster cache membership count mismatch for {pair}/{copy_name}/{cluster['cluster']}")
        indices.append(member)
    return indices


def operators_for_pair(pair: str, copy_name: str, n_cells: int) -> dict[int, np.ndarray]:
    """Open the five saved float64 operator arrays for one named copy."""
    result: dict[int, np.ndarray] = {}
    if pair == "P1":
        path_for = p1_operator_path
    else:
        scenario = "S0" if pair == "P2" else "S1"

        def path_for(local_copy: str, local_seed: int) -> Path:
            return t04_operator_path(scenario, local_copy, local_seed)

    for seed in SEEDS:
        path = path_for(copy_name, seed)
        validate_operator(path, n_cells)
        result[seed] = np.load(require_input(path), mmap_mode="r", allow_pickle=False)
    return result


def operator_fingerprints(
    latent: np.ndarray, matrices: dict[int, np.ndarray], cluster_indices: list[np.ndarray]
) -> np.ndarray:
    """Average the five saved centroid-operator eigenspectra for each cluster."""
    fingerprints = []
    z = np.asarray(latent, dtype=np.float64)
    for indices in cluster_indices:
        spectra = []
        for seed in SEEDS:
            centroid = centroid_operator(
                z[indices],
                np.asarray(matrices[seed][indices], dtype=np.float64),
                n_cells_centroid=30,
            ).astype(np.float64)
            symmetric = (centroid + centroid.T) / 2.0
            spectra.append(np.linalg.eigvalsh(symmetric)[::-1])
        fingerprints.append(np.mean(np.vstack(spectra), axis=0, dtype=np.float64))
    return np.vstack(fingerprints).astype(np.float64)


def standardised_operator_distance(left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pool-standardize 30 spectral coordinates and calculate Euclidean distance."""
    pooled = np.vstack([left, right]).astype(np.float64)
    mean = pooled.mean(axis=0, dtype=np.float64)
    scale = pooled.std(axis=0, dtype=np.float64)
    if np.any(scale == 0):
        raise RuntimeError("A pooled operator fingerprint coordinate has zero standard deviation")
    left_standardised = (left - mean) / scale
    right_standardised = (right - mean) / scale
    return cdist(left_standardised, right_standardised, metric="euclidean"), mean, scale, np.vstack([left_standardised, right_standardised])


def read_sorted_genes() -> list[str]:
    """Load K0's task-mandated sorted 2,000-gene list."""
    source = require_input(GENES_PATH)
    genes = [line.strip() for line in source.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(genes) != 2000 or genes != sorted(genes):
        raise RuntimeError("K0 genes.txt must contain exactly 2,000 alphabetically sorted genes")
    return genes


def p1_raw_halves() -> tuple[Any, Any]:
    """Reconstruct K3's exact raw-QC and K0 half split for the expression baseline."""
    source = sc.read_h5ad(require_input(RAW))
    source = source[source.obs["qc_pass"].to_numpy(dtype=bool)].copy()
    assignments = pd.read_csv(require_input(HALVES)).set_index("barcode")["half"]
    names = source.obs_names.astype(str)
    mapped = names.map(assignments)
    h1 = source[(names.isin(assignments.index)) & (mapped == "H1")].copy()
    h2 = source[(names.isin(assignments.index)) & (mapped == "H2")].copy()
    if h1.n_obs == 0 or h2.n_obs == 0:
        raise RuntimeError("K0 H1/H2 raw split could not be reconstructed")
    p1 = load_p1_latents()
    if not np.array_equal(h1.obs["cell_type"].astype(str).to_numpy(), p1["types_h1"].astype(str)):
        raise RuntimeError("P1 H1 raw cell-type order differs from K3's latent cache")
    if not np.array_equal(h2.obs["cell_type"].astype(str).to_numpy(), p1["types_h2"].astype(str)):
        raise RuntimeError("P1 H2 raw cell-type order differs from K3's latent cache")
    return h1, h2


def normalised_adata_matrix(adata: Any, genes: list[str]) -> np.ndarray:
    """Apply K3's own normalize_total(1e4), log1p baseline to selected raw counts."""
    work = adata[:, genes].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    matrix = work.X.toarray() if hasattr(work.X, "toarray") else np.asarray(work.X)
    return np.asarray(matrix)


def normalised_cached_counts(counts: np.ndarray) -> np.ndarray:
    """Apply T04's exact count-cache normalization route to one copy."""
    work = ad.AnnData(X=sp.csr_matrix(np.asarray(counts, dtype=np.float32)))
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    matrix = work.X.toarray() if sp.issparse(work.X) else np.asarray(work.X)
    return np.asarray(matrix)


def expression_for_pair(pair: str) -> tuple[np.ndarray, np.ndarray]:
    """Load each copy's own required 2,000-gene log-normalized baseline expression."""
    if pair == "P1":
        genes = read_sorted_genes()
        h1, h2 = p1_raw_halves()
        return normalised_adata_matrix(h1, genes), normalised_adata_matrix(h2, genes)
    scenario = "S0" if pair == "P2" else "S1"
    _, data = load_t04_case(scenario)
    return normalised_cached_counts(data["A"]), normalised_cached_counts(data["B"])


def expression_fingerprints(expression: np.ndarray, cluster_indices: list[np.ndarray]) -> np.ndarray:
    """Calculate one mean log-normalized expression vector for every cluster."""
    return np.vstack(
        [np.asarray(expression[indices], dtype=np.float64).mean(axis=0, dtype=np.float64) for indices in cluster_indices]
    )


def pearson_distance_matrix(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Calculate 1 - Pearson r for each pair of cluster mean expression vectors."""
    distance = np.empty((len(left), len(right)), dtype=np.float64)
    for row_index, row in enumerate(np.asarray(left, dtype=np.float64)):
        centered_row = row - row.mean(dtype=np.float64)
        row_norm = float(np.sqrt(np.sum(centered_row * centered_row)))
        if row_norm == 0:
            raise RuntimeError("A cluster expression fingerprint has zero variance")
        for col_index, column in enumerate(np.asarray(right, dtype=np.float64)):
            centered_column = column - column.mean(dtype=np.float64)
            column_norm = float(np.sqrt(np.sum(centered_column * centered_column)))
            if column_norm == 0:
                raise RuntimeError("A cluster expression fingerprint has zero variance")
            correlation = float(np.dot(centered_row, centered_column) / (row_norm * column_norm))
            distance[row_index, col_index] = 1.0 - correlation
    return distance


def matching(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Run the required rectangular Hungarian matching with copy 1 as rows."""
    if distance.ndim != 2 or min(distance.shape) == 0:
        return np.asarray([], dtype=int), np.asarray([], dtype=int)
    rows, columns = linear_sum_assignment(distance)
    return np.asarray(rows, dtype=int), np.asarray(columns, dtype=int)


def match_rows_for_fingerprint(
    pair: str,
    fingerprint: str,
    distance: np.ndarray,
    row_clusters: list[dict[str, Any]],
    col_clusters: list[dict[str, Any]],
    rows: np.ndarray,
    columns: np.ndarray,
) -> list[dict[str, Any]]:
    """Format chosen matching edges in the exact task-specified CSV schema."""
    output = []
    for row_index, col_index in zip(rows, columns):
        row_cluster = row_clusters[int(row_index)]
        col_cluster = col_clusters[int(col_index)]
        output.append(
            {
                "pair": pair,
                "fingerprint": fingerprint,
                "row_cluster": row_cluster["cluster"],
                "col_cluster": col_cluster["cluster"],
                "distance": float(distance[int(row_index), int(col_index)]),
                "row_label": row_cluster["label"],
                "col_label": col_cluster["label"],
                "correct": bool(row_cluster["label"] == col_cluster["label"]),
            }
        )
    return output


def accuracy(rows: list[dict[str, Any]]) -> float:
    """Return label-only accuracy after an already label-blind matching."""
    if not rows:
        return float("nan")
    return float(sum(bool(row["correct"]) for row in rows) / len(rows))


def permutation_accuracy(labels1: list[str], labels2: list[str], rng: np.random.Generator) -> float:
    """Draw exactly one direction-prescribed random rectangular matching accuracy."""
    n1, n2 = len(labels1), len(labels2)
    if n1 <= n2:
        selected = rng.permutation(n2)[:n1]
        correct = sum(labels1[index] == labels2[int(selected[index])] for index in range(n1))
        return float(correct / n1) if n1 else float("nan")
    selected = rng.permutation(n1)[:n2]
    correct = sum(labels1[int(selected[index])] == labels2[index] for index in range(n2))
    return float(correct / n2) if n2 else float("nan")


def chance_null(pair_results: dict[str, dict[str, Any]]) -> None:
    """Apply the one shared 1,000-draw null in prescribed P1/P2/P3 order."""
    rng = np.random.default_rng(NULL_SEED)
    for pair in PAIR_ORDER:
        result = pair_results[pair]
        labels1 = [str(cluster["label"]) for cluster in result["row_clusters"]]
        labels2 = [str(cluster["label"]) for cluster in result["col_clusters"]]
        values = [permutation_accuracy(labels1, labels2, rng) for _ in range(N_PERMUTATIONS)]
        null = np.asarray(values, dtype=np.float64)
        result["chance_null"] = null
        result["chance"] = float(null.mean(dtype=np.float64))
        result["p_op"] = float((1 + int(np.sum(null >= float(result["acc_op"])))) / (N_PERMUTATIONS + 1))


def fingerprints_path() -> Path:
    return OUT / "fingerprints.npz"


def scores_path() -> Path:
    return OUT / "scores.json"


def match_state_path() -> Path:
    return CACHE / "match_state.json"


def cluster_cache_hashes() -> dict[str, str]:
    """Hash all six cluster checkpoints so a resumed match cannot use stale groups."""
    result: dict[str, str] = {}
    for pair in PAIR_ORDER:
        for copy_name in COPY_ORDER[pair]:
            path = cluster_cache_path(pair, copy_name)
            if not path.exists():
                raise RuntimeError("Run --stage cluster before --stage match")
            result[f"{pair}_{copy_name}"] = sha256(path)
    return result


def valid_match_resume() -> dict[str, Any] | None:
    """Return saved gated values only when every persisted match input still agrees."""
    path = match_state_path()
    required_outputs = (OUT / "matches.csv", fingerprints_path(), scores_path())
    if not path.exists() or not all(candidate.exists() for candidate in required_outputs):
        return None
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
        if state.get("version") != "K4_match_v1":
            return None
        if state.get("input_manifest_sha256") != manifest_hash():
            return None
        if state.get("cluster_cache_sha256") != cluster_cache_hashes():
            return None
        gated = state.get("gated")
        if not isinstance(gated, dict) or set(gated) != set(GATED_MATCH_KEYS):
            return None
        return gated
    except (OSError, ValueError, KeyError):
        return None


def match_stage() -> None:
    """Construct operator/expression fingerprints, match them, and score labels afterward."""
    ensure_layout()
    require_active()
    resumed = valid_match_resume()
    if resumed is not None:
        update_gates(**resumed)
        if int(resumed["n_pairs_P1"]) < 8 or int(resumed["n_pairs_P2"]) < 8:
            stop_at("n_pairs_P1 < 8 or n_pairs_P2 < 8 → REVISE.")
        print(json.dumps(resumed, sort_keys=True), flush=True)
        return

    pair_results: dict[str, dict[str, Any]] = {}
    all_rows: list[dict[str, Any]] = []
    fingerprint_arrays: dict[str, np.ndarray] = {}
    for pair in PAIR_ORDER:
        z1, labels1, z2, labels2, source_hash = pair_latents_and_labels(pair)
        del labels1, labels2, source_hash
        copy1, copy2 = COPY_ORDER[pair]
        row_clusters = copy_cluster_metadata(pair, copy1)
        col_clusters = copy_cluster_metadata(pair, copy2)
        row_indices = copy_cluster_indices(pair, copy1, row_clusters)
        col_indices = copy_cluster_indices(pair, copy2, col_clusters)
        matrices1 = operators_for_pair(pair, copy1, len(z1))
        matrices2 = operators_for_pair(pair, copy2, len(z2))
        op1 = operator_fingerprints(z1, matrices1, row_indices)
        op2 = operator_fingerprints(z2, matrices2, col_indices)
        op_distance, op_mean, op_scale, op_standardised = standardised_operator_distance(op1, op2)
        op_rows, op_columns = matching(op_distance)
        chosen_operator = match_rows_for_fingerprint(
            pair, "operator", op_distance, row_clusters, col_clusters, op_rows, op_columns
        )

        expression1, expression2 = expression_for_pair(pair)
        exp1 = expression_fingerprints(expression1, row_indices)
        exp2 = expression_fingerprints(expression2, col_indices)
        expression_distance = pearson_distance_matrix(exp1, exp2)
        expression_rows, expression_columns = matching(expression_distance)
        chosen_expression = match_rows_for_fingerprint(
            pair, "expression", expression_distance, row_clusters, col_clusters, expression_rows, expression_columns
        )
        all_rows.extend(chosen_operator)
        all_rows.extend(chosen_expression)
        pair_results[pair] = {
            "row_clusters": row_clusters,
            "col_clusters": col_clusters,
            "operator_distance": op_distance,
            "expression_distance": expression_distance,
            "operator_rows": op_rows,
            "operator_columns": op_columns,
            "expression_rows": expression_rows,
            "expression_columns": expression_columns,
            "acc_op": accuracy(chosen_operator),
            "acc_expr": accuracy(chosen_expression),
        }
        fingerprint_arrays[f"{pair}_{copy1}_cluster"] = np.asarray([item["cluster"] for item in row_clusters], dtype=str)
        fingerprint_arrays[f"{pair}_{copy2}_cluster"] = np.asarray([item["cluster"] for item in col_clusters], dtype=str)
        fingerprint_arrays[f"{pair}_{copy1}_label"] = np.asarray([item["label"] for item in row_clusters], dtype=str)
        fingerprint_arrays[f"{pair}_{copy2}_label"] = np.asarray([item["label"] for item in col_clusters], dtype=str)
        fingerprint_arrays[f"{pair}_{copy1}_operator_fingerprint"] = op1
        fingerprint_arrays[f"{pair}_{copy2}_operator_fingerprint"] = op2
        fingerprint_arrays[f"{pair}_operator_pooled_mean"] = op_mean
        fingerprint_arrays[f"{pair}_operator_pooled_sd"] = op_scale
        fingerprint_arrays[f"{pair}_{copy1}_operator_standardised"] = op_standardised[: len(op1)]
        fingerprint_arrays[f"{pair}_{copy2}_operator_standardised"] = op_standardised[len(op1) :]
        fingerprint_arrays[f"{pair}_operator_distance"] = op_distance
        fingerprint_arrays[f"{pair}_expression_distance"] = expression_distance
        fingerprint_arrays[f"{pair}_operator_match_rows"] = op_rows
        fingerprint_arrays[f"{pair}_operator_match_columns"] = op_columns
        fingerprint_arrays[f"{pair}_expression_match_rows"] = expression_rows
        fingerprint_arrays[f"{pair}_expression_match_columns"] = expression_columns
        for matrices in (matrices1, matrices2):
            matrices.clear()

    chance_null(pair_results)
    columns = ["pair", "fingerprint", "row_cluster", "col_cluster", "distance", "row_label", "col_label", "correct"]
    atomic_csv(OUT / "matches.csv", pd.DataFrame(all_rows), columns)
    for pair in PAIR_ORDER:
        fingerprint_arrays[f"{pair}_chance_null"] = np.asarray(pair_results[pair]["chance_null"], dtype=np.float64)
    atomic_npz(fingerprints_path(), **fingerprint_arrays)

    scores: dict[str, Any] = {
        "null_seed": NULL_SEED,
        "permutations": N_PERMUTATIONS,
        "pair_order": list(PAIR_ORDER),
        "pairs": {},
    }
    for pair in PAIR_ORDER:
        result = pair_results[pair]
        scores["pairs"][pair] = {
            "n_row_clusters": int(len(result["row_clusters"])),
            "n_col_clusters": int(len(result["col_clusters"])),
            "n_pairs": int(len(result["operator_rows"])),
            "acc_op": float(result["acc_op"]),
            "acc_expr": float(result["acc_expr"]),
            "chance": float(result["chance"]),
            "p_op": float(result["p_op"]),
        }
    atomic_json(scores_path(), scores)
    gated = {
        "n_pairs_P1": int(scores["pairs"]["P1"]["n_pairs"]),
        "n_pairs_P2": int(scores["pairs"]["P2"]["n_pairs"]),
        "acc_op_P1": float(scores["pairs"]["P1"]["acc_op"]),
        "acc_op_P2": float(scores["pairs"]["P2"]["acc_op"]),
        "acc_op_P3": float(scores["pairs"]["P3"]["acc_op"]),
        "acc_expr_P1": float(scores["pairs"]["P1"]["acc_expr"]),
        "acc_expr_P2": float(scores["pairs"]["P2"]["acc_expr"]),
        "acc_expr_P3": float(scores["pairs"]["P3"]["acc_expr"]),
        "chance_P2": float(scores["pairs"]["P2"]["chance"]),
        "p_op_P2": float(scores["pairs"]["P2"]["p_op"]),
    }
    atomic_json(
        match_state_path(),
        {
            "version": "K4_match_v1",
            "input_manifest_sha256": manifest_hash(),
            "cluster_cache_sha256": cluster_cache_hashes(),
            "gated": gated,
        },
    )
    update_gates(**gated)
    if gated["n_pairs_P1"] < 8 or gated["n_pairs_P2"] < 8:
        stop_at("n_pairs_P1 < 8 or n_pairs_P2 < 8 → REVISE.")
    print(json.dumps(gated, sort_keys=True), flush=True)


def label_cluster_sort_key(item: tuple[int, dict[str, Any]]) -> tuple[str, str]:
    """Sort plotting rows/columns by label, using rank cluster ID as a tiebreaker."""
    index, cluster = item
    del index
    return str(cluster["label"]), str(cluster["cluster"])


def plot_anchor_matrix(axis: plt.Axes, pair: str, values: dict[str, np.ndarray]) -> None:
    """Plot one label-ordered operator distance matrix with selected edges as pluses."""
    copy1, copy2 = COPY_ORDER[pair]
    row_clusters = [
        {"cluster": cluster, "label": label}
        for cluster, label in zip(values[f"{pair}_{copy1}_cluster"].astype(str), values[f"{pair}_{copy1}_label"].astype(str))
    ]
    col_clusters = [
        {"cluster": cluster, "label": label}
        for cluster, label in zip(values[f"{pair}_{copy2}_cluster"].astype(str), values[f"{pair}_{copy2}_label"].astype(str))
    ]
    row_order = [index for index, item in sorted(enumerate(row_clusters), key=label_cluster_sort_key)]
    col_order = [index for index, item in sorted(enumerate(col_clusters), key=label_cluster_sort_key)]
    distance = values[f"{pair}_operator_distance"]
    ordered = distance[np.ix_(row_order, col_order)]
    image = axis.imshow(ordered, aspect="auto", cmap="viridis")
    row_position = {original: position for position, original in enumerate(row_order)}
    col_position = {original: position for position, original in enumerate(col_order)}
    for row_index, col_index in zip(
        values[f"{pair}_operator_match_rows"].astype(int),
        values[f"{pair}_operator_match_columns"].astype(int),
    ):
        axis.scatter(
            col_position[int(col_index)],
            row_position[int(row_index)],
            marker="+",
            s=150,
            linewidths=1.8,
            color="white",
        )
    axis.set_xticks(np.arange(len(col_order)))
    axis.set_yticks(np.arange(len(row_order)))
    axis.set_xticklabels(
        [f"{col_clusters[index]['label']}\n{col_clusters[index]['cluster']}" for index in col_order],
        rotation=60,
        ha="right",
        fontsize=7,
    )
    axis.set_yticklabels(
        [f"{row_clusters[index]['label']}\n{row_clusters[index]['cluster']}" for index in row_order],
        fontsize=7,
    )
    axis.set_xlabel(f"{copy2} clusters (label order)")
    axis.set_ylabel(f"{copy1} clusters (label order)")
    axis.set_title(f"{pair}: operator distance")
    axis.figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04, label="Euclidean distance")


def anchors_figure() -> None:
    """Create the sole requested K4 figure from saved fingerprints and scores."""
    if not fingerprints_path().exists() or not scores_path().exists():
        raise RuntimeError("Run --stage match before --stage report")
    with np.load(fingerprints_path(), allow_pickle=False) as archive:
        values = {key: np.asarray(archive[key]) for key in archive.files}
    scores = json.loads(scores_path().read_text(encoding="utf-8"))
    figure, axes = plt.subplots(2, 2, figsize=(16, 13))
    for pair, axis in zip(PAIR_ORDER, axes.flat[:3]):
        plot_anchor_matrix(axis, pair, values)
    bar_axis = axes.flat[3]
    positions = np.arange(len(PAIR_ORDER))
    width = 0.24
    operator = [float(scores["pairs"][pair]["acc_op"]) for pair in PAIR_ORDER]
    expression = [float(scores["pairs"][pair]["acc_expr"]) for pair in PAIR_ORDER]
    chance = [float(scores["pairs"][pair]["chance"]) for pair in PAIR_ORDER]
    bar_axis.bar(positions - width, operator, width, label="acc_op", color="C3")
    bar_axis.bar(positions, expression, width, label="acc_expr", color="C0")
    bar_axis.bar(positions + width, chance, width, label="chance", color="0.55")
    bar_axis.set_xticks(positions, PAIR_ORDER)
    bar_axis.set_ylim(0, 1.05)
    bar_axis.set_ylabel("label-match accuracy")
    bar_axis.set_title("Anchoring accuracy and shared chance null")
    bar_axis.grid(axis="y", alpha=0.25)
    bar_axis.legend(frameon=False)
    figure.suptitle("K4 — Anchoring: operator fingerprints versus expression", y=0.995)
    figure.tight_layout()
    atomic_figure(OUT / "fig_anchors.png", figure)


def decide(gated: dict[str, Any]) -> str:
    """Implement the pre-registered outcome decision rule in its stated order."""
    if float(gated["p_op_P2"]) > 0.01 or float(gated["acc_op_P2"]) < float(gated["chance_P2"]) + 0.20:
        return "no_anchor"
    if (
        float(gated["acc_op_P1"]) >= 0.80
        and float(gated["acc_op_P2"]) >= 0.80
        and float(gated["acc_op_P2"]) >= float(gated["acc_expr_P2"]) + 0.10
    ):
        return "anchors_better"
    if float(gated["acc_op_P1"]) >= 0.80 and float(gated["acc_op_P2"]) >= 0.80:
        return "anchors_equal"
    return "weak"


def artifact_paths() -> list[Path]:
    """Return every K4 result artifact other than hashable declared deliverables."""
    paths = [
        input_manifest_path(),
        run_meta_path(),
        OUT / "symmetry_stats.csv",
        OUT / "cluster_drops.json",
        scores_path(),
        gates_path(),
        fingerprints_path(),
        match_state_path(),
    ]
    for pair in PAIR_ORDER:
        for copy_name in COPY_ORDER[pair]:
            paths.append(cluster_cache_path(pair, copy_name))
    return [path for path in paths if path.exists()]


def file_entry(path: Path) -> dict[str, str]:
    """Format one report-folder-relative checksum entry."""
    return {"path": str(path.relative_to(OUT)), "sha256": sha256(path)}


def complete_deliverable_paths() -> list[Path]:
    """List all hashable non-self-referential required deliverable files."""
    return [
        requirements_path(),
        OUT / "clusters.csv",
        OUT / "matches.csv",
        OUT / "fig_anchors.png",
        OUT / "scripts" / "k4_run.py",
    ]


def stopped_deliverable_paths() -> list[Path]:
    """Limit stopped reports to the completed non-self-referential stages."""
    paths = [requirements_path(), OUT / "scripts" / "k4_run.py"]
    for path in (OUT / "clusters.csv", OUT / "matches.csv"):
        if path.exists():
            paths.append(path)
    return paths


def report_gated_values(gates: dict[str, Any]) -> dict[str, Any]:
    """Select only K4 pre-registered numbers, retaining the allowed stopped subset."""
    keys = ["symmetry_max", *GATED_MATCH_KEYS]
    if not gates.get("stopped"):
        keys.append("k4_outcome")
    return {key: gates[key] for key in keys if key in gates}


def regeneration_commands(gated: dict[str, Any], deliverables: list[Path]) -> dict[str, str]:
    """Attach a single exact stage command to every result and gated number."""
    base = 'cd "$HOME/scopalign-cloud" && .venv/bin/python reports/K4/scripts/k4_run.py'
    commands: dict[str, str] = {
        "report.json": base + " --stage report",
        "requirements.lock": base + " --stage check",
        "clusters.csv": base + " --stage cluster",
        "matches.csv": base + " --stage match",
        "fig_anchors.png": base + " --stage report",
        "scripts/k4_run.py": base + " --stage all",
    }
    for path in deliverables:
        commands.setdefault(str(path.relative_to(OUT)), base + " --stage all")
    commands["symmetry_max"] = base + " --stage check"
    for key in GATED_MATCH_KEYS:
        if key in gated:
            commands[key] = base + " --stage match"
    if "k4_outcome" in gated:
        commands["k4_outcome"] = base + " --stage report"
    return commands


def report_stage() -> None:
    """Generate the required figure and schema-conformant completed or stopped report."""
    ensure_layout()
    meta = make_meta()
    write_requirements()
    gates = get_gates()
    if "symmetry_max" not in gates:
        raise RuntimeError("Run --stage check before --stage report")
    stopped = bool(gates.get("stopped", False))
    if not stopped:
        required = ("symmetry_max", *GATED_MATCH_KEYS)
        missing = [key for key in required if key not in gates]
        if missing:
            raise RuntimeError(f"Run --stage match before --stage report; missing {missing}")
        outcome = decide(gates)
        update_gates(k4_outcome=outcome)
        gates = get_gates()
        anchors_figure()
    gated = report_gated_values(gates)
    deliverable_paths = stopped_deliverable_paths() if stopped else complete_deliverable_paths()
    missing_deliverables = [str(path) for path in deliverable_paths if not path.exists()]
    if missing_deliverables:
        raise RuntimeError(f"Required K4 output is missing: {missing_deliverables}")
    report: dict[str, Any] = {
        "task_id": "K4",
        "executor": "manus",
        "branch": "task/K4",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": BASE_SEED,
        "started": meta["started"],
        "finished": utcnow(),
        "status": "stopped" if stopped else "complete",
        "gated": gated,
        "deliverables": [file_entry(path) for path in deliverable_paths],
        "artifacts": [file_entry(path) for path in artifact_paths()],
        "regeneration": regeneration_commands(gated, deliverable_paths),
        "deviations": [],
        "deviations_requested": [],
        "observations": [
            "K4 — Anchoring: can operator fingerprints match cell groups across batches without labels?",
            "Operator fingerprints use the saved float64 operators and saved checkpoint symmetry values only.",
            "Cluster labels are majority coarse cell_type values calculated after Leiden clustering for scoring and display.",
        ],
        "known_failures": KNOWN_FAILURES,
        "lockfile_sha256": sha256(requirements_path()),
    }
    if stopped:
        report["stop_rule"] = str(gates["stop_rule"])
    atomic_json(OUT / "report.json", report)
    print(json.dumps({"stage": "report", **gated, "status": report["status"]}, sort_keys=True), flush=True)


def staged_paths() -> list[str]:
    """Return currently staged repository paths for K4-only commit validation."""
    text = subprocess.check_output(
        ["git", "-C", str(ROOT), "diff", "--cached", "--name-only"], text=True
    )
    return [line for line in text.splitlines() if line]


def checkpoint_stage(stage: str) -> None:
    """Commit and push only reports/K4 after every completed --stage all checkpoint."""
    if os.environ.get("K4_AUTO_CHECKPOINT", "1") != "1":
        return
    branch = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "--abbrev-ref", "HEAD"], text=True
    ).strip()
    if branch != "task/K4":
        raise RuntimeError(f"K4 checkpoint requires task/K4, found {branch}")
    preexisting = staged_paths()
    non_k4 = [path for path in preexisting if not path.startswith("reports/K4/")]
    if non_k4:
        raise RuntimeError(f"Refusing to commit pre-staged paths outside reports/K4: {non_k4}")
    subprocess.run(["git", "-C", str(ROOT), "add", "--", "reports/K4"], check=True)
    staged = staged_paths()
    non_k4 = [path for path in staged if not path.startswith("reports/K4/")]
    if non_k4:
        raise RuntimeError(f"K4 checkpoint staged a path outside reports/K4: {non_k4}")
    if not staged:
        return
    subprocess.run(["git", "-C", str(ROOT), "commit", "-m", f"task: K4 {stage} stage"], check=True)
    subprocess.run(["git", "-C", str(ROOT), "push", "origin", "task/K4"], check=True)


def main() -> None:
    """Dispatch the four resumable K4 stages and their all-stage checkpoints."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("check", "cluster", "match", "report", "all"))
    stage = parser.parse_args().stage

    if stage in ("check", "all"):
        check_stage()
        if stage == "all":
            checkpoint_stage("check")
    if stage == "check":
        return

    if stage == "all" and is_stopped():
        report_stage()
        checkpoint_stage("report")
        return
    if stage in ("cluster", "all"):
        cluster_stage()
        if stage == "all":
            checkpoint_stage("cluster")
    if stage == "cluster":
        return

    if stage == "all" and is_stopped():
        report_stage()
        checkpoint_stage("report")
        return
    if stage in ("match", "all"):
        match_stage()
        if stage == "all":
            checkpoint_stage("match")
    if stage == "match":
        return

    if stage in ("report", "all"):
        report_stage()
        if stage == "all":
            checkpoint_stage("report")


if __name__ == "__main__":
    main()
