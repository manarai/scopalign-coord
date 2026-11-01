#!/usr/bin/env python3
"""K3: gene-network space from centroid operators, with cached T04 inputs."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scjdo2.compare import centroid_operator
from scjdo2.energy_drift import EnergyDriftAdapter, _EnergyMLP
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "K3"
CACHE = OUT / "cache"
K0 = ROOT / "reports" / "K0"
T04_ROOT = ROOT / "reports" / "T04"
T04_SCRIPT = T04_ROOT / "scripts" / "t04_run.py"
RAW = ROOT / "data" / "pbmc_10k_v3" / "pbmc_10k_v3.h5ad"
HALVES = ROOT / "data" / "pbmc_10k_v3" / "k0_halves.csv"
GENES_PATH = K0 / "genes.txt"
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
BASE_SEED = 20261014
ALT_SEEDS = SEEDS[1:]
D = 30
N_CENTROID = 30
NULL_SEED = 20261300
KNOWN = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(value)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json(path: Path, value: dict) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_npy(path: Path, value: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("wb") as handle:
            np.save(handle, value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_torch(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        torch.save(value, temporary)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def import_t04():
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("t04_k3_readonly", T04_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {T04_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T04 = import_t04()


def gates_path() -> Path:
    return OUT / "gates.json"


def get_gates() -> dict:
    return json.loads(gates_path().read_text()) if gates_path().exists() else {}


def update_gates(**values) -> dict:
    all_values = get_gates()
    all_values.update(values)
    atomic_json(gates_path(), all_values)
    return all_values


def stopped() -> bool:
    return bool(get_gates().get("stopped", False))


def requirements_path() -> Path:
    return OUT / "requirements.lock"


def write_requirements() -> None:
    if not requirements_path().exists():
        output = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
        atomic_text(requirements_path(), output)


def make_meta() -> dict:
    path = OUT / "run_meta.json"
    if path.exists():
        return json.loads(path.read_text())
    value = {
        "started": utcnow(), "seed": BASE_SEED, "drift_seeds": list(SEEDS),
        "null_rng": NULL_SEED,
        "null_traversal": ["network_kind operator then correlation", "comparison tech, genebatch, type, bio", "sorted type", "operator seed then split"],
        "operator_centroid_cells": N_CENTROID,
    }
    atomic_json(path, value)
    return value


def genes() -> list[str]:
    values = GENES_PATH.read_text().splitlines()
    if len(values) != 2000 or values != sorted(values):
        raise RuntimeError("K0 gene list must contain 2,000 sorted genes")
    return values


def raw_qc():
    value = sc.read_h5ad(RAW)
    return value[value.obs["qc_pass"].to_numpy(dtype=bool)].copy()


def p1_halves():
    source = raw_qc()
    assignments = pd.read_csv(HALVES).set_index("barcode")["half"]
    names = source.obs_names.astype(str)
    mapped = names.map(assignments)
    h1 = source[(names.isin(assignments.index)) & (mapped == "H1")].copy()
    h2 = source[(names.isin(assignments.index)) & (mapped == "H2")].copy()
    if h1.n_obs == 0 or h2.n_obs == 0:
        raise RuntimeError("K0 H1/H2 halves could not be reconstructed")
    return h1, h2


def normalised_dense(adata) -> np.ndarray:
    item = adata.copy()
    sc.pp.normalize_total(item, target_sum=1e4)
    sc.pp.log1p(item)
    matrix = item.X.toarray() if hasattr(item.X, "toarray") else np.asarray(item.X)
    return np.asarray(matrix)


def p1_cache_path() -> Path:
    return CACHE / "p1_latents.npz"


def rederive_p1():
    selected_genes = genes()
    h1, h2 = p1_halves()
    x1 = normalised_dense(h1[:, selected_genes].copy())
    x2 = normalised_dense(h2[:, selected_genes].copy())
    pca1 = PCA(n_components=D, svd_solver="full", random_state=BASE_SEED)
    pca2 = PCA(n_components=D, svd_solver="full", random_state=BASE_SEED)
    z1 = np.asarray(pca1.fit_transform(x1), dtype=np.float32)
    z2 = np.asarray(pca2.fit_transform(x2), dtype=np.float32)
    return {
        "z_h1": z1, "z_h2": z2,
        "w_h1": np.asarray(pca1.components_), "w_h2": np.asarray(pca2.components_),
        "types_h1": h1.obs["cell_type"].astype(str).to_numpy(dtype=str),
        "types_h2": h2.obs["cell_type"].astype(str).to_numpy(dtype=str),
        "x_h1": x1, "x_h2": x2,
    }


def save_p1_cache(values: dict) -> None:
    path = p1_cache_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("wb") as handle:
            np.savez_compressed(handle, z_h1=values["z_h1"], z_h2=values["z_h2"],
                                w_h1=values["w_h1"], w_h2=values["w_h2"],
                                types_h1=values["types_h1"], types_h2=values["types_h2"])
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def load_p1_cache() -> dict:
    with np.load(p1_cache_path(), allow_pickle=False) as arrays:
        return {name: arrays[name] for name in arrays.files}


def relative_latent_check(values: dict) -> dict:
    committed1 = np.load(K0 / "latent_pca30_h1.npy")
    committed2 = np.load(K0 / "latent_pca30_h2.npy")
    differences = [np.abs(values["z_h1"] - committed1), np.abs(values["z_h2"] - committed2)]
    maximum = float(max(float(x.max()) for x in differences))
    denominator = float(max(np.max(np.abs(committed1)), np.max(np.abs(committed2))))
    relative = maximum / denominator if denominator else float("inf")
    correlations = []
    for current, committed in ((values["z_h1"], committed1), (values["z_h2"], committed2)):
        for col in range(D):
            correlations.append(abs(float(np.corrcoef(current[:, col], committed[:, col])[0, 1])))
    result = {"latent_check_max": maximum, "latent_check_denominator": denominator,
              "latent_check_rel": relative, "latent_component_abs_r_min": float(min(correlations))}
    atomic_json(OUT / "latent_check.json", result)
    return result


def eligible_types(types_left: np.ndarray, types_right: np.ndarray) -> tuple[list[str], list[dict]]:
    shared = sorted(set(types_left.astype(str)) & set(types_right.astype(str)))
    rows, keep = [], []
    for label in shared:
        n_left = int((types_left.astype(str) == label).sum())
        n_right = int((types_right.astype(str) == label).sum())
        include = n_left >= 100 and n_right >= 100
        rows.append({"cell_type": label, "n_left": n_left, "n_right": n_right, "eligible": include})
        if include:
            keep.append(label)
    return keep, rows


def drift_target(half: str, seed: int) -> Path:
    return CACHE / "drifts" / f"P1_{half}_{seed}.pt"


def operator_target(half: str, seed: int) -> Path:
    return CACHE / "operators" / f"P1_{half}_{seed}.npy"


def validate_drift(path: Path) -> dict:
    saved = torch.load(path, map_location="cpu", weights_only=False)
    required = {"state_dict", "config", "symmetry"}
    if not required.issubset(saved):
        raise RuntimeError(f"Incomplete cached drift checkpoint: {path}")
    return saved


def p1_fit_worker(item: tuple[str, int, np.ndarray, str]) -> dict:
    half, seed, z, target_text = item
    torch.set_num_threads(1)
    target = Path(target_text)
    if target.exists():
        saved = validate_drift(target)
    else:
        config = T04.drift_config(seed)
        drift = EnergyDriftAdapter(D, config=config).fit(np.asarray(z, dtype=np.float32), verbose=True)
        symmetry = drift.symmetry_check(np.asarray(z, dtype=np.float32), n_cells=100)
        saved = {"state_dict": drift._model.state_dict(), "config": asdict(drift.config),
                 "symmetry": {key: float(value) for key, value in symmetry.items()}}
        atomic_torch(target, saved)
    return {"half": half, "seed": int(seed), "path": str(target.relative_to(ROOT)),
            "median_asym": float(saved["symmetry"]["median_asym"]),
            "max_asym": float(saved["symmetry"]["max_asym"])}


def model_from_checkpoint(path: Path):
    return T04.build_model(validate_drift(path))


def ensure_p1_operator(half: str, seed: int, z: np.ndarray) -> Path:
    target = operator_target(half, seed)
    expected = (len(z), D, D)
    if target.exists():
        loaded = np.load(target, mmap_mode="r")
        if loaded.shape == expected and loaded.dtype == np.float64:
            return target
        raise RuntimeError(f"Invalid cached P1 operator array: {target}")
    matrices = T04.individual_operator(model_from_checkpoint(drift_target(half, seed)), np.asarray(z, dtype=np.float64))
    matrices = (np.asarray(matrices, dtype=np.float64) + np.swapaxes(np.asarray(matrices, dtype=np.float64), 1, 2)) / 2.0
    if matrices.shape != expected:
        raise RuntimeError(f"Unexpected P1 operator shape {matrices.shape} for {half}/{seed}")
    atomic_npy(target, matrices)
    return target


def t04_input_paths() -> list[Path]:
    paths = [RAW, HALVES, GENES_PATH, K0 / "latent_pca30_h1.npy", K0 / "latent_pca30_h2.npy",
             T04_ROOT / "cache" / "data" / "S0.npz", T04_ROOT / "cache" / "embeddings" / "S0.npz",
             T04_ROOT / "cache" / "data" / "S1.npz", T04_ROOT / "cache" / "embeddings" / "S1.npz"]
    for scenario, sides in (("S0", ("A", "B")), ("S1", ("B",))):
        for side in sides:
            for seed in SEEDS:
                paths += [T04_ROOT / "cache" / "drifts" / f"{scenario}_{side}_{seed}.pt",
                          T04_ROOT / "cache" / "operators" / f"{scenario}_{side}_{seed}.npy"]
    return paths


def write_input_manifest() -> None:
    values = []
    for path in t04_input_paths():
        if not path.exists():
            raise RuntimeError(f"Required K3 input missing: {path}")
        values.append({"path": str(path.relative_to(ROOT)), "sha256": sha(path), "bytes": path.stat().st_size})
    atomic_json(OUT / "input_manifest.json", {"inputs": values})


def load_t04(scenario: str) -> dict:
    embedding_path = T04_ROOT / "cache" / "embeddings" / f"{scenario}.npz"
    data_path = T04_ROOT / "cache" / "data" / f"{scenario}.npz"
    with np.load(embedding_path, allow_pickle=False) as arrays:
        embedding = {name: arrays[name] for name in arrays.files}
    with np.load(data_path, allow_pickle=False) as arrays:
        data = {name: arrays[name] for name in arrays.files}
    n_a = len(embedding["z_a"])
    if len(embedding["z_b"]) + n_a != len(embedding["types"]):
        raise RuntimeError(f"{scenario} embedding labels do not match latent rows")
    for side, z_name, w_name, type_name in (("A", "z_a", "w_a", "types_a"), ("B", "z_b", "w_b", "types_b")):
        if embedding[z_name].shape[1] != D or embedding[w_name].shape != (D, 2000):
            raise RuntimeError(f"{scenario} {side} has unexpected latent/loading shape")
        side_types = embedding["types"][:n_a] if side == "A" else embedding["types"][n_a:]
        if len(data[type_name]) != len(side_types) or not np.array_equal(data[type_name].astype(str), side_types.astype(str)):
            raise RuntimeError(f"{scenario} {side} type labels do not align between caches")
        if data[side].shape[0] != len(side_types) or data[side].shape[1] != 2000:
            raise RuntimeError(f"{scenario} {side} expression counts do not align")
    return {"embedding": embedding, "data": data, "n_a": n_a}


def load_t04_operators(scenario: str, side: str) -> dict[int, np.ndarray]:
    values = {}
    for seed in SEEDS:
        path = T04_ROOT / "cache" / "operators" / f"{scenario}_{side}_{seed}.npy"
        array = np.load(path, mmap_mode="r")
        if array.ndim != 3 or array.shape[1:] != (D, D):
            raise RuntimeError(f"Unexpected T04 operator shape: {path}")
        values[seed] = array
    return values


def fits_stage() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    make_meta(); write_requirements(); write_input_manifest()
    p1 = rederive_p1()
    save_p1_cache(p1)
    latent = relative_latent_check(p1)
    p1_types, p1_rows = eligible_types(p1["types_h1"], p1["types_h2"])
    s0 = load_t04("S0")
    p2_types, p2_rows = eligible_types(s0["data"]["types_a"].astype(str), s0["data"]["types_b"].astype(str))
    if p1_types != p2_types:
        raise RuntimeError("P2 type eligibility does not equal the specified P1 type set")
    type_rows = []
    for row in p1_rows:
        type_rows.append({"pair": "P1", **row})
    for row in p2_rows:
        type_rows.append({"pair": "P2", **row})
    pd.DataFrame(type_rows).to_csv(OUT / "type_counts.csv", index=False)

    todo = [(half, seed, p1[f"z_{half}"] , str(drift_target(half, seed)))
            for half in ("h1", "h2") for seed in SEEDS]
    results = []
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=3, mp_context=context) as pool:
        futures = [pool.submit(p1_fit_worker, item) for item in todo]
        for future in as_completed(futures):
            results.append(future.result())
    results = sorted(results, key=lambda row: (row["half"], row["seed"]))
    pd.DataFrame(results).to_csv(OUT / "fit_symmetry.csv", index=False)
    for half in ("h1", "h2"):
        for seed in SEEDS:
            ensure_p1_operator(half, seed, p1[f"z_{half}"])
    symmetry = float(max(row["median_asym"] for row in results))
    gates = update_gates(n_types_P1=int(len(p1_types)), n_types_P2=int(len(p2_types)),
                         latent_check_rel=float(latent["latent_check_rel"]), symmetry_max=symmetry)
    if gates["n_types_P1"] < 3 or gates["latent_check_rel"] > 1e-4 or gates["symmetry_max"] > 1e-4:
        update_gates(stopped=True, stop_rule="n_types_P1 < 3, latent_check_rel > 1e-4 or symmetry_max > 1e-4 → REVISE.")
    print(json.dumps({key: get_gates()[key] for key in ("n_types_P1", "n_types_P2", "latent_check_rel", "symmetry_max")}, sort_keys=True), flush=True)


def op_descriptor(z: np.ndarray, matrices: dict[int, np.ndarray], w: np.ndarray, indices: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    selected = np.asarray(indices, dtype=int)
    centroid = centroid_operator(np.asarray(z[selected], dtype=np.float64), np.asarray(matrices[seed][selected], dtype=np.float64), n_cells_centroid=N_CENTROID).astype(np.float64)
    centroid = (centroid + centroid.T) / 2.0
    return centroid, np.asarray(w, dtype=np.float64)


def op_distance(left: tuple[np.ndarray, np.ndarray], right: tuple[np.ndarray, np.ndarray]) -> float:
    c1, w1 = left; c2, w2 = right
    g11 = w1 @ w1.T
    g22 = w2 @ w2.T
    g12 = w1 @ w2.T
    g21 = g12.T
    norm1 = float(np.trace(c1.T @ g11 @ c1 @ g11))
    norm2 = float(np.trace(c2.T @ g22 @ c2 @ g22))
    cross = float(np.trace(c1.T @ g12 @ c2 @ g21))
    difference = max(0.0, norm1 + norm2 - 2.0 * cross)
    return float(np.sqrt(difference / norm2)) if norm2 > 0 else float("inf")


def corr_descriptor(matrix: np.ndarray) -> tuple[np.ndarray, int]:
    values = np.asarray(matrix, dtype=np.float64)
    centered = values - values.mean(axis=0, keepdims=True)
    scale = np.sqrt(np.sum(centered ** 2, axis=0))
    good = scale > 0
    standardised = np.zeros_like(centered, dtype=np.float64)
    standardised[:, good] = centered[:, good] / scale[good]
    return standardised, int((~good).sum())


def gram_norm_sq(left: np.ndarray, right: np.ndarray) -> float:
    total = 0.0
    for start in range(0, len(left), 256):
        block = left[start:start + 256] @ right.T
        total += float(np.sum(block * block))
    return total


def corr_distance(left: np.ndarray, right: np.ndarray) -> float:
    norm1 = gram_norm_sq(left, left)
    norm2 = gram_norm_sq(right, right)
    cross = gram_norm_sq(left, right)
    difference = max(0.0, norm1 + norm2 - 2.0 * cross)
    return float(np.sqrt(difference / norm2)) if norm2 > 0 else float("inf")


def operator_null(z: np.ndarray, matrices: dict[int, np.ndarray], w: np.ndarray, mask: np.ndarray, rng: np.random.Generator) -> float:
    indices = np.flatnonzero(mask)
    values = []
    for seed in ALT_SEEDS:
        for _ in range(5):
            permutation = rng.permutation(indices)
            half = len(permutation) // 2
            values.append(op_distance(op_descriptor(z, matrices, w, permutation[:half], BASE_SEED),
                                      op_descriptor(z, matrices, w, permutation[half:2 * half], seed)))
    return float(np.median(values))


def corr_null(expression: np.ndarray, mask: np.ndarray, rng: np.random.Generator) -> float:
    indices = np.flatnonzero(mask)
    values = []
    for _ in range(20):
        permutation = rng.permutation(indices)
        half = len(permutation) // 2
        left, _ = corr_descriptor(expression[permutation[:half]])
        right, _ = corr_descriptor(expression[permutation[half:2 * half]])
        values.append(corr_distance(left, right))
    return float(np.median(values))


def observed_op(z1, m1, w1, mask1, z2, m2, w2, mask2) -> float:
    return op_distance(op_descriptor(z1, m1, w1, np.flatnonzero(mask1), BASE_SEED),
                       op_descriptor(z2, m2, w2, np.flatnonzero(mask2), BASE_SEED))


def observed_corr(x1, mask1, x2, mask2) -> tuple[float, int, int]:
    left, zero1 = corr_descriptor(x1[mask1])
    right, zero2 = corr_descriptor(x2[mask2])
    return corr_distance(left, right), zero1, zero2


def p1_expression() -> dict:
    selected_genes = genes()
    h1, h2 = p1_halves()
    return {"x_h1": normalised_dense(h1[:, selected_genes].copy()), "x_h2": normalised_dense(h2[:, selected_genes].copy())}


def compare_stage() -> None:
    if stopped():
        raise RuntimeError("K3 has a pre-registered validity stop; do not run compare")
    p1 = load_p1_cache()
    p1_x = p1_expression()
    s0 = load_t04("S0")
    s1 = load_t04("S1")
    p1_types = [row["cell_type"] for row in pd.read_csv(OUT / "type_counts.csv").query("pair == 'P1' and eligible").to_dict("records")]
    p1_m = {half: {seed: np.load(operator_target(half, seed), mmap_mode="r") for seed in SEEDS} for half in ("h1", "h2")}
    s0_m = {side: load_t04_operators("S0", side) for side in ("A", "B")}
    s1_m_b = load_t04_operators("S1", "B")
    n_a0, n_a1 = s0["n_a"], s1["n_a"]
    s0_types_a, s0_types_b = s0["data"]["types_a"].astype(str), s0["data"]["types_b"].astype(str)
    s1_types_b = s1["data"]["types_b"].astype(str)
    responders_b = s1["data"]["responder_b"].astype(bool)
    x_s0_a, x_s0_b = normalised_dense(ad.AnnData(np.asarray(s0["data"]["A"]))), normalised_dense(ad.AnnData(np.asarray(s0["data"]["B"])))
    x_s1_b = normalised_dense(ad.AnnData(np.asarray(s1["data"]["B"])))
    rng = np.random.default_rng(NULL_SEED)
    rows, zero_rows = [], []

    def add_pair(pair, comparison, label, op_observed, op_null, corr_observed, corr_null_value, zero_left, zero_right):
        rows.append({"pair": pair, "comparison": comparison, "network_kind": "operator", "type": label,
                     "observed": op_observed, "null": op_null, "ratio": op_observed / op_null})
        rows.append({"pair": pair, "comparison": comparison, "network_kind": "correlation", "type": label,
                     "observed": corr_observed, "null": corr_null_value, "ratio": corr_observed / corr_null_value})
        zero_rows.extend([
            {"pair": pair, "comparison": comparison, "type": label, "side": "left", "n_zero_variance": zero_left},
            {"pair": pair, "comparison": comparison, "type": label, "side": "right", "n_zero_variance": zero_right},
        ])

    # Fixed null-draw traversal: operator then correlation, comparison tech, genebatch, type, bio, sorted type.
    for kind in ("operator", "correlation"):
        for label in p1_types:
            h1_mask, h2_mask = p1["types_h1"].astype(str) == label, p1["types_h2"].astype(str) == label
            if kind == "operator":
                observed = observed_op(p1["z_h1"], p1_m["h1"], p1["w_h1"], h1_mask, p1["z_h2"], p1_m["h2"], p1["w_h2"], h2_mask)
                null = operator_null(p1["z_h2"], p1_m["h2"], p1["w_h2"], h2_mask, rng)
                rows.append({"pair": "P1", "comparison": "tech", "network_kind": kind, "type": label, "observed": observed, "null": null, "ratio": observed / null})
            else:
                observed, z1, z2 = observed_corr(p1_x["x_h1"], h1_mask, p1_x["x_h2"], h2_mask)
                null = corr_null(p1_x["x_h2"], h2_mask, rng)
                rows.append({"pair": "P1", "comparison": "tech", "network_kind": kind, "type": label, "observed": observed, "null": null, "ratio": observed / null})
                zero_rows.extend([{"pair": "P1", "comparison": "tech", "type": label, "side": "left", "n_zero_variance": z1}, {"pair": "P1", "comparison": "tech", "type": label, "side": "right", "n_zero_variance": z2}])
        for label in p1_types:
            a_mask, b_mask = s0_types_a == label, s0_types_b == label
            if kind == "operator":
                observed = observed_op(s0["embedding"]["z_a"], s0_m["A"], s0["embedding"]["w_a"], a_mask, s0["embedding"]["z_b"], s0_m["B"], s0["embedding"]["w_b"], b_mask)
                null = operator_null(s0["embedding"]["z_b"], s0_m["B"], s0["embedding"]["w_b"], b_mask, rng)
                rows.append({"pair": "P2", "comparison": "genebatch", "network_kind": kind, "type": label, "observed": observed, "null": null, "ratio": observed / null})
            else:
                observed, z1, z2 = observed_corr(x_s0_a, a_mask, x_s0_b, b_mask)
                null = corr_null(x_s0_b, b_mask, rng)
                rows.append({"pair": "P2", "comparison": "genebatch", "network_kind": kind, "type": label, "observed": observed, "null": null, "ratio": observed / null})
                zero_rows.extend([{"pair": "P2", "comparison": "genebatch", "type": label, "side": "left", "n_zero_variance": z1}, {"pair": "P2", "comparison": "genebatch", "type": label, "side": "right", "n_zero_variance": z2}])
        mono, tcell = "Monocytes", "T cells"
        h2_mono, h2_t = p1["types_h2"].astype(str) == mono, p1["types_h2"].astype(str) == tcell
        if not (h2_mono.any() and h2_t.any()):
            raise RuntimeError("P1 requires Monocytes and T cells for type comparison")
        if kind == "operator":
            observed = observed_op(p1["z_h2"], p1_m["h2"], p1["w_h2"], h2_mono, p1["z_h2"], p1_m["h2"], p1["w_h2"], h2_t)
            null = operator_null(p1["z_h2"], p1_m["h2"], p1["w_h2"], h2_t, rng)
            rows.append({"pair": "P1", "comparison": "type", "network_kind": kind, "type": "Monocytes_vs_T_cells", "observed": observed, "null": null, "ratio": observed / null})
        else:
            observed, z1, z2 = observed_corr(p1_x["x_h2"], h2_mono, p1_x["x_h2"], h2_t)
            null = corr_null(p1_x["x_h2"], h2_t, rng)
            rows.append({"pair": "P1", "comparison": "type", "network_kind": kind, "type": "Monocytes_vs_T_cells", "observed": observed, "null": null, "ratio": observed / null})
            zero_rows.extend([{"pair": "P1", "comparison": "type", "type": "Monocytes_vs_T_cells", "side": "left", "n_zero_variance": z1}, {"pair": "P1", "comparison": "type", "type": "Monocytes_vs_T_cells", "side": "right", "n_zero_variance": z2}])
        resp = (s1_types_b == mono) & responders_b
        nonresp = (s1_types_b == mono) & ~responders_b
        if not (resp.any() and nonresp.any()):
            raise RuntimeError("S1 B requires responder and non-responder monocytes")
        if kind == "operator":
            observed = observed_op(s1["embedding"]["z_b"], s1_m_b, s1["embedding"]["w_b"], resp, s1["embedding"]["z_b"], s1_m_b, s1["embedding"]["w_b"], nonresp)
            null = operator_null(s1["embedding"]["z_b"], s1_m_b, s1["embedding"]["w_b"], nonresp, rng)
            rows.append({"pair": "S1", "comparison": "bio", "network_kind": kind, "type": "Monocytes_responders_vs_nonresponders", "observed": observed, "null": null, "ratio": observed / null})
        else:
            observed, z1, z2 = observed_corr(x_s1_b, resp, x_s1_b, nonresp)
            null = corr_null(x_s1_b, nonresp, rng)
            rows.append({"pair": "S1", "comparison": "bio", "network_kind": kind, "type": "Monocytes_responders_vs_nonresponders", "observed": observed, "null": null, "ratio": observed / null})
            zero_rows.extend([{"pair": "S1", "comparison": "bio", "type": "Monocytes_responders_vs_nonresponders", "side": "left", "n_zero_variance": z1}, {"pair": "S1", "comparison": "bio", "type": "Monocytes_responders_vs_nonresponders", "side": "right", "n_zero_variance": z2}])

    table = pd.DataFrame(rows, columns=["pair", "comparison", "network_kind", "type", "observed", "null", "ratio"])
    table.to_csv(OUT / "ratios.csv", index=False)
    pd.DataFrame(zero_rows).to_csv(OUT / "zero_variance.csv", index=False)
    values = {}
    for kind, prefix in (("operator", "net"), ("correlation", "corr")):
        for comparison in ("tech", "genebatch", "type", "bio"):
            subset = table[(table.network_kind == kind) & (table.comparison == comparison)]
            values[f"{prefix}_ratio_{comparison}"] = float(np.median(subset.ratio.to_numpy(dtype=float)))
    update_gates(**values)
    print(json.dumps(values, sort_keys=True), flush=True)


def figure_stage(values: dict) -> None:
    labels = ["tech", "genebatch", "type", "bio"]
    operator = [values[f"net_ratio_{label}"] for label in labels]
    correlation = [values[f"corr_ratio_{label}"] for label in labels]
    positions = np.arange(len(labels)); width = 0.36
    fig, axis = plt.subplots(figsize=(9, 5))
    axis.bar(positions - width / 2, operator, width, label="operator network", color="C3")
    axis.bar(positions + width / 2, correlation, width, label="correlation network", color="C0")
    axis.axhline(1.5, color="0.25", ls="--", lw=1, label="1.5")
    axis.axhline(3.0, color="0.25", ls=":", lw=1.2, label="3")
    axis.set_xticks(positions, labels)
    axis.set_ylabel("observed / refit-or-sampling null")
    axis.set_title("K3 gene-network contrast ratios")
    axis.legend(frameon=False, ncol=4)
    axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "fig_network_ratios.png", dpi=180)
    plt.close(fig)


def decide(values: dict) -> str:
    if values["net_ratio_tech"] > 2.0 or values["net_ratio_genebatch"] > 2.0:
        return "batch_sensitive"
    if values["net_ratio_type"] < 3.0 or values["net_ratio_bio"] < 1.5:
        return "uninformative"
    if (values["corr_ratio_tech"] <= values["net_ratio_tech"] and
        values["corr_ratio_genebatch"] <= values["net_ratio_genebatch"] and
        values["corr_ratio_type"] >= values["net_ratio_type"] and
        values["corr_ratio_bio"] >= values["net_ratio_bio"]):
        return "no_better_than_correlation"
    return "promising"


def report_stage() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = make_meta(); write_requirements()
    gates = get_gates()
    base = 'cd "$HOME/scopalign-cloud" && .venv/bin/python reports/K3/scripts/k3_run.py'
    if gates.get("stopped"):
        gated = {key: gates[key] for key in ("n_types_P1", "n_types_P2", "latent_check_rel", "symmetry_max") if key in gates}
        files = ["requirements.lock", "scripts/k3_run.py"]
        report = {"task_id": "K3", "executor": "manus", "branch": "task/K3",
                  "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
                  "seed": BASE_SEED, "started": meta["started"], "finished": utcnow(), "status": "stopped",
                  "stop_rule": gates["stop_rule"], "gated": gated,
                  "deliverables": [{"path": name, "sha256": sha(OUT / name)} for name in files],
                  "artifacts": [{"path": name, "sha256": sha(OUT / name)} for name in ("input_manifest.json", "latent_check.json", "type_counts.csv", "fit_symmetry.csv", "run_meta.json") if (OUT / name).exists()],
                  "regeneration": {**{name: base + " --stage fits" for name in files}, **{name: base + " --stage fits" for name in gated}},
                  "deviations": [], "deviations_requested": [],
                  "observations": ["K3 stopped at its pre-registered validity rule before network comparison."],
                  "known_failures": KNOWN, "lockfile_sha256": sha(requirements_path())}
        atomic_json(OUT / "report.json", report)
        print(json.dumps({"stage": "report", "status": "stopped", **gated}, sort_keys=True), flush=True)
        return
    required = ("n_types_P1", "n_types_P2", "latent_check_rel", "symmetry_max", "net_ratio_tech", "net_ratio_genebatch", "net_ratio_type", "net_ratio_bio", "corr_ratio_tech", "corr_ratio_genebatch", "corr_ratio_type", "corr_ratio_bio")
    if any(key not in gates for key in required):
        raise RuntimeError("Run fits and compare before report")
    gated = {key: gates[key] for key in required}
    gated["k3_outcome"] = decide(gated)
    update_gates(**gated)
    figure_stage(gated)
    files = ["requirements.lock", "ratios.csv", "fig_network_ratios.png", "scripts/k3_run.py"]
    artifacts = ["input_manifest.json", "latent_check.json", "type_counts.csv", "fit_symmetry.csv", "zero_variance.csv", "run_meta.json"]
    fit_gates = {"n_types_P1", "n_types_P2", "latent_check_rel", "symmetry_max"}
    report = {"task_id": "K3", "executor": "manus", "branch": "task/K3",
              "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
              "seed": BASE_SEED, "started": meta["started"], "finished": utcnow(), "status": "complete", "gated": gated,
              "deliverables": [{"path": name, "sha256": sha(OUT / name)} for name in files],
              "artifacts": [{"path": name, "sha256": sha(OUT / name)} for name in artifacts if (OUT / name).exists()],
              "regeneration": {**{name: base + " --stage all" for name in files}, **{name: base + (" --stage report" if name == "k3_outcome" else (" --stage fits" if name in fit_gates else " --stage compare")) for name in gated}},
              "deviations": [], "deviations_requested": [],
              "observations": ["Operator and correlation-network distances were evaluated without saving 2000 by 2000 network arrays."],
              "known_failures": KNOWN, "lockfile_sha256": sha(requirements_path())}
    atomic_json(OUT / "report.json", report)
    print(json.dumps({"stage": "report", **gated}, sort_keys=True), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("fits", "compare", "report", "all"))
    stage = parser.parse_args().stage
    if stage in ("fits", "all"):
        fits_stage()
    if stage in ("compare", "all") and not stopped():
        compare_stage()
    if stage in ("report", "all"):
        report_stage()


if __name__ == "__main__":
    main()
