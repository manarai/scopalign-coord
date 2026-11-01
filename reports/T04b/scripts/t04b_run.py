#!/usr/bin/env python3
"""T04b: semi-synthetic batch × biology sweep using T04's fixed design."""
import argparse
import multiprocessing as mp
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import time
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import scvi
import harmonypy as hm
import threadpoolctl

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "T04b"
CACHE = OUT / "cache"
T04_ROOT = ROOT / "reports" / "T04"
T04_SCRIPT = T04_ROOT / "scripts" / "t04_run.py"
T04_SCORES = T04_ROOT / "scores.csv"
D = 30
FS = (0, 1, 2, 4, 8)
SIGMAS = (0.25, 0.5, 1.0)
BATCH_TYPES = ("global", "celltype")
METHODS = ("none", "joint", "harmony", "scvi", "gene_rot")
SWEEP_WORKERS = 2
TORCH_THREADS_PER_SWEEP_WORKER = 3


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def import_t04():
    """Import T04 by path without writing into the read-only T04 report."""
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("t04_reuse_only", T04_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import T04 runner at {T04_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T04 = import_t04()


def f_label(f: int) -> str:
    return f"f{int(f)}"


def sigma_label(sigma: float) -> str:
    return {0.25: "s0p25", 0.5: "s0p5", 1.0: "s1p0"}[float(sigma)]


def key(batch_type: str, sigma: float, f: int) -> str:
    return f"{batch_type}__{sigma_label(sigma)}__{f_label(f)}"


def condition_rows():
    """Fixed sweep processing order: the two T04 reference conditions first."""
    reference = [("global", 0.5, 0), ("global", 0.5, 4)]
    remaining = [(batch, sigma, f) for batch in BATCH_TYPES for sigma in SIGMAS for f in FS
                 if (batch, sigma, f) not in reference]
    return reference + remaining


def data_draw_rows():
    """Fixed batch-Poisson draw order after the two legacy T04 draws and A draw."""
    reference = {("global", 0.5, 0), ("global", 0.5, 4)}
    return [(batch, sigma, f) for batch in BATCH_TYPES for sigma in SIGMAS for f in FS
            if (batch, sigma, f) not in reference]


def data_path(batch_type: str, sigma: float, f: int) -> Path:
    return CACHE / "data" / f"{key(batch_type, sigma, f)}.npz"


def oracle_path(f: int) -> Path:
    # Positive-f batch-free oracle is independent of batch type. The two f=0
    # aliases below are stored only to make the requested six oracle archives.
    return CACHE / "oracles" / f"oracle_{f_label(f)}.npz"


def zero_oracle_path(batch_type: str) -> Path:
    return CACHE / "oracles" / f"oracle_f0_{batch_type}.npz"


def condition_result_path(batch_type: str, sigma: float, f: int) -> Path:
    return CACHE / "sweep" / f"{key(batch_type, sigma, f)}.csv"


def continuation_path() -> Path:
    return OUT / "sweep_continuation.json"


def gates_path() -> Path:
    return OUT / "gates.json"


def get_gates() -> dict:
    return json.loads(gates_path().read_text()) if gates_path().exists() else {}


def update_gates(**values) -> dict:
    values_all = get_gates()
    values_all.update(values)
    gates_path().write_text(json.dumps(values_all, indent=2, sort_keys=True) + "\n")
    return values_all


def make_meta() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    (OUT / ".gitignore").write_text("cache/\n__pycache__/\ngates.json\nrun.log\nrun.pid\n")
    path = OUT / "run_meta.json"
    if not path.exists():
        path.write_text(json.dumps({
            "started": utcnow(),
            "platform": platform.platform(),
            "blas": threadpoolctl.threadpool_info(),
            "openblas_num_threads_env": os.environ.get("OPENBLAS_NUM_THREADS"),
            "omp_num_threads_env": os.environ.get("OMP_NUM_THREADS"),
            "t04_script_sha256": sha(T04_SCRIPT),
            "t04_scores_sha256": sha(T04_SCORES),
            "condition_processing_order": [{"batch_type": b, "sigma": s, "f": f} for b, s, f in condition_rows()],
        }, indent=2, sort_keys=True) + "\n")
    return json.loads(path.read_text())


def write_requirements() -> None:
    with (OUT / "requirements.lock").open("w") as handle:
        subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, stdout=handle)


def save_npz_if_missing(path: Path, **arrays) -> None:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **arrays)


def load_npz(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as arrays:
        return {name: arrays[name] for name in arrays.files}


def source_metadata(genes: list[str]) -> dict:
    return {
        "source_sha256": T04.sha(T04.SOURCE),
        "source_fingerprint": T04.fingerprint(T04.SOURCE),
        "halves_sha256": T04.sha(T04.HALVES),
        "halves_fingerprint": T04.fingerprint(T04.HALVES),
        "genes_sha256": T04.sha(T04.GENES_PATH),
        "t04_script_sha256": sha(T04_SCRIPT),
        "t04_scores_sha256": sha(T04_SCORES),
        "n_genes": int(len(genes)),
    }


def biological_matrix(b_raw: np.ndarray, a_raw: np.ndarray, mono_a: np.ndarray,
                      responder_idx: np.ndarray, bio_positions: np.ndarray, f: int) -> np.ndarray:
    out = b_raw.copy()
    if f == 0:
        return out
    mu = a_raw[mono_a][:, bio_positions].mean(axis=0, dtype=np.float64)
    depth_b = b_raw.sum(axis=1, dtype=np.float64)
    nbar = float(depth_b.mean())
    rates = float(f) * (mu[None, :] + 0.1) * (depth_b[responder_idx, None] / nbar)
    # Fresh seed per f: f=4 has precisely T04's only biology Poisson draw.
    increment = np.random.default_rng(20261102).poisson(rates).astype(np.int64)
    out[responder_idx[:, None], bio_positions[None, :]] += increment
    return out


def data_stage() -> None:
    make_meta()
    write_requirements()
    genes = T04.read_genes()
    a, b = T04.source_halves()
    a_raw, b_raw = T04.as_dense_counts(a, genes), T04.as_dense_counts(b, genes)
    types_a = a.obs["cell_type"].astype(str).to_numpy(dtype=str)
    types_b = b.obs["cell_type"].astype(str).to_numpy(dtype=str)
    mono_a, mono_b = types_a == "Monocytes", types_b == "Monocytes"
    b_mono_idx = np.flatnonzero(mono_b)
    responder_idx = np.sort(np.random.default_rng(20261100).choice(b_mono_idx, len(b_mono_idx) // 2, replace=False))
    responder = np.zeros(len(b), dtype=bool)
    responder[responder_idx] = True
    bio_positions = np.asarray([genes.index(gene) for gene in T04.BIO_GENES if gene in genes], dtype=np.int64)
    if len(bio_positions) != len(T04.BIO_GENES):
        raise RuntimeError("All 22 T04 interferon genes must be in genes.txt")
    b_biology = {f: biological_matrix(b_raw, a_raw, mono_a, responder_idx, bio_positions, f) for f in FS}

    # Exact legacy T04 draw stream: standard-normal eps then f=0 B, f=4 B,
    # and one A resample. The two named condition archives must match T04.
    batch_rng = np.random.default_rng(20261101)
    epsilon = batch_rng.normal(0.0, 1.0, len(genes))
    factor_ref = np.exp(0.5 * epsilon)
    b_ref_f0 = batch_rng.poisson(b_raw.astype(np.float64) * factor_ref[None, :] * 0.7).astype(np.int64)
    b_ref_f4 = batch_rng.poisson(b_biology[4].astype(np.float64) * factor_ref[None, :] * 0.7).astype(np.int64)
    a_once = batch_rng.poisson(a_raw).astype(np.int64)

    t04_s0, t04_s1, t04_oracle = T04.load_case("S0"), T04.load_case("S1"), T04.load_oracle()
    reference_arrays_match = bool(np.array_equal(a_once, t04_s0["A"]) and np.array_equal(a_once, t04_s1["A"])
                                  and np.array_equal(b_ref_f0, t04_s0["B"]) and np.array_equal(b_ref_f4, t04_s1["B"])
                                  and np.array_equal(b_biology[4], t04_oracle["B"]))
    if not reference_arrays_match:
        raise RuntimeError("T04 reference count construction did not reproduce cached S0/S1/oracle arrays")

    type_names = sorted(set(types_b.astype(str)))
    type_index = {name: index for index, name in enumerate(type_names)}
    b_type_rows = np.asarray([type_index[name] for name in types_b.astype(str)], dtype=np.int64)
    eta = np.random.default_rng(20261103).normal(0.0, 1.0, (len(type_names), len(genes)))
    type_factor = np.exp(0.3 * eta[b_type_rows])
    shared = {
        "types_a": types_a,
        "types_b": types_b,
        "responder_b": responder,
        "barcodes_a": a.obs_names.astype(str).to_numpy(dtype=str),
        "barcodes_b": b.obs_names.astype(str).to_numpy(dtype=str),
    }

    # Write the legacy reference-condition payloads first. Other Poisson draws
    # continue the same batch RNG stream in the recorded canonical order.
    save_npz_if_missing(data_path("global", 0.5, 0), A=a_once, B=b_ref_f0, **shared)
    save_npz_if_missing(data_path("global", 0.5, 4), A=a_once, B=b_ref_f4, **shared)
    for batch_type, sigma, f in data_draw_rows():
        factor = np.exp(float(sigma) * epsilon)
        if batch_type == "celltype":
            factor = factor[None, :] * type_factor
        else:
            factor = factor[None, :]
        b_draw = batch_rng.poisson(b_biology[f].astype(np.float64) * factor * 0.7).astype(np.int64)
        save_npz_if_missing(data_path(batch_type, sigma, f), A=a_once, B=b_draw, **shared)

    # Six no-batch oracle archives requested by T04b: four unique positive-f
    # denominators plus f=0 aliases for the two batch types.
    for f in (1, 2, 4, 8):
        save_npz_if_missing(oracle_path(f), A=a_once, B=b_biology[f], **shared)
    for batch_type in BATCH_TYPES:
        save_npz_if_missing(zero_oracle_path(batch_type), A=a_once, B=b_biology[0], **shared)

    condition_records = []
    for batch_type, sigma, f in condition_rows():
        path = data_path(batch_type, sigma, f)
        condition_records.append({"batch_type": batch_type, "sigma": float(sigma), "f": int(f),
                                  "path": str(path.relative_to(OUT)), "sha256": sha(path)})
    oracle_records = []
    for path in [oracle_path(f) for f in (1, 2, 4, 8)] + [zero_oracle_path(b) for b in BATCH_TYPES]:
        oracle_records.append({"path": str(path.relative_to(OUT)), "sha256": sha(path)})
    manifest = {"conditions": condition_records, "oracles": oracle_records}
    (OUT / "data_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    meta = source_metadata(genes)
    meta.update({
        "n_conditions": len(condition_records), "n_oracles": len(oracle_records),
        "n_bio_genes": int(len(bio_positions)), "n_responders": int(responder.sum()),
        "n_cells_a": int(len(a)), "n_cells_b": int(len(b)),
        "n_monocytes_a": int(mono_a.sum()), "n_monocytes_b": int(mono_b.sum()),
        "responder_seed": 20261100, "biology_seed": 20261102, "batch_seed": 20261101,
        "celltype_factor_seed": 20261103, "depth_factor": 0.7,
        "bio_genes": list(T04.BIO_GENES), "batch_types": list(BATCH_TYPES),
        "sigma_values": list(SIGMAS), "biology_f_values": list(FS),
        "batch_draw_order": ["epsilon_standard_normal", "T04_global_sigma0.5_f0_B",
                             "T04_global_sigma0.5_f4_B", "A_once", "remaining_conditions_in_data_draw_order"],
        "data_draw_order": [{"batch_type": btype, "sigma": sigma, "f": f} for btype, sigma, f in data_draw_rows()],
        "celltype_eta_type_order": type_names,
        "reference_arrays_match_t04": reference_arrays_match,
    })
    (OUT / "data_meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"stage": "data", "n_conditions": len(condition_records), "n_oracles": len(oracle_records),
                      "n_responders": int(responder.sum()), "reference_arrays_match_t04": reference_arrays_match}, sort_keys=True), flush=True)


def oracle_bio_keep(f: int) -> float:
    oracle = load_npz(oracle_path(f))
    dataset = np.asarray(["A"] * len(oracle["A"]) + ["B"] * len(oracle["B"]), dtype=str)
    types = np.concatenate([oracle["types_a"].astype(str), oracle["types_b"].astype(str)])
    responder = np.concatenate([np.zeros(len(oracle["A"]), dtype=bool), oracle["responder_b"].astype(bool)])
    matrix = np.vstack([T04.dense_log_counts(oracle["A"]), T04.dense_log_counts(oracle["B"])])
    rep = PCA(n_components=D, svd_solver="full", random_state=20261014).fit_transform(matrix).astype(np.float64)
    # dataset intentionally retained for the exact T04 data layout but not used by bio_keep.
    _ = dataset
    return float(T04.bio_keep(rep, types, responder))


def metrics_for_condition(batch_type: str, sigma: float, f: int) -> pd.DataFrame:
    condition = load_npz(data_path(batch_type, sigma, f))
    a_counts, b_counts = condition["A"], condition["B"]
    xa, xb = T04.dense_log_counts(a_counts), T04.dense_log_counts(b_counts)
    pca_a = PCA(n_components=D, svd_solver="full", random_state=20261014)
    pca_b = PCA(n_components=D, svd_solver="full", random_state=20261014)
    za = pca_a.fit_transform(xa).astype(np.float64)
    zb = pca_b.fit_transform(xb).astype(np.float64)
    joint = PCA(n_components=D, svd_solver="full", random_state=20261014).fit_transform(np.vstack([xa, xb])).astype(np.float64)
    dataset = np.asarray(["A"] * len(za) + ["B"] * len(zb), dtype=str)
    meta = pd.DataFrame({"dataset": dataset})
    harmony = hm.run_harmony(joint, meta, "dataset", random_state=20261014, max_iter_harmony=20).Z_corr.T.astype(np.float64)
    scvi.settings.seed = 20261014
    scvi_data = T04.raw_scvi_adata(a_counts, b_counts, T04.read_genes())
    scvi.model.SCVI.setup_anndata(scvi_data, batch_key="dataset")
    scvi_model = scvi.model.SCVI(scvi_data, n_latent=D)
    scvi_model.train(max_epochs=400, accelerator="cpu", devices=1)
    scvi_latent = scvi_model.get_latent_representation().astype(np.float64)
    u, _, vt = np.linalg.svd(pca_b.components_.astype(np.float64) @ pca_a.components_.astype(np.float64).T, full_matrices=False)
    q_load = (u @ vt).astype(np.float64)
    reps = {
        "none": np.vstack([za, zb]).astype(np.float64),
        "joint": joint,
        "harmony": harmony,
        "scvi": scvi_latent,
        "gene_rot": np.vstack([za @ q_load.T, zb]).astype(np.float64),
    }
    types = np.concatenate([condition["types_a"].astype(str), condition["types_b"].astype(str)])
    responder = np.concatenate([np.zeros(len(za), dtype=bool), condition["responder_b"].astype(bool)])
    oracle = oracle_bio_keep(f) if f > 0 else float("nan")
    rows = []
    for method in METHODS:
        score = T04.score_rep(reps[method], dataset, types, responder)
        rows.append({
            "batch_type": batch_type, "sigma": float(sigma), "f": int(f), "method": method,
            "mix_bg": float(score["mix_bg"]),
            "retention": float(score["bio_score"] / oracle) if f > 0 else float("nan"),
            "bio_false": float(score["bio_score"]) if f == 0 else float("nan"),
            "sil": float(score["sil"]),
            "oracle_bio_keep": float(oracle) if f > 0 else float("nan"),
        })
    return pd.DataFrame(rows)


def checkpoint_is_complete(path: Path, batch_type: str, sigma: float, f: int) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        rows = pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return False
    required = {"batch_type", "sigma", "f", "method", "mix_bg", "retention", "bio_false", "sil"}
    if not required.issubset(rows.columns) or len(rows) != len(METHODS):
        return False
    return (set(rows["method"].astype(str)) == set(METHODS)
            and rows["batch_type"].astype(str).eq(batch_type).all()
            and np.allclose(rows["sigma"].to_numpy(dtype=float), float(sigma), rtol=0.0, atol=0.0)
            and rows["f"].to_numpy(dtype=int).tolist() == [int(f)] * len(METHODS))


def checkpoint_records() -> list[dict]:
    records = []
    for batch_type, sigma, f in condition_rows():
        target = condition_result_path(batch_type, sigma, f)
        if checkpoint_is_complete(target, batch_type, sigma, f):
            records.append({"batch_type": batch_type, "sigma": float(sigma), "f": int(f),
                            "path": str(target.relative_to(OUT)), "sha256": sha(target)})
    return records


def atomic_write_csv(table: pd.DataFrame, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            table.to_csv(handle, index=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
        directory_fd = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary.exists():
            temporary.unlink()


_THREAD_LIMITER = None


def init_sweep_worker() -> None:
    """Set T04's three-Torch-thread setting in each clean spawned worker."""
    import torch

    torch.set_num_threads(TORCH_THREADS_PER_SWEEP_WORKER)
    torch.set_num_interop_threads(TORCH_THREADS_PER_SWEEP_WORKER)


def configure_spawn_environment() -> None:
    # T04's serial run had no explicit native-BLAS environment caps: NumPy's
    # OpenBLAS pool therefore used the host default (six threads). Preserve
    # that exact numerical regime; torch's own thread setting is fixed above.
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.pop(name, None)


def preflight_condition_worker(spec: tuple[str, float, int]) -> pd.DataFrame:
    """Rerun a completed reference condition without publishing a checkpoint."""
    return metrics_for_condition(*spec)


def sweep_condition_worker(spec: tuple[str, float, int]) -> dict:
    """Evaluate one independent condition and atomically save its checkpoint."""
    batch_type, sigma, f = spec
    target = condition_result_path(batch_type, sigma, f)
    if checkpoint_is_complete(target, batch_type, sigma, f):
        return {"batch_type": batch_type, "sigma": sigma, "f": f, "skipped": True, "elapsed_seconds": 0.0}
    started = time.perf_counter()
    table = metrics_for_condition(batch_type, sigma, f)
    atomic_write_csv(table, target)
    return {"batch_type": batch_type, "sigma": sigma, "f": f, "skipped": False,
            "elapsed_seconds": round(time.perf_counter() - started, 3), "pid": os.getpid()}


def complete_checkpoints() -> list[tuple[str, float, int]]:
    finished = []
    for batch_type, sigma, f in condition_rows():
        path = condition_result_path(batch_type, sigma, f)
        if checkpoint_is_complete(path, batch_type, sigma, f):
            finished.append((batch_type, sigma, f))
    return finished


def collect_sweep() -> pd.DataFrame:
    paths = [condition_result_path(batch_type, sigma, f) for batch_type, sigma, f in condition_rows()]
    if not all(path.exists() for path in paths):
        raise RuntimeError("Not all 30 per-condition sweep checkpoints exist")
    table = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    table = table.sort_values(["batch_type", "sigma", "f", "method"]).reset_index(drop=True)
    if len(table) != 150 or table.groupby(["batch_type", "sigma", "f"])["method"].nunique().ne(5).any():
        raise RuntimeError("Sweep cache does not contain five method rows for every condition")
    table.to_csv(OUT / "sweep.csv", index=False)
    return table


def table_value(table: pd.DataFrame, method: str, metric: str) -> float:
    selected = table.loc[table["method"].astype(str) == method, metric]
    if len(selected) != 1:
        raise RuntimeError(f"Expected one row for {method}/{metric}, found {len(selected)}")
    return float(selected.iloc[0])


def parallel_reference_preflight() -> dict:
    """Require spawned workers to reproduce the two existing T04 conditions."""
    references = (("global", 0.5, 0), ("global", 0.5, 4))
    configure_spawn_environment()
    context = mp.get_context("spawn")
    results: dict[tuple[str, float, int], pd.DataFrame] = {}
    with ProcessPoolExecutor(max_workers=SWEEP_WORKERS, mp_context=context, initializer=init_sweep_worker) as pool:
        futures = {pool.submit(preflight_condition_worker, spec): spec for spec in references}
        for future in as_completed(futures):
            results[futures[future]] = future.result()
    rows = []
    strict_methods = ("none", "joint", "harmony", "gene_rot")
    for batch_type, sigma, f in references:
        saved = pd.read_csv(condition_result_path(batch_type, sigma, f))
        rerun = results[(batch_type, sigma, f)]
        metrics = ("mix_bg", "bio_false", "sil") if f == 0 else ("mix_bg", "retention", "sil")
        for method in METHODS:
            for metric in metrics:
                prior = table_value(saved, method, metric)
                current = table_value(rerun, method, metric)
                rows.append({"batch_type": batch_type, "sigma": sigma, "f": f, "method": method,
                             "metric": metric, "saved": prior, "rerun": current,
                             "abs_difference": abs(prior - current)})
    strict_diffs = [row["abs_difference"] for row in rows if row["method"] in strict_methods]
    scvi_diffs = [row["abs_difference"] for row in rows if row["method"] == "scvi"]
    maximum = float(max(strict_diffs))
    if maximum > 1e-6:
        raise RuntimeError(f"Spawn-worker strict reference preflight failed: maximum absolute difference {maximum}")
    return {"conditions": [{"batch_type": b, "sigma": s, "f": f} for b, s, f in references],
            "comparisons": rows, "strict_methods_max_abs": maximum,
            "scvi_max_abs": float(max(scvi_diffs))}


def json_atomic_write(path: Path, value: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def previous_serial_runner_sha256() -> str | None:
    try:
        serial = subprocess.check_output(
            ["git", "-C", str(ROOT), "show", "HEAD^:reports/T04b/scripts/t04b_run.py"])
    except subprocess.CalledProcessError:
        return None
    return hashlib.sha256(serial).hexdigest()


def create_parallel_continuation(pending: list[tuple[str, float, int]]) -> dict:
    value = {
        "created": utcnow(),
        "pre_parallel_runner_sha256": previous_serial_runner_sha256(),
        "parallel_runner_sha256": sha(Path(__file__)),
        "requirements_lock_sha256": sha(OUT / "requirements.lock"),
        "platform": platform.platform(),
        "workers": SWEEP_WORKERS,
        "start_method": "spawn",
        "threads_per_worker": TORCH_THREADS_PER_SWEEP_WORKER,
        "thread_environment": {name: os.environ.get(name) for name in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
        "pre_stop_valid_checkpoints": checkpoint_records(),
        "invalid_or_recomputed_targets": [],
        "canonical_submission_order": [{"batch_type": b, "sigma": s, "f": f} for b, s, f in condition_rows()],
        "submitted_conditions": [{"batch_type": b, "sigma": s, "f": f} for b, s, f in pending],
        "completion_order": [],
        "unchanged_method_configuration": {
            "data_archives": "cached condition matrices created by --stage data",
            "scvi_seed": 20261014, "scvi_n_latent": D, "scvi_max_epochs": 400,
            "scvi_accelerator": "cpu", "scvi_devices": 1,
            "harmony_seed": 20261014, "harmony_max_iter": 20,
            "scores": ["mix_bg", "retention", "bio_false", "sil"],
        },
    }
    value["preflight"] = parallel_reference_preflight()
    json_atomic_write(continuation_path(), value)
    return value


def reproduction_check(table: pd.DataFrame) -> tuple[float, float, pd.DataFrame]:
    t04 = pd.read_csv(T04_SCORES)
    mapping = [("global", 0.5, 0, "S0", "bio_false"), ("global", 0.5, 4, "S1", "retention")]
    rows = []
    for batch_type, sigma, f, scenario, biology_column in mapping:
        current = table[(table.batch_type == batch_type) & (table.sigma == sigma) & (table.f == f)]
        expected = t04[t04.scenario == scenario]
        for method in METHODS:
            left, right = current[current.method == method], expected[expected.method == method]
            if len(left) != 1 or len(right) != 1:
                raise RuntimeError(f"Missing T04 reproduction row for {batch_type}/{sigma}/{f}/{method}")
            for metric in ("mix_bg", biology_column, "sil"):
                observed, target = float(left.iloc[0][metric]), float(right.iloc[0][metric])
                if not np.isfinite(observed) or not np.isfinite(target):
                    raise RuntimeError(f"Non-finite reproduction metric {metric} for {method}")
                rows.append({"batch_type": batch_type, "sigma": sigma, "f": f, "scenario": scenario,
                             "method": method, "metric": metric, "t04b": observed, "t04": target,
                             "abs_difference": abs(observed - target)})
    detail = pd.DataFrame(rows).sort_values(["f", "method", "metric"])
    strict = float(detail[detail.method.isin(("none", "joint", "harmony", "gene_rot"))]["abs_difference"].max())
    scvi_abs = float(detail[detail.method == "scvi"]["abs_difference"].max())
    return strict, scvi_abs, detail


def aggregate_gates(table: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    counts = table[["batch_type", "sigma", "f"]].drop_duplicates()
    if len(counts) != 30:
        raise RuntimeError(f"Expected 30 completed conditions, found {len(counts)}")
    repro, repro_scvi, reproduction = reproduction_check(table)
    reproduction.to_csv(OUT / "reproduction.csv", index=False)
    positive = table[table.f > 0].pivot(index=["batch_type", "sigma", "f"], columns="method", values=["mix_bg", "retention"])
    positive.columns = [f"{metric}_{method}" for metric, method in positive.columns]
    positive = positive.reset_index()
    if len(positive) != 24:
        raise RuntimeError("Expected 24 biology-positive conditions")
    for method in ("gene_rot", "harmony", "scvi"):
        positive[f"removes_{method}"] = positive[f"mix_bg_{method}"] >= 0.80
    positive["gain"] = positive["retention_gene_rot"] - np.maximum(positive["retention_harmony"], positive["retention_scvi"])
    positive["rigid_win"] = positive["removes_gene_rot"] & (positive["gain"] >= 0.10)
    zero = table[table.f == 0]
    false_max = {method: float(zero[zero.method == method]["bio_false"].max()) for method in ("gene_rot", "harmony", "scvi")}
    gates = {
        "n_conditions": int(len(counts)),
        "repro_t04_max_abs": repro,
        "repro_t04_scvi_abs": repro_scvi,
        "frac_rigid_removes_batch": float(positive["removes_gene_rot"].mean()),
        "frac_harmony_removes_batch": float(positive["removes_harmony"].mean()),
        "frac_scvi_removes_batch": float(positive["removes_scvi"].mean()),
        "frac_rigid_wins": float(positive["rigid_win"].mean()),
        "frac_rigid_wins_global": float(positive.loc[positive.batch_type == "global", "rigid_win"].mean()),
        "frac_rigid_wins_celltype": float(positive.loc[positive.batch_type == "celltype", "rigid_win"].mean()),
        "median_gain": float(positive["gain"].median()),
        "bio_false_max_gene_rot": false_max["gene_rot"],
        "bio_false_max_harmony": false_max["harmony"],
        "bio_false_max_scvi": false_max["scvi"],
    }
    return gates, positive


def sweep_stage() -> None:
    make_meta()
    if not (OUT / "data_meta.json").exists():
        raise RuntimeError("Run --stage data first")
    total = len(condition_rows())
    pending = []
    for spec in condition_rows():
        target = condition_result_path(*spec)
        if target.exists() and not checkpoint_is_complete(target, *spec):
            raise RuntimeError(f"Incomplete condition checkpoint requires review: {target}")
        if not checkpoint_is_complete(target, *spec):
            pending.append(spec)
    if pending:
        continuation = (json.loads(continuation_path().read_text()) if continuation_path().exists()
                        else create_parallel_continuation(pending))
        configure_spawn_environment()
        context = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=SWEEP_WORKERS, mp_context=context, initializer=init_sweep_worker) as pool:
            futures = [pool.submit(sweep_condition_worker, spec) for spec in pending]
            for future in as_completed(futures):
                result = future.result()
                done = len(complete_checkpoints())
                continuation["completion_order"].append(result)
                json_atomic_write(continuation_path(), continuation)
                print(json.dumps({"stage": "sweep", "event": "condition_complete", "completed_conditions": done,
                                  "total_conditions": total, "workers": SWEEP_WORKERS, **result},
                                 sort_keys=True), flush=True)
    table = collect_sweep()
    values, positive = aggregate_gates(table)
    update_gates(**values)
    wins = positive[["batch_type", "sigma", "f", "mix_bg_gene_rot", "mix_bg_harmony", "mix_bg_scvi",
                     "retention_gene_rot", "retention_harmony", "retention_scvi", "gain", "removes_gene_rot", "rigid_win"]]
    wins.to_csv(OUT / "win_conditions.csv", index=False)
    print(json.dumps({"stage": "sweep", **values}, sort_keys=True), flush=True)


def line_figure(table: pd.DataFrame, value: str, output: Path, ylabel: str) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    method_colors = {"none": "0.35", "joint": "C1", "harmony": "C2", "scvi": "C4", "gene_rot": "C3"}
    for row, batch_type in enumerate(BATCH_TYPES):
        for col, sigma in enumerate(SIGMAS):
            axis = axes[row, col]
            selected = table[(table.batch_type == batch_type) & (table.sigma == sigma)]
            if value == "retention":
                selected = selected[selected.f > 0]
            for method in METHODS:
                part = selected[selected.method == method].sort_values("f")
                axis.plot(part.f, part[value], marker="o", lw=1.8, label=method, color=method_colors[method])
            if value == "mix_bg":
                axis.axhline(0.8, color="0.25", ls=":", lw=1)
            axis.set_title(f"{batch_type}; σ={sigma}")
            axis.set_xlabel("biology size f")
            axis.set_ylabel(ylabel)
            axis.grid(alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(METHODS), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def wins_figure(positive: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharex=True, sharey=True)
    colors = {"win": "C3", "lose": "C1", "no_mix": "0.65"}
    for row, batch_type in enumerate(BATCH_TYPES):
        for col, sigma in enumerate(SIGMAS):
            axis = axes[row, col]
            sub = positive[(positive.batch_type == batch_type) & (positive.sigma == sigma)].sort_values("f")
            status = np.where(~sub["removes_gene_rot"], "no_mix", np.where(sub["rigid_win"], "win", "lose"))
            axis.scatter(sub.f, np.ones(len(sub)), s=420, c=[colors[value] for value in status], marker="s")
            for f, item, gain in zip(sub.f, status, sub.gain):
                axis.text(f, 1, f"f={int(f)}\n{item}\n{gain:.2f}", ha="center", va="center", fontsize=8)
            axis.set_title(f"{batch_type}; σ={sigma}")
            axis.set_xticks([1, 2, 4, 8])
            axis.set_yticks([])
            axis.set_xlim(0.3, 8.7)
    handles = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=colors[name], markersize=12, label=name)
               for name in ("win", "lose", "no_mix")]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(OUT / "fig_wins.png", dpi=180)
    plt.close(fig)


def report_stage() -> None:
    make_meta()
    write_requirements()
    if not (OUT / "sweep.csv").exists():
        raise RuntimeError("Run --stage sweep before --stage report")
    table = pd.read_csv(OUT / "sweep.csv")
    values, positive = aggregate_gates(table)
    if values["frac_rigid_wins"] >= 0.75 and values["bio_false_max_gene_rot"] <= 0.10:
        outcome = "general"
    elif values["frac_rigid_wins"] >= 0.25:
        outcome = "limited"
    else:
        outcome = "none"
    values["t04b_outcome"] = outcome
    update_gates(**values)
    line_figure(table, "retention", OUT / "fig_retention.png", "retention")
    line_figure(table, "mix_bg", OUT / "fig_mixing.png", "mix_bg")
    wins_figure(positive)

    files = ["requirements.lock", "sweep.csv", "fig_retention.png", "fig_mixing.png", "fig_wins.png", "scripts/t04b_run.py"]
    command = 'cd "$HOME/scopalign-cloud" && .venv/bin/python reports/T04b/scripts/t04b_run.py'
    regeneration = {path: command + " --stage all" for path in files}
    for name in values:
        regeneration[name] = command + (" --stage report" if name == "t04b_outcome" else " --stage sweep")
    artifacts = []
    for name in ("data_meta.json", "data_manifest.json", "reproduction.csv", "win_conditions.csv", "run_meta.json",
                 "sweep_continuation.json"):
        path = OUT / name
        if path.exists():
            artifacts.append({"path": name, "sha256": sha(path)})
    meta = make_meta()
    report = {
        "task_id": "T04b",
        "executor": "manus",
        "branch": "task/T04b",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": 20261014,
        "started": meta["started"],
        "finished": utcnow(),
        "status": "complete",
        "gated": values,
        "deliverables": [{"path": name, "sha256": sha(OUT / name)} for name in files],
        "artifacts": artifacts,
        "regeneration": regeneration,
        "deviations": [{
            "what": "Completed independent sweep conditions with two spawn-worker processes and three Torch CPU threads per worker, retaining atomic per-condition CSV checkpoints.",
            "why": "The serial runner used three of six available CPU cores and did not complete. Each worker preserves the pre-registered cached condition matrix, seeds, methods, 400-epoch CPU scVI fit, and scoring procedure; conditions have no cross-condition state."
        }],
        "deviations_requested": [],
        "observations": [
            "Cached 30 condition matrices and six batch-free oracle archives before the sweep.",
            "The two global sigma=0.5 reference matrices equal the corresponding T04 cached count matrices.",
            "Per-condition sweep checkpoints are stored under the gitignored cache directory.",
        ],
        "known_failures": T04.KNOWN,
        "lockfile_sha256": sha(OUT / "requirements.lock"),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"stage": "report", **values}, sort_keys=True), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("data", "sweep", "report", "all"))
    stage = parser.parse_args().stage
    stages = ("data", "sweep", "report") if stage == "all" else (stage,)
    for current in stages:
        if current == "data":
            data_stage()
        elif current == "sweep":
            sweep_stage()
        elif current == "report":
            report_stage()


if __name__ == "__main__":
    main()
