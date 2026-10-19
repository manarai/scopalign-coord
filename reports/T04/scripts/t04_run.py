#!/usr/bin/env python3
"""T04: semi-synthetic known batch and known biology benchmark."""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scipy.stats import ortho_group
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.compare import centroid_operator
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP
import scvi
import threadpoolctl

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "T04"
CACHE = OUT / "cache"
DATA = ROOT / "data" / "pbmc_10k_v3"
SOURCE = DATA / "pbmc_10k_v3.h5ad"
HALVES = DATA / "k0_halves.csv"
GENES_PATH = ROOT / "reports" / "K0" / "genes.txt"
D = 30
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
SCENARIOS = ("S0", "S1")
METHODS = ("none", "joint", "harmony", "scvi", "gene_rot", "scopalign")
ROTATION_METHODS = ("none", "gene_rot", "scopalign")
SCENARIO_REPLICATE = {"S0": 0, "S1": 1}
BIO_GENES = (
    "ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
    "RSAD2", "IFI44", "IFI44L", "HERC5", "STAT1", "IRF7", "XAF1", "LY6E", "EPSTI1", "IFITM3", "PLSCR1",
)
KNOWN = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
WORKER_PROBLEM = None


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def data_cache(scenario):
    return CACHE / "data" / f"{scenario}.npz"


def emb_cache(scenario):
    return CACHE / "embeddings" / f"{scenario}.npz"


def drift_path(scenario, side, seed):
    return CACHE / "drifts" / f"{scenario}_{side}_{seed}.pt"


def operator_path(scenario, side, seed):
    return CACHE / "operators" / f"{scenario}_{side}_{seed}.npy"


def gates_path():
    return OUT / "gates.json"


def requirements_path():
    return OUT / "requirements.lock"


def get_gates():
    return json.loads(gates_path().read_text()) if gates_path().exists() else {}


def put_gates(values):
    gates_path().write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")


def update_gates(**values):
    gate = get_gates()
    gate.update(values)
    put_gates(gate)
    return gate


def fingerprint(path):
    text = subprocess.check_output([sys.executable, "checker/fingerprint.py", str(path)], cwd=ROOT, text=True)
    return text.split()[0]


def write_requirements():
    subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, stdout=requirements_path().open("w"))


def make_meta():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    ignore = "cache/\n__pycache__/\ngates.json\nrun.log\nrun.pid\n"
    (OUT / ".gitignore").write_text(ignore)
    path = OUT / "run_meta.json"
    if not path.exists():
        path.write_text(json.dumps({
            "started": utcnow(),
            "platform": platform.platform(),
            "blas": threadpoolctl.threadpool_info(),
            "openblas_num_threads_env": os.environ.get("OPENBLAS_NUM_THREADS"),
            "omp_num_threads_env": os.environ.get("OMP_NUM_THREADS"),
            "default_openblas_settings": True,
        }, indent=2, sort_keys=True) + "\n")
    return json.loads(path.read_text())


def read_genes():
    genes = [line.strip() for line in GENES_PATH.read_text().splitlines() if line.strip()]
    if len(genes) != 2000 or genes != sorted(genes):
        raise RuntimeError("K0 genes.txt must contain exactly 2000 sorted genes")
    return genes


def source_halves():
    source = sc.read_h5ad(SOURCE)
    source = source[source.obs["qc_pass"].to_numpy(dtype=bool)].copy()
    halves = pd.read_csv(HALVES).set_index("barcode")["half"]
    names = source.obs_names.astype(str)
    a_mask = names.isin(halves.index) & (names.map(halves) == "H1")
    b_mask = names.isin(halves.index) & (names.map(halves) == "H2")
    a, b = source[a_mask].copy(), source[b_mask].copy()
    if a.n_obs != b.n_obs:
        raise RuntimeError("K0 halves are not equal sized")
    return a, b


def as_dense_counts(adata, genes):
    matrix = adata[:, genes].X
    matrix = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
    return np.asarray(matrix, dtype=np.int64)


def dense_log_counts(counts):
    work = ad.AnnData(X=sp.csr_matrix(np.asarray(counts, dtype=np.float32)))
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    return work.X.toarray() if sp.issparse(work.X) else np.asarray(work.X)


def load_case(scenario):
    with np.load(data_cache(scenario), allow_pickle=False) as arrays:
        return {key: arrays[key] for key in arrays.files}


def load_embedding(scenario):
    with np.load(emb_cache(scenario), allow_pickle=False) as arrays:
        return {key: arrays[key] for key in arrays.files}


def write_embedding(scenario, values):
    path = emb_cache(scenario)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **values)


def assert_active():
    gate = get_gates()
    if gate.get("stopped"):
        raise RuntimeError("T04 has a pre-registered stop; no later stage may run")


def data_stage():
    make_meta()
    genes = read_genes()
    a, b = source_halves()
    if any(gene not in a.var_names or gene not in b.var_names for gene in BIO_GENES):
        found = [gene for gene in BIO_GENES if gene in a.var_names and gene in b.var_names]
    else:
        found = list(BIO_GENES)
    a_counts = as_dense_counts(a, genes)
    b_counts = as_dense_counts(b, genes)
    types_a = a.obs["cell_type"].astype(str).to_numpy(dtype=str)
    types_b = b.obs["cell_type"].astype(str).to_numpy(dtype=str)
    mono_a = types_a == "Monocytes"
    mono_b = types_b == "Monocytes"
    b_mono_idx = np.flatnonzero(mono_b)
    if len(b_mono_idx) < 2:
        raise RuntimeError("No adequate B monocytes for specified responders")
    responder_idx = np.sort(np.random.default_rng(20261100).choice(b_mono_idx, len(b_mono_idx) // 2, replace=False))
    responder = np.zeros(len(b), dtype=bool)
    responder[responder_idx] = True
    bio_positions = np.asarray([genes.index(gene) for gene in BIO_GENES if gene in genes], dtype=np.int64)
    b_biology = b_counts.copy()
    if len(bio_positions) == len(BIO_GENES):
        mu = a_counts[mono_a][:, bio_positions].mean(axis=0, dtype=np.float64)
        depth_b = b_counts.sum(axis=1, dtype=np.float64)
        nbar = float(depth_b.mean())
        rates = 4.0 * (mu[None, :] + 0.1) * (depth_b[responder_idx, None] / nbar)
        b_biology[responder_idx[:, None], bio_positions[None, :]] += np.random.default_rng(20261102).poisson(rates).astype(np.int64)
    batch_rng = np.random.default_rng(20261101)
    batch_gene = np.exp(batch_rng.normal(0.0, 0.5, len(genes)))
    # The single generator follows the prescribed scenario order S0 then S1; A is drawn once after both B draws.
    b_s0 = batch_rng.poisson(b_counts.astype(np.float64) * batch_gene[None, :] * 0.7).astype(np.int64)
    b_s1 = batch_rng.poisson(b_biology.astype(np.float64) * batch_gene[None, :] * 0.7).astype(np.int64)
    a_resampled = batch_rng.poisson(a_counts).astype(np.int64)
    oracle_b = b_biology.astype(np.int64)
    shared = {
        "types_a": types_a, "types_b": types_b, "responder_b": responder,
        "barcodes_a": a.obs_names.astype(str).to_numpy(dtype=str), "barcodes_b": b.obs_names.astype(str).to_numpy(dtype=str),
    }
    for scenario, b_matrix in (("S0", b_s0), ("S1", b_s1)):
        payload = {"A": a_resampled, "B": b_matrix, **shared}
        path = data_cache(scenario)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **payload)
    oracle_path = CACHE / "data" / "oracle_S1.npz"
    np.savez_compressed(oracle_path, A=a_resampled, B=oracle_b, **shared)
    meta = {
        "source_sha256": sha(SOURCE), "source_fingerprint": fingerprint(SOURCE),
        "halves_sha256": sha(HALVES), "halves_fingerprint": fingerprint(HALVES),
        "genes_sha256": sha(GENES_PATH), "n_genes": len(genes), "n_bio_genes": len(bio_positions),
        "n_responders": int(responder.sum()), "n_cells_a": int(len(a)), "n_cells_b": int(len(b)),
        "n_monocytes_a": int(mono_a.sum()), "n_monocytes_b": int(mono_b.sum()),
        "responder_seed": 20261100, "biology_seed": 20261102, "batch_seed": 20261101,
        "batch_draw_order": ["batch_gene", "S0_B", "S1_B", "A_once"],
        "bio_genes": list(BIO_GENES),
    }
    (OUT / "data_meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    write_requirements()
    gate = update_gates(n_bio_genes=int(len(bio_positions)), n_responders=int(responder.sum()))
    if gate["n_bio_genes"] != 22:
        update_gates(stopped=True, stop_rule="n_bio_genes ≠ 22 → REVISE.")
    print(json.dumps({"n_bio_genes": gate["n_bio_genes"], "n_responders": gate["n_responders"]}, sort_keys=True), flush=True)


def raw_scvi_adata(a_counts, b_counts, genes):
    matrix = sp.vstack([sp.csr_matrix(a_counts), sp.csr_matrix(b_counts)]).astype(np.float32)
    obs = pd.DataFrame({"dataset": ["A"] * len(a_counts) + ["B"] * len(b_counts)})
    var = pd.DataFrame(index=pd.Index(genes, dtype=str))
    return ad.AnnData(X=matrix, obs=obs, var=var)


def load_oracle():
    with np.load(CACHE / "data" / "oracle_S1.npz", allow_pickle=False) as arrays:
        return {key: arrays[key] for key in arrays.files}


def tested_types(types, dataset, threshold=50):
    left, right = types[dataset == "A"], types[dataset == "B"]
    return sorted(kind for kind in set(left) & set(right) if np.sum(left == kind) >= threshold and np.sum(right == kind) >= threshold)


def mix_bg(rep, dataset, responder):
    keep = ~responder
    local_rep, local_ds = rep[keep], dataset[keep]
    neighbors = NearestNeighbors(n_neighbors=31).fit(local_rep).kneighbors(local_rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(local_ds[neighbors] != local_ds[:, None], axis=1)) / 0.5)


def bio_keep(rep, types, responder):
    monocytes = types == "Monocytes"
    local_rep, local_responder = rep[monocytes], responder[monocytes]
    if int(local_responder.sum()) == 0 or len(local_rep) < 31:
        return float("nan")
    neighbors = NearestNeighbors(n_neighbors=31).fit(local_rep).kneighbors(local_rep, return_distance=False)[:, 1:]
    query = np.flatnonzero(local_responder)
    p = float(np.mean(np.mean(local_responder[neighbors[query]], axis=1)))
    p0 = float(local_responder.sum() / len(local_responder))
    return float((p - p0) / (1.0 - p0))


def sil_score(rep, types, dataset):
    kinds = tested_types(types, dataset)
    selected = np.isin(types, kinds)
    return float(silhouette_score(rep[selected], types[selected], metric="euclidean"))


def score_rep(rep, dataset, types, responder):
    return {
        "mix_bg": mix_bg(rep, dataset, responder),
        "bio_score": bio_keep(rep, types, responder),
        "sil": sil_score(rep, types, dataset),
    }


def scores_table():
    path = OUT / "scores.csv"
    return pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["scenario", "method", "mix_bg", "bio_keep", "retention", "bio_false", "sil"])


def store_score_rows(rows):
    existing = scores_table()
    incoming = pd.DataFrame(rows)
    full = pd.concat([existing, incoming], ignore_index=True)
    full = full.drop_duplicates(["scenario", "method"], keep="last").sort_values(["scenario", "method"])
    full.to_csv(OUT / "scores.csv", index=False)


def embed_stage():
    make_meta()
    assert_active()
    genes = read_genes()
    oracle = load_oracle()
    oracle_dataset = np.asarray(["A"] * len(oracle["A"]) + ["B"] * len(oracle["B"]), dtype=str)
    oracle_types = np.concatenate([oracle["types_a"].astype(str), oracle["types_b"].astype(str)])
    oracle_responder = np.concatenate([np.zeros(len(oracle["A"]), dtype=bool), oracle["responder_b"].astype(bool)])
    oracle_log = np.vstack([dense_log_counts(oracle["A"]), dense_log_counts(oracle["B"])])
    oracle_rep = PCA(n_components=D, svd_solver="full", random_state=20261014).fit_transform(oracle_log).astype(np.float64)
    oracle_bio = bio_keep(oracle_rep, oracle_types, oracle_responder)
    all_rows = []
    for scenario in SCENARIOS:
        case = load_case(scenario)
        a_counts, b_counts = case["A"], case["B"]
        xa, xb = dense_log_counts(a_counts), dense_log_counts(b_counts)
        pca_a = PCA(n_components=D, svd_solver="full", random_state=20261014)
        pca_b = PCA(n_components=D, svd_solver="full", random_state=20261014)
        za = pca_a.fit_transform(xa).astype(np.float64)
        zb = pca_b.fit_transform(xb).astype(np.float64)
        stacked = np.vstack([xa, xb])
        joint = PCA(n_components=D, svd_solver="full", random_state=20261014).fit_transform(stacked).astype(np.float64)
        dataset = np.asarray(["A"] * len(za) + ["B"] * len(zb), dtype=str)
        meta = pd.DataFrame({"dataset": dataset})
        harmony = hm.run_harmony(joint, meta, "dataset", random_state=20261014, max_iter_harmony=20).Z_corr.T.astype(np.float64)
        scvi.settings.seed = 20261014
        scvi_data = raw_scvi_adata(a_counts, b_counts, genes)
        scvi.model.SCVI.setup_anndata(scvi_data, batch_key="dataset")
        scvi_model = scvi.model.SCVI(scvi_data, n_latent=D)
        scvi_model.train(max_epochs=400, accelerator="cpu", devices=1)
        scvi_latent = scvi_model.get_latent_representation().astype(np.float64)
        u, _, vt = np.linalg.svd(pca_b.components_.astype(np.float64) @ pca_a.components_.astype(np.float64).T, full_matrices=False)
        q_load = (u @ vt).astype(np.float64)
        rep = {
            "none": np.vstack([za, zb]).astype(np.float64), "joint": joint, "harmony": harmony,
            "scvi": scvi_latent, "gene_rot": np.vstack([za @ q_load.T, zb]).astype(np.float64),
        }
        types = np.concatenate([case["types_a"].astype(str), case["types_b"].astype(str)])
        responder = np.concatenate([np.zeros(len(za), dtype=bool), case["responder_b"].astype(bool)])
        values = {"z_a": za, "z_b": zb, "w_a": pca_a.components_.astype(np.float64), "w_b": pca_b.components_.astype(np.float64),
                  "q_load": q_load, "dataset": dataset, "types": types, "responder": responder, **rep}
        write_embedding(scenario, values)
        for method in ("none", "joint", "harmony", "scvi", "gene_rot"):
            score = score_rep(rep[method], dataset, types, responder)
            all_rows.append({"scenario": scenario, "method": method, "mix_bg": score["mix_bg"],
                             "bio_keep": score["bio_score"] if scenario == "S1" else np.nan,
                             "retention": score["bio_score"] / oracle_bio if scenario == "S1" else np.nan,
                             "bio_false": score["bio_score"] if scenario == "S0" else np.nan, "sil": score["sil"]})
    store_score_rows(all_rows)
    gate = update_gates(bio_keep_oracle=float(oracle_bio))
    if gate["bio_keep_oracle"] < 0.3:
        update_gates(stopped=True, stop_rule="bio_keep_oracle < 0.3 → REVISE.")
    print(json.dumps({"bio_keep_oracle": float(oracle_bio), "methods_completed": 5}, sort_keys=True), flush=True)


def drift_config(seed):
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=seed, batch_size=256,
                             loss="esm", lr=5e-4, beta=1.0, esm_laplacian="hutchinson", device="cpu")


def fit_worker(item):
    torch.set_num_threads(1)
    scenario, side, seed, z, target = item
    target = Path(target)
    if target.exists():
        saved = torch.load(target, map_location="cpu", weights_only=False)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        drift = EnergyDriftAdapter(D, config=drift_config(seed)).fit(np.asarray(z, dtype=np.float32), verbose=True)
        symmetry = drift.symmetry_check(np.asarray(z, dtype=np.float32), n_cells=100)
        saved = {"state_dict": drift._model.state_dict(), "config": asdict(drift.config),
                 "symmetry": {key: float(value) for key, value in symmetry.items()}}
        torch.save(saved, target)
    return {"scenario": scenario, "side": side, "seed": int(seed), "path": str(target.relative_to(ROOT)),
            "median_asym": float(saved["symmetry"]["median_asym"]), "max_asym": float(saved["symmetry"]["max_asym"])}


def fits_stage():
    make_meta()
    assert_active()
    todo, rows = [], []
    for scenario in SCENARIOS:
        embedding = load_embedding(scenario)
        for side, z in (("A", embedding["z_a"]), ("B", embedding["z_b"])):
            for seed in SEEDS:
                target = drift_path(scenario, side, seed)
                if target.exists():
                    saved = torch.load(target, map_location="cpu", weights_only=False)
                    rows.append({"scenario": scenario, "side": side, "seed": int(seed), "path": str(target.relative_to(ROOT)),
                                 "median_asym": float(saved["symmetry"]["median_asym"]), "max_asym": float(saved["symmetry"]["max_asym"])})
                else:
                    todo.append((scenario, side, seed, np.asarray(z, dtype=np.float32), str(target)))
    if todo:
        workers = max(1, (os.cpu_count() or 1) - 1)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows.extend(pool.map(fit_worker, todo))
    table = pd.DataFrame(rows).drop_duplicates(["scenario", "side", "seed"], keep="last").sort_values(["scenario", "side", "seed"])
    table.to_csv(OUT / "fit_symmetry.csv", index=False)
    symmetry_max = float(table["median_asym"].max())
    gate = update_gates(symmetry_max=symmetry_max)
    if symmetry_max > 1e-4:
        update_gates(stopped=True, stop_rule="symmetry_max > 1e-4 → REVISE.")
    print(json.dumps({"n_fits": int(len(table)), "symmetry_max": symmetry_max}, sort_keys=True), flush=True)


def build_model(saved):
    model = _EnergyMLP(D, int(saved["config"]["hidden"]), int(saved["config"]["depth"]))
    model.load_state_dict(saved["state_dict"])
    model = model.to(dtype=torch.float64, device="cpu").eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model


def models_for(scenario, side):
    return [build_model(torch.load(drift_path(scenario, side, seed), map_location="cpu", weights_only=False)) for seed in SEEDS]


def fields(model, values, need_j, track):
    z = values.detach().clone().requires_grad_(True)
    energy = model(z).reshape(-1)
    u = -torch.autograd.grad(energy.sum(), z, create_graph=(need_j or track), retain_graph=True)[0]
    if not need_j:
        return (u, None) if track else (u.detach(), None)
    rows = []
    for index in range(D):
        rows.append(torch.autograd.grad(u[:, index].sum(), z, create_graph=track, retain_graph=(track or index < D - 1))[0])
    jac = torch.stack(rows, dim=1)
    return (u, jac) if track else (u.detach(), jac.detach())


def ensemble(models, values, need_j, track):
    velocities, jacobians = [], []
    for model in models:
        velocity, jacobian = fields(model, values, need_j=need_j, track=track)
        velocities.append(velocity)
        if need_j:
            jacobians.append(jacobian)
    velocity = torch.stack(velocities).mean(dim=0)
    jacobian = torch.stack(jacobians).mean(dim=0) if need_j else None
    return velocity, jacobian


class E5Problem:
    def __init__(self, models_a, models_b, z_a, z_b, replicate):
        self.models_a, self.models_b = models_a, models_b
        self.a = torch.as_tensor(np.asarray(z_a, dtype=np.float64), dtype=torch.float64)
        self.b = torch.as_tensor(np.asarray(z_b, dtype=np.float64), dtype=torch.float64)
        rng = np.random.default_rng(20261040 + replicate)
        self.ia = rng.choice(len(z_a), min(256, len(z_a)), replace=False)
        self.ib = rng.choice(len(z_b), min(256, len(z_b)), replace=False)
        self.a_curv, self.b_curv = self.a[self.ia], self.b[self.ib]
        self.ua, _ = ensemble(self.models_a, self.a, need_j=False, track=False)
        self.ub, _ = ensemble(self.models_b, self.b, need_j=False, track=False)
        _, self.ja = ensemble(self.models_a, self.a_curv, need_j=True, track=False)
        _, self.jb = ensemble(self.models_b, self.b_curv, need_j=True, track=False)
        self.den_ua, self.den_ub = torch.sum(self.ua * self.ua), torch.sum(self.ub * self.ub)
        self.den_ja, self.den_jb = torch.sum(self.ja * self.ja), torch.sum(self.jb * self.jb)

    def loss(self, q, track):
        ub_qa, _ = ensemble(self.models_b, self.a @ q.T, need_j=False, track=track)
        ua_qtb, _ = ensemble(self.models_a, self.b @ q, need_j=False, track=track)
        d1 = torch.sum((ub_qa - self.ua @ q.T) ** 2) / self.den_ua
        d2 = torch.sum((ua_qtb - self.ub @ q) ** 2) / self.den_ub
        flow = (d1 + d2) / 2.0
        _, jb_qa = ensemble(self.models_b, self.a_curv @ q.T, need_j=True, track=track)
        _, ja_qtb = ensemble(self.models_a, self.b_curv @ q, need_j=True, track=track)
        qjaqt = torch.einsum("ij,njk,lk->nil", q, self.ja, q)
        qtjbq = torch.einsum("ji,njk,kl->nil", q, self.jb, q)
        h1 = torch.sum((jb_qa - qjaqt) ** 2) / self.den_ja
        h2 = torch.sum((ja_qtb - qtjbq) ** 2) / self.den_jb
        curvature = (h1 + h2) / 2.0
        return flow + curvature, flow, curvature


class E5Closure:
    def __init__(self, problem, q0, skew):
        self.problem, self.q0, self.skew = problem, q0, skew

    def __call__(self):
        self.skew.grad = None
        q = self.q0 @ torch.matrix_exp(self.skew - self.skew.T)
        value, _, _ = self.problem.loss(q, track=True)
        value.backward()
        return value


def starts():
    yield 0, np.eye(D, dtype=np.float64)
    for index in range(1, 20):
        yield index, np.asarray(ortho_group.rvs(D, random_state=20263000 + index), dtype=np.float64)


def init_start_worker(problem):
    global WORKER_PROBLEM
    torch.set_num_threads(1)
    WORKER_PROBLEM = problem


def optimize_start(item):
    torch.set_num_threads(1)
    start, q0_array = item
    q0 = torch.as_tensor(q0_array, dtype=torch.float64)
    skew = torch.zeros((D, D), dtype=torch.float64, requires_grad=True)
    optimizer = torch.optim.LBFGS([skew], lr=1.0, max_iter=500, tolerance_grad=1e-10,
                                  tolerance_change=1e-12, line_search_fn="strong_wolfe")
    optimizer.step(E5Closure(WORKER_PROBLEM, q0, skew))
    q = (q0 @ torch.matrix_exp(skew - skew.T)).detach()
    value, flow, curvature = WORKER_PROBLEM.loss(q, track=False)
    return {"start": int(start), "L_hat": float(value.detach()), "D": float(flow.detach()),
            "H": float(curvature.detach()), "Q": q.cpu().numpy()}


def estimate_e5(problem):
    workers = min(3, max(1, (os.cpu_count() or 1) - 1))
    with ProcessPoolExecutor(max_workers=workers, initializer=init_start_worker, initargs=(problem,)) as pool:
        result = list(pool.map(optimize_start, starts()))
    return min(result, key=lambda row: (row["L_hat"], row["start"])), result


def score_scopalign(scenario, q):
    embedding = load_embedding(scenario)
    rep = np.vstack([embedding["z_a"].astype(np.float64) @ q.T, embedding["z_b"].astype(np.float64)])
    values = {key: embedding[key] for key in embedding}
    values["scopalign"] = rep
    write_embedding(scenario, values)
    score = score_rep(rep, embedding["dataset"].astype(str), embedding["types"].astype(str), embedding["responder"].astype(bool))
    row = {"scenario": scenario, "method": "scopalign", "mix_bg": score["mix_bg"],
           "bio_keep": score["bio_score"] if scenario == "S1" else np.nan,
           "retention": score["bio_score"] / float(get_gates()["bio_keep_oracle"]) if scenario == "S1" else np.nan,
           "bio_false": score["bio_score"] if scenario == "S0" else np.nan, "sil": score["sil"]}
    store_score_rows([row])


def estimate_stage():
    make_meta()
    assert_active()
    q_values, start_rows = {}, []
    for scenario in SCENARIOS:
        embedding = load_embedding(scenario)
        problem = E5Problem(models_for(scenario, "A"), models_for(scenario, "B"), embedding["z_a"], embedding["z_b"], SCENARIO_REPLICATE[scenario])
        best, results = estimate_e5(problem)
        q_values[f"Q_{scenario}"] = best["Q"].astype(np.float64)
        for row in results:
            start_rows.append({"scenario": scenario, "start": row["start"], "L_hat": row["L_hat"], "D": row["D"], "H": row["H"], "selected": bool(row["start"] == best["start"])})
        score_scopalign(scenario, best["Q"])
    np.savez_compressed(OUT / "Q_hats.npz", **q_values)
    pd.DataFrame(start_rows).sort_values(["scenario", "start"]).to_csv(OUT / "estimate_starts.csv", index=False)
    print(json.dumps({"n_starts": int(len(start_rows)), "start_workers": int(min(3, max(1, (os.cpu_count() or 1) - 1)))}, sort_keys=True), flush=True)


def individual_operator(model, z):
    chunks = []
    for start in range(0, len(z), 128):
        values = torch.as_tensor(np.asarray(z[start:start + 128], dtype=np.float64), dtype=torch.float64)
        _, jacobian = fields(model, values, need_j=True, track=False)
        matrix = jacobian.cpu().numpy().astype(np.float64)
        chunks.append((matrix + np.swapaxes(matrix, 1, 2)) / 2.0)
    return np.concatenate(chunks, axis=0)


def operator(scenario, side, seed, z):
    path = operator_path(scenario, side, seed)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path, individual_operator(build_model(torch.load(drift_path(scenario, side, seed), map_location="cpu", weights_only=False)), z))
    return np.load(path)


def residual(left, right):
    diff = np.asarray(left, dtype=np.float64) - np.asarray(right, dtype=np.float64)
    denominator = float(np.sum(np.asarray(right, dtype=np.float64) ** 2))
    return float(np.sqrt(np.sum(diff ** 2) / denominator)) if denominator else float("inf")


def flag_groups(types_a, types_b, responder_b):
    groups = []
    for kind in sorted((set(types_a) & set(types_b)) - {"Monocytes"}):
        mask_a, mask_b = types_a == kind, types_b == kind
        if int(mask_a.sum()) >= 100 and int(mask_b.sum()) >= 100:
            groups.append((kind, mask_a, mask_b))
    mono_a, mono_b = types_a == "Monocytes", types_b == "Monocytes"
    for name, mask_b in (("Monocytes responders", mono_b & responder_b), ("Monocytes nonresponders", mono_b & ~responder_b)):
        if int(mono_a.sum()) >= 100 and int(mask_b.sum()) >= 100:
            groups.append((name, mono_a, mask_b))
    return groups


def centroid(z, matrices, mask):
    return centroid_operator(z[mask], matrices[mask], n_cells_centroid=30).astype(np.float64)


def group_null(z_a, z_b, matrices_a, matrices_b, mask_a, mask_b, rng):
    values = []
    for z, matrices, mask in ((z_a, matrices_a, mask_a), (z_b, matrices_b, mask_b)):
        indices = np.flatnonzero(mask)
        for seed in SEEDS[1:]:
            permutation = rng.permutation(indices)
            half = len(permutation) // 2
            first, second = permutation[:half], permutation[half:2 * half]
            left = centroid_operator(z[first], matrices[20261014][first], n_cells_centroid=30).astype(np.float64)
            right = centroid_operator(z[second], matrices[seed][second], n_cells_centroid=30).astype(np.float64)
            values.append(residual(left, right))
    return float(np.median(values))


def flag_stage():
    make_meta()
    assert_active()
    with np.load(OUT / "Q_hats.npz", allow_pickle=False) as values:
        qs = {scenario: values[f"Q_{scenario}"].astype(np.float64) for scenario in SCENARIOS}
    rows = []
    for scenario in SCENARIOS:
        embedding = load_embedding(scenario)
        z_a, z_b = embedding["z_a"].astype(np.float64), embedding["z_b"].astype(np.float64)
        types = embedding["types"].astype(str)
        responder = embedding["responder"].astype(bool)
        types_a, types_b = types[:len(z_a)], types[len(z_a):]
        responder_b = responder[len(z_a):]
        matrices = {side: {seed: operator(scenario, side, seed, z) for seed in SEEDS} for side, z in (("A", z_a), ("B", z_b))}
        groups = flag_groups(types_a, types_b, responder_b)
        if not groups:
            raise RuntimeError("No valid T04 operator flag groups")
        rng = np.random.default_rng(20261022)
        for name, mask_a, mask_b in groups:
            refit = group_null(z_a, z_b, matrices["A"], matrices["B"], mask_a, mask_b, rng)
            c_a = centroid(z_a, matrices["A"][20261014], mask_a)
            c_b = centroid(z_b, matrices["B"][20261014], mask_b)
            for method, q in (("none", np.eye(D)), ("gene_rot", embedding["q_load"].astype(np.float64)), ("scopalign", qs[scenario])):
                observed = residual(q @ c_a @ q.T, c_b)
                rows.append({"scenario": scenario, "method": method, "group": name, "observed_residual": observed,
                             "refit_null": refit, "flag_ratio": float(observed / refit) if refit else float("inf")})
    table = pd.DataFrame(rows).sort_values(["scenario", "method", "group"])
    table.to_csv(OUT / "flags.csv", index=False)
    print(json.dumps({"n_flag_rows": int(len(table)), "n_groups_s1": int((table.scenario == "S1").sum() // 3)}, sort_keys=True), flush=True)


def make_umap(rep):
    work = ad.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32))
    work.obsm["view"] = np.asarray(rep, dtype=np.float64)
    sc.pp.neighbors(work, n_neighbors=15, use_rep="view", random_state=20261014)
    sc.tl.umap(work, min_dist=0.5, random_state=20261014)
    return np.asarray(work.obsm["X_umap"], dtype=np.float64)


def view_stage(scenario):
    embedding = load_embedding(scenario)
    scores = scores_table().set_index(["scenario", "method"])
    dataset, types, responder = embedding["dataset"].astype(str), embedding["types"].astype(str), embedding["responder"].astype(bool)
    reps = {method: embedding[method].astype(np.float64) for method in METHODS}
    fig, axes = plt.subplots(len(METHODS), 3, figsize=(15, 4 * len(METHODS)))
    palette = plt.get_cmap("tab10")
    for row, method in enumerate(METHODS):
        xy = make_umap(reps[method])
        score = scores.loc[(scenario, method)]
        title_score = f"mix={float(score.mix_bg):.3f}, " + (f"retention={float(score.retention):.3f}" if scenario == "S1" else f"bio_false={float(score.bio_false):.3f}")
        for label, color in (("A", "C0"), ("B", "C1")):
            mask = dataset == label
            axes[row, 0].scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=0.7, color=color, label=label, linewidths=0)
        axes[row, 0].set_title(f"{method}: half; {title_score}")
        axes[row, 0].legend(markerscale=2, fontsize=7)
        for index, kind in enumerate(sorted(set(types))):
            mask = types == kind
            axes[row, 1].scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=0.7, color=palette(index % 10), label=kind, linewidths=0)
        axes[row, 1].set_title(f"{method}: cell type; {title_score}")
        axes[row, 1].legend(markerscale=2, fontsize=6, ncol=2)
        axes[row, 2].scatter(xy[:, 0], xy[:, 1], s=3, alpha=0.25, color="0.65", linewidths=0)
        axes[row, 2].scatter(xy[responder, 0], xy[responder, 1], s=6, alpha=0.85, color="crimson", label="responders" if scenario == "S1" else "pseudo-responders", linewidths=0)
        axes[row, 2].set_title(f"{method}: responder status; {title_score}")
        axes[row, 2].legend(markerscale=2, fontsize=7)
        for column in range(3):
            axes[row, column].set_xticks([])
            axes[row, column].set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / f"fig_umap_{scenario}.png", dpi=180)
    plt.close(fig)


def flags_figure():
    flags = pd.read_csv(OUT / "flags.csv")
    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=True)
    for axis, scenario in zip(axes, SCENARIOS):
        table = flags[flags.scenario == scenario]
        groups = sorted(table.group.unique())
        x = np.arange(len(groups))
        width = 0.25
        for offset, method in enumerate(ROTATION_METHODS):
            values = [float(table[(table.method == method) & (table.group == group)].flag_ratio.iloc[0]) for group in groups]
            axis.bar(x + (offset - 1) * width, values, width=width, label=method)
        axis.axhline(1.5, color="0.3", ls="--", lw=1)
        axis.axhline(2.0, color="0.3", ls=":", lw=1)
        axis.set_xticks(x)
        axis.set_xticklabels(groups, rotation=30, ha="right")
        axis.set_title(scenario)
        axis.set_ylabel("flag ratio")
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_flags.png", dpi=180)
    plt.close(fig)


def outcome(gated):
    methods = METHODS
    removes = {method: gated[f"mix_bg_{method}"] >= 0.80 for method in methods}
    keeps = {method: gated[f"retention_{method}"] >= 0.8 for method in methods}
    no_false = {method: gated[f"bio_false_{method}"] <= 0.1 for method in methods}
    flags = {method: gated[f"flag_resp_{method}"] >= 2.0 and gated[f"flag_other_max_{method}"] <= 1.5 and gated[f"flag_s0_max_{method}"] <= 1.5 for method in ROTATION_METHODS}
    if removes["scopalign"] and keeps["scopalign"] and no_false["scopalign"] and ((removes["harmony"] and not keeps["harmony"]) or (removes["scvi"] and not keeps["scvi"]) or flags["scopalign"]):
        return "supported"
    if (not keeps["scopalign"]) or ((gated["mix_bg_scopalign"] < gated["mix_bg_harmony"] - 0.1 and gated["mix_bg_scopalign"] < gated["mix_bg_scvi"] - 0.1) and not flags["scopalign"]):
        return "not_supported"
    return "mixed"


def report_stage():
    make_meta()
    gate = get_gates()
    write_requirements()
    stopped = bool(gate.get("stopped"))
    if stopped:
        gated = {key: value for key, value in gate.items() if key in {"n_bio_genes", "n_responders", "bio_keep_oracle", "symmetry_max"}}
    else:
        score = scores_table().set_index(["scenario", "method"])
        flags = pd.read_csv(OUT / "flags.csv")
        gated = {"n_bio_genes": int(gate["n_bio_genes"]), "n_responders": int(gate["n_responders"]), "bio_keep_oracle": float(gate["bio_keep_oracle"])}
        for method in METHODS:
            gated[f"mix_bg_{method}"] = float(score.loc[("S1", method), "mix_bg"])
            gated[f"retention_{method}"] = float(score.loc[("S1", method), "retention"])
            gated[f"bio_false_{method}"] = float(score.loc[("S0", method), "bio_false"])
        for method in ROTATION_METHODS:
            s1 = flags[(flags.scenario == "S1") & (flags.method == method)]
            s0 = flags[(flags.scenario == "S0") & (flags.method == method)]
            resp = float(s1.loc[s1.group == "Monocytes responders", "flag_ratio"].iloc[0])
            other = float(s1.loc[s1.group != "Monocytes responders", "flag_ratio"].max())
            gated[f"flag_resp_{method}"] = resp
            gated[f"flag_other_max_{method}"] = other
            gated[f"flag_s0_max_{method}"] = float(s0.flag_ratio.max())
        gated["symmetry_max"] = float(gate["symmetry_max"])
        gated["t04_outcome"] = outcome(gated)
        view_stage("S1")
        view_stage("S0")
        flags_figure()
    files = ["requirements.lock", "data_meta.json", "scores.csv", "flags.csv", "fig_umap_S1.png", "fig_umap_S0.png", "fig_flags.png", "scripts/t04_run.py"]
    if stopped:
        files = [path for path in files if (OUT / path).exists()]
    command = 'cd "$HOME/scopalign-cloud" && .venv/bin/python reports/T04/scripts/t04_run.py'
    regeneration = {path: command + " --stage all" for path in files}
    stage_for = {"n_bio_genes": "data", "n_responders": "data", "bio_keep_oracle": "embed", "symmetry_max": "fits", "t04_outcome": "report"}
    for key in gated:
        if key.startswith(("mix_bg_", "retention_", "bio_false_")):
            stage_for[key] = "estimate" if key.endswith("scopalign") else "embed"
        elif key.startswith("flag_"):
            stage_for[key] = "flag"
        regeneration[key] = command + " --stage " + stage_for[key]
    meta = make_meta()
    report = {
        "task_id": "T04", "executor": "manus", "branch": "task/T04",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": 20261014, "started": meta["started"], "finished": utcnow(), "status": "stopped" if stopped else "complete",
        "gated": gated, "deliverables": [{"path": path, "sha256": sha(OUT / path)} for path in files],
        "artifacts": [{"path": "Q_hats.npz", "sha256": sha(OUT / "Q_hats.npz")}, {"path": "run_meta.json", "sha256": sha(OUT / "run_meta.json")},
                      {"path": "fit_symmetry.csv", "sha256": sha(OUT / "fit_symmetry.csv")}, {"path": "estimate_starts.csv", "sha256": sha(OUT / "estimate_starts.csv")}],
        "regeneration": regeneration, "deviations": [], "deviations_requested": [],
        "observations": [f"Input SHA-256: {path.relative_to(ROOT)} = {sha(path)}." for path in (SOURCE, HALVES, GENES_PATH)],
        "known_failures": KNOWN, "lockfile_sha256": sha(requirements_path()),
    }
    if stopped:
        report["stop_rule"] = str(gate["stop_rule"])
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(gated, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("data", "embed", "fits", "estimate", "flag", "report", "all"))
    stage = parser.parse_args().stage
    stages = ("data", "embed", "fits", "estimate", "flag", "report") if stage == "all" else (stage,)
    for item in stages:
        if item == "data":
            data_stage()
        elif item == "embed":
            embed_stage()
        elif item == "fits":
            fits_stage()
        elif item == "estimate":
            estimate_stage()
        elif item == "flag":
            flag_stage()
        elif item == "report":
            report_stage()
        if get_gates().get("stopped") and item != "report":
            if stage == "all":
                report_stage()
            return


if __name__ == "__main__":
    main()
