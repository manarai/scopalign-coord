#!/usr/bin/env python3
"""C0 cloud bootstrap validation: data contents, PCA, and one K2 drift fit."""
import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import distribution
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "C0"
DATA = ROOT / "data"
K0 = ROOT / "reports" / "K0"
K2 = ROOT / "reports" / "K2"
KNOWN = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
RAW = {
    DATA / "pbmc_v2_v3/raw/pbmc_1k_v2_filtered_feature_bc_matrix.h5": "344673d4eca5af4b4f8343ef85c039fa5b24fea75a7d714a54427548b7b60213",
    DATA / "pbmc_v2_v3/raw/pbmc_1k_v3_filtered_feature_bc_matrix.h5": "8191f576550c1b449d03441b9eb098ee9f73fa82513d171ba87d31d551e3ffda",
    DATA / "pbmc_10k_v3/raw/pbmc_10k_v3_filtered_feature_bc_matrix.h5": "ebc5dedc938830e20f8e1aafb893b0a9e0bf88584f3ef2d00b232dd277a188af",
}
FINGERPRINTS = {
    DATA / "pbmc_v2_v3/pbmc_v2.h5ad": "a6af64ce541a5c1d8e6af594c6ca7e421af8725f2e1e36cabea407dd720efb10",
    DATA / "pbmc_v2_v3/pbmc_v3.h5ad": "bcb3f6cd14e807e2367838c58142ef3d2aa341f985a1a5bc8c88b55b75656d62",
    DATA / "pbmc_v2_v3/pbmc_v3_thin.h5ad": "83679e2d9eaadb7f1834d7469de704234dccb9e8d7136de13e6014be307084ba",
    DATA / "pbmc_v2_v3/k2_cells.csv": "d59d4e47489f1f50e4839e7d9d25f8e827eaa45a99b710598f5e0396de132a07",
    DATA / "pbmc_10k_v3/pbmc_10k_v3.h5ad": "fd80192876343811a32822d8e560f0eff453e91dad936890374313e3405c4e3d",
    DATA / "pbmc_10k_v3/A_matched.h5ad": "3e361683ba93d2d60c511a72de98869ba03c902d8ec2cdfda36d9a39aff2263b",
    DATA / "pbmc_10k_v3/B_matched.h5ad": "94ffda9b16dd9bc7e764c36d7472db2226a97b21e487242863157f68943596c2",
    DATA / "pbmc_10k_v3/k0_halves.csv": "41fdb7ba192b953938ac62191ce09fdccb3b2733f6c5340e251fd439dba0fd16",
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path):
    value = subprocess.check_output([sys.executable, "checker/fingerprint.py", str(path)], cwd=ROOT, text=True)
    return value.split()[0]


def dense_norm(adata, genes):
    missing = [gene for gene in genes if gene not in adata.var_names]
    if missing:
        raise RuntimeError(f"K0 PCA check missing {len(missing)} requested genes")
    work = adata[:, genes].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    return work.X.toarray() if hasattr(work.X, "toarray") else np.asarray(work.X)


def k0_latent_check():
    source = ad.read_h5ad(DATA / "pbmc_10k_v3/pbmc_10k_v3.h5ad")
    source = source[source.obs["qc_pass"].to_numpy(dtype=bool)].copy()
    halves = pd.read_csv(DATA / "pbmc_10k_v3/k0_halves.csv").set_index("barcode")["half"]
    names = source.obs_names.astype(str)
    h1 = source[names.isin(halves.index) & (names.map(halves) == "H1")].copy()
    h2 = source[names.isin(halves.index) & (names.map(halves) == "H2")].copy()
    genes = [item.strip() for item in (K0 / "genes.txt").read_text().splitlines() if item.strip()]
    if len(genes) != 2000:
        raise RuntimeError("K0 genes.txt must contain 2000 genes")
    p1 = PCA(n_components=30, svd_solver="full", random_state=20261014)
    p2 = PCA(n_components=30, svd_solver="full", random_state=20261014)
    z1 = p1.fit_transform(dense_norm(h1, genes)).astype(np.float32)
    z2 = p2.fit_transform(dense_norm(h2, genes)).astype(np.float32)
    c1 = np.load(K0 / "latent_pca30_h1.npy").astype(np.float32)
    c2 = np.load(K0 / "latent_pca30_h2.npy").astype(np.float32)
    return float(max(np.max(np.abs(z1 - c1)), np.max(np.abs(z2 - c2))))


def k2_cfg(seed):
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=seed, batch_size=256, loss="esm", lr=5e-4, beta=1.0, esm_laplacian="hutchinson", device="cpu")


def k2_tangent_basis(z):
    neighbors = NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    basis = np.empty((len(z), 30, 10), dtype=np.float64)
    for index in range(len(z)):
        _, vectors = np.linalg.eigh(np.cov(z[neighbors[index]], rowvar=False))
        basis[index] = vectors[:, -10:]
    return basis


def k2_drift_repro_check():
    latent = np.load(K2 / "latent_pca30_a.npy").astype(np.float32)
    basis = k2_tangent_basis(latent)
    target = np.load(K2 / "ops_a.npz", allow_pickle=False)["lam_seeds"][0].astype(np.float64)
    cache = OUT / "cloud_k2_drift.npz"
    if cache.is_file():
        lam = np.load(cache, allow_pickle=False)["lam_seeds"].astype(np.float64)
        return float(np.max(np.abs(lam - target)))
    drift = EnergyDriftAdapter(30, config=k2_cfg(20261014)).fit(latent, verbose=True)
    jacobian = drift.jacobian_autodiff(latent).astype(np.float64)
    jacobian = (jacobian + np.swapaxes(jacobian, 1, 2)) / 2.0
    restricted = np.einsum("ndi,ndm,nmj->nij", basis, jacobian, basis)
    restricted = (restricted + np.swapaxes(restricted, 1, 2)) / 2.0
    lam = np.linalg.eigvalsh(restricted).astype(np.float64)
    np.savez_compressed(cache, lam_seeds=lam, config=json.dumps(asdict(drift.config), sort_keys=True))
    return float(np.max(np.abs(lam - target)))


def installed_scjdo2_sha():
    direct = Path(distribution("scjdo2")._path) / "direct_url.json"
    if direct.is_file():
        obj = json.loads(direct.read_text())
        return str(obj.get("vcs_info", {}).get("commit_id", ""))
    return ""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    raw_rows = []
    for path, expected in RAW.items():
        got = sha(path)
        raw_rows.append((path.relative_to(ROOT).as_posix(), got, expected, int(got == expected)))
    fp_rows = []
    for path, expected in FINGERPRINTS.items():
        got = fingerprint(path)
        fp_rows.append((path.relative_to(ROOT).as_posix(), got, expected, int(got == expected)))
    with (OUT / "fingerprints.txt").open("w") as handle:
        for rel, got, expected, match in fp_rows:
            handle.write(f"{got}\t{rel}\texpected={expected}\tmatch={match}\n")
    latent_max = k0_latent_check()
    drift_max = k2_drift_repro_check()
    scjdo_sha = installed_scjdo2_sha()
    gated = {
        "n_raw_sha_match": int(sum(row[3] for row in raw_rows)),
        "n_fingerprints_match": int(sum(row[3] for row in fp_rows)),
        "latent_check_max": float(latent_max),
        "drift_repro_max_abs": float(drift_max),
        "scjdo2_sha": scjdo_sha,
        "platform": platform.platform(),
        "torch_version": torch.__version__,
    }
    lock = OUT / "requirements.lock"
    subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=lock.open("w"), check=True)
    script = "scripts/c0_check.py"
    base = 'cd "$HOME/scopalign-cloud" && .venv/bin/python reports/C0/scripts/c0_check.py'
    report = {
        "task_id": "C0", "executor": "manus", "branch": "task/C0",
        "commit": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "seed": 20261014, "started": started, "finished": datetime.now(timezone.utc).isoformat(), "status": "complete",
        "gated": gated,
        "deliverables": [
            {"path": "requirements.lock", "sha256": sha(lock)},
            {"path": "fingerprints.txt", "sha256": sha(OUT / "fingerprints.txt")},
            {"path": script, "sha256": sha(OUT / script)},
            {"path": "report.json", "sha256": "0" * 64},
        ],
        "regeneration": {name: base for name in gated} | {"requirements.lock": base, "fingerprints.txt": base, script: base, "report.json": base},
        "deviations": [
            {"what": "Installed CPython 3.11.16 with uv.", "why": "The base Ubuntu image supplied Python 3.12 only and its repositories did not provide python3.11; C0 requires a Python 3.11 venv."},
            {"what": "Installed the exact non-scJDO lock entries with --no-deps.", "why": "The pinned lock contains an internal resolver conflict: dynamo-release==1.5.3 declares anndata<0.11 while the same lock pins anndata==0.12.19; --no-deps retained every pinned lock version."},
            {"what": "Used torch==2.14.1+cpu instead of the CUDA-linked PyPI torch==2.14.1 Linux wheel.", "why": "The pinned PyPI Linux wheel could not import without libcudart.so.13 in this CPU-only cloud sandbox; the matching CPU build retains base version 2.14.1."},
            {"what": "Ran reports/T00/scripts/t00_stage.py without --stage arguments.", "why": "The committed T00 script exposes only its default prepare command, which performs both its requested raw downloads and h5ad builds."},
        ],
        "deviations_requested": [],
        "observations": [
            "All raw and content fingerprint comparisons are written to the generated outputs.",
            "cloud_k2_drift.npz is a non-deliverable local diagnostic produced by the cross-platform drift check.",
        ],
        "known_failures": KNOWN, "lockfile_sha256": sha(lock),
        "artifacts": [{"path": "cloud_k2_drift.npz", "sha256": sha(OUT / "cloud_k2_drift.npz")}],
    }
    # report.json cannot checksum itself before it is written; omit self-reference from deliverables.
    report["deliverables"] = report["deliverables"][:-1]
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(gated, sort_keys=True))


if __name__ == "__main__":
    main()
