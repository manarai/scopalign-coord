#!/usr/bin/env python
"""Stage the pre-registered 10x PBMC v2/v3 count matrices for T00."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path(__file__).resolve().parents[3]
REPORT = ROOT / "reports" / "T00"
DATA = ROOT / "data" / "pbmc_v2_v3"
RAW = DATA / "raw"
URLS = {
    "v2": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5",
    "v3": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.stat().st_size > 0:
        return
    partial = output.with_suffix(output.suffix + ".partial")
    subprocess.run([
        "curl", "--fail", "--location", "--retry", "5", "--retry-all-errors",
        "--continue-at", "-", "--output", str(partial), url,
    ], check=True)
    partial.replace(output)


def raw_path(chemistry: str) -> Path:
    return RAW / Path(URLS[chemistry]).name


def write_checksums() -> None:
    lines = []
    for chemistry in ("v2", "v3"):
        path = raw_path(chemistry)
        if not path.exists():
            raise FileNotFoundError(path)
        lines.append(f"{sha256(path)}  {path.relative_to(ROOT)}")
    (REPORT / "data_checksums.txt").write_text("\n".join(lines) + "\n")


def integer_counts(adata, chemistry: str) -> None:
    matrix = adata.X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.allclose(values, np.round(values)):
        raise RuntimeError(f"{chemistry}: 10x matrix contains non-integer raw counts")
    adata.X = matrix.astype(np.int32)


def celltypist_labels(adata) -> np.ndarray:
    # Keep model-cache writes within the task data directory.
    model_home = DATA / "celltypist_home"
    model_home.mkdir(parents=True, exist_ok=True)
    os.environ["HOME"] = str(model_home)
    import celltypist

    celltypist.models.download_models(model="Immune_All_High.pkl")
    label_input = adata.copy()
    sc.pp.normalize_total(label_input, target_sum=1e4)
    sc.pp.log1p(label_input)
    prediction = celltypist.annotate(
        label_input, model="Immune_All_High.pkl", majority_voting=False
    )
    labels = prediction.predicted_labels
    if isinstance(labels, pd.DataFrame):
        labels = labels["predicted_labels"]
    labels = pd.Series(labels, index=labels.index, dtype="string").reindex(adata.obs_names)
    if labels.isna().any():
        raise RuntimeError("CellTypist returned missing predicted labels")
    return labels.astype(str).to_numpy()


def build(chemistry: str) -> None:
    adata = sc.read_10x_h5(raw_path(chemistry), genome=None, gex_only=True)
    if "gene_ids" not in adata.var.columns:
        raise RuntimeError(f"{chemistry}: read_10x_h5 did not retain var['gene_ids']")
    adata.var_names_make_unique()
    integer_counts(adata, chemistry)
    adata.var["mt"] = np.asarray(adata.var_names.astype(str).str.startswith("MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata.obs["qc_pass"] = (
        (adata.obs["n_genes_by_counts"] >= 200)
        & (adata.obs["pct_counts_mt"] < 20)
    ).astype(bool)
    adata.obs["cell_type"] = celltypist_labels(adata)
    adata.write_h5ad(DATA / f"pbmc_{chemistry}.h5ad", compression="gzip")


def prepare() -> None:
    for chemistry, url in URLS.items():
        download(url, raw_path(chemistry))
    write_checksums()
    for chemistry in URLS:
        build(chemistry)
    (REPORT / "label_map.csv").write_text("label_v2,label_v3,label\n")


def count() -> None:
    result = {}
    for chemistry in ("v2", "v3"):
        adata = sc.read_h5ad(DATA / f"pbmc_{chemistry}.h5ad")
        result[f"n_cells_{chemistry}_raw"] = int(adata.n_obs)
    print(json.dumps(result, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", action="store_true")
    args = parser.parse_args()
    if args.count:
        count()
    else:
        prepare()


if __name__ == "__main__":
    main()
