#!/usr/bin/env python
"""Pre-registered shared-detection positional-layer comparison for T00c."""
import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse
from sklearn.decomposition import FactorAnalysis, PCA

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "pbmc_v2_v3"
OUT = ROOT / "reports" / "T00c"
SEED = 20261014
SEEDS = (20261014, 20261015, 20261016)
D, K, FLOOR = 30, 40, 0.05


def qc_data(dataset):
    a = sc.read_h5ad(DATA / f"pbmc_{dataset}.h5ad")
    return a[a.obs["qc_pass"].to_numpy(dtype=bool)].copy()


def nonzero_counts(a):
    x = a.X
    return np.asarray((x > 0).sum(axis=0)).ravel()


def detected_genes(data):
    keep = {}
    for dataset, a in data.items():
        keep[dataset] = set(a.var_names[nonzero_counts(a) >= 3].astype(str))
    return sorted(keep["v2"] & keep["v3"])


def select_genes(write=False):
    data = {dataset: qc_data(dataset) for dataset in ("v2", "v3")}
    detected = detected_genes(data)
    if len(detected) < 2000:
        raise RuntimeError(f"Only {len(detected)} genes pass the detection filter")
    data = {dataset: a[:, detected].copy() for dataset, a in data.items()}
    both = ad.concat([data["v2"], data["v3"]], axis=0, join="inner", merge="same", index_unique=None)
    sc.pp.highly_variable_genes(both, flavor="seurat_v3", n_top_genes=2000, batch_key=None)
    hvg = sorted(both.var_names[both.var["highly_variable"].to_numpy()].astype(str))
    if len(hvg) != 2000:
        raise RuntimeError(f"Expected 2000 HVGs, got {len(hvg)}")
    if write:
        (OUT / "genes.txt").write_text("\n".join(hvg) + "\n")
    return data, detected, hvg


def load_inputs():
    genes = (OUT / "genes.txt").read_text().splitlines()
    if len(genes) != 2000 or genes != sorted(genes):
        raise RuntimeError("genes.txt must contain 2000 alphabetically sorted symbols")
    data = {dataset: qc_data(dataset) for dataset in ("v2", "v3")}
    return {dataset: a[:, genes].copy() for dataset, a in data.items()}


def norm(a):
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)


def ranks(z):
    z = np.asarray(z, dtype=np.float32)
    n, d = z.shape
    if d != D:
        raise RuntimeError(f"Expected d={D}, got {d}")
    out = np.zeros(n, dtype=np.int32)
    for c in range(n):
        distances = np.linalg.norm(z - z[c], axis=1)
        distances[c] = np.inf
        idx = np.argpartition(distances, K)[:K]
        ev = np.abs(np.linalg.eigvalsh(np.cov(z[idx], rowvar=False, bias=False)))
        out[c] = int((ev > FLOOR * float(ev.max())).sum())
    return out


def path(kind, embed, dataset, seed=None):
    suffix = f"_{seed}" if seed is not None else ""
    return OUT / f"{kind}_{embed}_{dataset}{suffix}.npy"


def fit(embed, dataset, seed):
    a = load_inputs()[dataset]
    suffix = seed if embed == "ldvae" else None
    noise = None
    if embed == "pca":
        if seed != SEED:
            raise ValueError(f"PCA seed is fixed at {SEED}")
        z = PCA(D, random_state=SEED).fit_transform(norm(a))
    elif embed == "fa":
        if seed != SEED:
            raise ValueError(f"FA seed is fixed at {SEED}")
        model = FactorAnalysis(D, random_state=SEED, max_iter=2000)
        z = model.fit_transform(norm(a))
        noise = model.noise_variance_
    elif embed == "ldvae":
        if seed not in SEEDS:
            raise ValueError(f"LDVAE seed must be one of {SEEDS}")
        scvi.settings.seed = seed
        scvi.model.LinearSCVI.setup_anndata(a, batch_key=None)
        model = scvi.model.LinearSCVI(a, n_latent=D)
        model.train(max_epochs=250, accelerator="cpu", devices=1)
        z = model.get_latent_representation()
    else:
        raise ValueError(embed)
    z = np.asarray(z, dtype=np.float32)
    r = ranks(z)
    np.save(path("latent", embed, dataset, suffix), z)
    np.save(path("rank", embed, dataset, suffix), r)
    if noise is not None:
        np.save(OUT / f"fa_noise_{dataset}.npy", np.asarray(noise, dtype=np.float64))
    print(json.dumps({"embed": embed, "dataset": dataset, "seed": seed, "n_cells": int(len(r)), "median_rank": float(np.median(r))}, sort_keys=True))


def rows():
    out = []
    for embed in ("pca", "fa"):
        for dataset in ("v2", "v3"):
            r = np.load(path("rank", embed, dataset))
            out.append(dict(embed=embed, dataset=dataset, seed=SEED, median_rank=float(np.median(r)), frac_at_ceiling=float((r >= 28).mean())))
    for seed in SEEDS:
        for dataset in ("v2", "v3"):
            r = np.load(path("rank", "ldvae", dataset, seed))
            out.append(dict(embed="ldvae", dataset=dataset, seed=seed, median_rank=float(np.median(r)), frac_at_ceiling=float((r >= 28).mean())))
    return out


def plot():
    bins = np.arange(-0.5, D + 1.5, 1)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for axis, embed in zip(axes, ("pca", "fa", "ldvae")):
        for dataset, color in (("v2", "C0"), ("v3", "C1")):
            values = np.concatenate([np.load(path("rank", embed, dataset, seed)) for seed in SEEDS]) if embed == "ldvae" else np.load(path("rank", embed, dataset))
            axis.hist(values, bins=bins, alpha=0.55, color=color, label=dataset)
        axis.set_title(embed.upper())
        axis.set_xlabel("tangent rank")
        axis.set_xticks(range(0, D + 1, 5))
    axes[0].set_ylabel("cells")
    axes[-1].legend(title="dataset")
    fig.tight_layout()
    fig.savefig(OUT / "fig_rank_dist.png", dpi=180)
    plt.close(fig)


def summarise():
    frame = pd.DataFrame(rows(), columns=["embed", "dataset", "seed", "median_rank", "frac_at_ceiling"])
    frame.to_csv(OUT / "rank_summary.csv", index=False)
    plot()
    data, detected, hvg = select_genes(write=False)
    med = frame.loc[frame.embed != "ldvae"].pivot(index="embed", columns="dataset", values="median_rank")
    fa = np.concatenate([np.load(path("rank", "fa", dataset)) for dataset in ("v2", "v3")])
    ldvae_20261014 = np.concatenate([np.load(path("rank", "ldvae", dataset, SEED)) for dataset in ("v2", "v3")])
    zero = {dataset: int((nonzero_counts(data[dataset][:, hvg]) == 0).sum()) for dataset in ("v2", "v3")}
    noise = {dataset: int((np.load(OUT / f"fa_noise_{dataset}.npy") < 1e-6).sum()) for dataset in ("v2", "v3")}
    gated = {
        "n_genes_detected_both": int(len(detected)),
        "n_hvg_zero_v2": zero["v2"], "n_hvg_zero_v3": zero["v3"],
        "n_fa_tiny_noise_v2": noise["v2"], "n_fa_tiny_noise_v3": noise["v3"],
        "rank_pca_v2": float(med.loc["pca", "v2"]), "rank_pca_v3": float(med.loc["pca", "v3"]),
        "rank_fa_v2": float(med.loc["fa", "v2"]), "rank_fa_v3": float(med.loc["fa", "v3"]),
        "rank_ldvae_v2": float(np.median(frame.loc[(frame.embed == "ldvae") & (frame.dataset == "v2"), "median_rank"])),
        "rank_ldvae_v3": float(np.median(frame.loc[(frame.embed == "ldvae") & (frame.dataset == "v3"), "median_rank"])),
        "frac_fa_at_ceiling": float((fa >= 28).mean()),
        "frac_ldvae_at_ceiling": float((ldvae_20261014 >= 28).mean()),
    }
    print(json.dumps(gated, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", action="store_true")
    parser.add_argument("--embed", choices=("pca", "fa", "ldvae"))
    parser.add_argument("--dataset", choices=("v2", "v3"))
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--summarise", action="store_true")
    args = parser.parse_args()
    if args.genes:
        if args.embed or args.dataset or args.summarise:
            parser.error("--genes cannot be combined with another action")
        _, detected, hvg = select_genes(write=True)
        print(json.dumps({"n_genes_detected_both": len(detected), "n_hvgs": len(hvg)}, sort_keys=True))
    elif args.summarise:
        if args.embed or args.dataset:
            parser.error("--summarise cannot be combined with --embed or --dataset")
        summarise()
    elif args.embed and args.dataset:
        fit(args.embed, args.dataset, args.seed)
    else:
        parser.error("choose --genes, --summarise, or --embed with --dataset")


if __name__ == "__main__":
    main()
