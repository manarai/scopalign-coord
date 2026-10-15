#!/usr/bin/env python
"""Pre-registered PCA, FA, and LDVAE tangent-rank calculation for T00b."""
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
OUT = ROOT / "reports" / "T00b"
SEED = 20261014
SEEDS = (20261014, 20261015, 20261016)
D, K, FLOOR = 30, 40, 0.05


def load_inputs():
    data = {}
    for ds in ("v2", "v3"):
        a = sc.read_h5ad(DATA / f"pbmc_{ds}.h5ad")
        data[ds] = a[a.obs["qc_pass"].to_numpy(dtype=bool)].copy()
    genes = np.intersect1d(data["v2"].var_names.astype(str), data["v3"].var_names.astype(str))
    if len(genes) < 2000:
        raise RuntimeError("Fewer than 2000 intersected genes")
    data = {ds: a[:, genes].copy() for ds, a in data.items()}
    both = ad.concat([data["v2"], data["v3"]], axis=0, join="inner", merge="same", index_unique=None)
    sc.pp.highly_variable_genes(both, n_top_genes=2000, flavor="seurat_v3", batch_key=None)
    hvg = both.var_names[both.var["highly_variable"].to_numpy()]
    if len(hvg) != 2000:
        raise RuntimeError(f"Expected 2000 HVGs, got {len(hvg)}")
    return {ds: a[:, hvg].copy() for ds, a in data.items()}


def log_normalised(a):
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)


def ranks(z):
    z = np.asarray(z, dtype=np.float32)
    n, d = z.shape
    if d != D:
        raise RuntimeError(f"Expected d={D}, got {d}")
    k = min(K, n - 1)
    out = np.zeros(n, dtype=np.int32)
    for c in range(n):
        idx = np.argpartition(np.linalg.norm(z - z[c], axis=1), k)[:k]
        sigma = np.cov(z[idx], rowvar=False, bias=False) + 1e-6 * np.eye(d)
        ev = np.abs(np.linalg.eigvalsh(sigma))
        out[c] = (ev > FLOOR * float(ev.max())).sum()
    return out


def path(kind, embed, dataset, seed=None):
    suffix = f"_{seed}" if seed is not None else ""
    return OUT / f"{kind}_{embed}_{dataset}{suffix}.npy"


def fit(embed, dataset, seed):
    a = load_inputs()[dataset]
    suffix = seed if embed == "ldvae" else None
    if embed == "pca":
        if seed != SEED:
            raise ValueError(f"PCA seed is fixed at {SEED}")
        z = PCA(D, random_state=SEED).fit_transform(log_normalised(a))
    elif embed == "fa":
        if seed != SEED:
            raise ValueError(f"FA seed is fixed at {SEED}")
        z = FactorAnalysis(D, random_state=SEED, max_iter=2000).fit_transform(log_normalised(a))
    elif embed == "ldvae":
        if seed not in SEEDS:
            raise ValueError(f"LDVAE seed must be one of {SEEDS}")
        scvi.settings.seed = seed
        scvi.model.LinearSCVI.setup_anndata(a, batch_key=None)
        m = scvi.model.LinearSCVI(a, n_latent=D)
        m.train(max_epochs=250, accelerator="cpu", devices=1)
        z = m.get_latent_representation()
    else:
        raise ValueError(embed)
    z = np.asarray(z, dtype=np.float32)
    r = ranks(z)
    np.save(path("latent", embed, dataset, suffix), z)
    np.save(path("rank", embed, dataset, suffix), r)
    print(json.dumps({"embed": embed, "dataset": dataset, "seed": seed, "n_cells": int(len(r)), "median_rank": float(np.median(r))}, sort_keys=True))


def rows():
    out = []
    for embed in ("pca", "fa"):
        for ds in ("v2", "v3"):
            r = np.load(path("rank", embed, ds))
            out.append(dict(embed=embed, dataset=ds, seed=SEED, median_rank=float(np.median(r)), frac_at_ceiling=float((r >= 28).mean())))
    for seed in SEEDS:
        for ds in ("v2", "v3"):
            r = np.load(path("rank", "ldvae", ds, seed))
            out.append(dict(embed="ldvae", dataset=ds, seed=seed, median_rank=float(np.median(r)), frac_at_ceiling=float((r >= 28).mean())))
    return out


def plot():
    bins = np.arange(-0.5, D + 1.5, 1)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for ax, embed in zip(axs, ("pca", "fa", "ldvae")):
        for ds, color in (("v2", "C0"), ("v3", "C1")):
            values = np.concatenate([np.load(path("rank", embed, ds, s)) for s in SEEDS]) if embed == "ldvae" else np.load(path("rank", embed, ds))
            ax.hist(values, bins=bins, alpha=.55, color=color, label=ds)
        ax.set_title(embed.upper())
        ax.set_xlabel("tangent rank")
        ax.set_xticks(range(0, D + 1, 5))
    axs[0].set_ylabel("cells")
    axs[-1].legend(title="dataset")
    fig.tight_layout()
    fig.savefig(OUT / "fig_rank_dist.png", dpi=180)
    plt.close(fig)


def summarise():
    frame = pd.DataFrame(rows(), columns=["embed", "dataset", "seed", "median_rank", "frac_at_ceiling"])
    frame.to_csv(OUT / "rank_summary.csv", index=False)
    plot()
    med = frame.loc[frame.embed != "ldvae"].pivot(index="embed", columns="dataset", values="median_rank")
    fa = np.concatenate([np.load(path("rank", "fa", ds)) for ds in ("v2", "v3")])
    gated = {
        "rank_pca_v2": float(med.loc["pca", "v2"]), "rank_pca_v3": float(med.loc["pca", "v3"]),
        "rank_fa_v2": float(med.loc["fa", "v2"]), "rank_fa_v3": float(med.loc["fa", "v3"]),
        "rank_ldvae_v2": float(np.median(frame.loc[(frame.embed == "ldvae") & (frame.dataset == "v2"), "median_rank"])),
        "rank_ldvae_v3": float(np.median(frame.loc[(frame.embed == "ldvae") & (frame.dataset == "v3"), "median_rank"])),
        "frac_fa_at_ceiling": float((fa >= 28).mean()),
    }
    print(json.dumps(gated, sort_keys=True))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embed", choices=("pca", "fa", "ldvae"))
    p.add_argument("--dataset", choices=("v2", "v3"))
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--summarise", action="store_true")
    a = p.parse_args()
    if a.summarise:
        if a.embed or a.dataset:
            p.error("--summarise cannot be combined with --embed or --dataset")
        summarise()
    elif a.embed and a.dataset:
        fit(a.embed, a.dataset, a.seed)
    else:
        p.error("--embed and --dataset are required unless --summarise is used")


if __name__ == "__main__":
    main()
