#!/usr/bin/env python
import argparse
import hashlib
import json
import os
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'T02b'
K2 = ROOT / 'reports' / 'K2'
T01D = ROOT / 'reports' / 'T01d'
DATA = ROOT / 'data' / 'pbmc_v2_v3'
CELLTYPIST_HOME = DATA / 'celltypist_home'
os.environ['CELLTYPIST_FOLDER'] = str(CELLTYPIST_HOME)

import celltypist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.stats import ortho_group
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig

LAT = {'a': K2 / 'latent_pca30_a.npy', 'b': K2 / 'latent_pca30_b.npy'}
OPS = {'a': K2 / 'ops_a.npz', 'b': K2 / 'ops_b.npz'}
V2 = DATA / 'pbmc_v2.h5ad'
V3 = DATA / 'pbmc_v3.h5ad'
V3_THIN = DATA / 'pbmc_v3_thin.h5ad'
CELLS = DATA / 'k2_cells.csv'
FINE = DATA / 'fine_labels.csv'
GENES = T01D / 'genes.txt'
INPUTS = (LAT['a'], LAT['b'], OPS['a'], OPS['b'], CELLS, V2, V3, V3_THIN, GENES)
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
D, R = 30, 10
DEVIATIONS = [{'what':'Repeated the ten prescribed deterministic drift fits after correcting the group-label array serialization used by group_means.npz.', 'why':'The first fit-stage attempt wrote groups with object dtype, which prevented allow_pickle=False loading; the corrected rerun used the same inputs, configuration, seeds and outputs.'}]
KNOWN = [
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility',
]

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def gates_path(): return OUT / 'gates.json'
def mat_path(dataset, seed): return OUT / f'M_{dataset}_{seed}.npy'
def drift_path(dataset, seed): return OUT / f'drift_{dataset}_{seed}.pt'
def group_means_path(): return OUT / 'group_means.npz'
def loadings_path(): return OUT / 'loadings.npz'
def group_table_path(): return OUT / 'group_table.csv'
def get_gates(): return json.loads(gates_path().read_text())
def put_gates(g): gates_path().write_text(json.dumps(g, sort_keys=True) + '\n')

def load_ops(dataset):
    with np.load(OPS[dataset], allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def load_matched_data():
    cells = pd.read_csv(CELLS)
    a0 = sc.read_h5ad(V2)
    a0 = a0[a0.obs['qc_pass'].to_numpy(dtype=bool)].copy()
    b0 = sc.read_h5ad(V3_THIN)
    bar_a = cells.loc[cells.dataset == 'A', 'barcode'].astype(str).to_numpy()
    bar_b = cells.loc[cells.dataset == 'B', 'barcode'].astype(str).to_numpy()
    a, b = a0[bar_a].copy(), b0[bar_b].copy()
    for dataset, ad, bars in (('a', a, bar_a), ('b', b, bar_b)):
        op = load_ops(dataset)
        if not np.array_equal(ad.obs_names.astype(str).to_numpy(), bars):
            raise RuntimeError(f'k2_cells ordering mismatch for {dataset}')
        if not np.array_equal(ad.obs_names.astype(str).to_numpy(), op['barcode'].astype(str)):
            raise RuntimeError(f'K2 barcode mismatch for {dataset}')
        if not np.array_equal(ad.obs['cell_type'].astype(str).to_numpy(), op['cell_type'].astype(str)):
            raise RuntimeError(f'K2 cell-type mismatch for {dataset}')
    return a, b

def normalize_for_celltypist(ad):
    ad = ad.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad

def label_original(path):
    ad = normalize_for_celltypist(sc.read_h5ad(path))
    result = celltypist.annotate(ad, model='Immune_All_Low.pkl', majority_voting=False)
    labels = result.predicted_labels['predicted_labels'].astype(str)
    labels.index = labels.index.astype(str)
    return labels

def labels_stage():
    CELLTYPIST_HOME.mkdir(parents=True, exist_ok=True)
    celltypist.models.download_models(model='Immune_All_Low.pkl')
    full_a, full_b = label_original(V2), label_original(V3)
    cells = pd.read_csv(CELLS)
    rows = []
    for dataset, labels, label_name in (('A', full_a, 'a'), ('B', full_b, 'b')):
        bars = cells.loc[cells.dataset == dataset, 'barcode'].astype(str).to_numpy()
        mapped = labels.reindex(bars)
        if mapped.isna().any():
            raise RuntimeError(f'missing CellTypist labels for {label_name}: {int(mapped.isna().sum())}')
        rows.extend({'barcode': barcode, 'dataset': dataset, 'fine_label': fine} for barcode, fine in zip(bars, mapped.to_numpy()))
    fine = pd.DataFrame(rows, columns=['barcode', 'dataset', 'fine_label'])
    fine.to_csv(FINE, index=False)
    counts_a = fine.loc[fine.dataset == 'A', 'fine_label'].value_counts()
    counts_b = fine.loc[fine.dataset == 'B', 'fine_label'].value_counts()
    groups = sorted(set(counts_a.index.astype(str)) | set(counts_b.index.astype(str)))
    table = pd.DataFrame([{
        'fine_label': group, 'n_a': int(counts_a.get(group, 0)), 'n_b': int(counts_b.get(group, 0)),
        'used': bool(int(counts_a.get(group, 0)) >= 40 and int(counts_b.get(group, 0)) >= 40)
    } for group in groups])
    table.to_csv(group_table_path(), index=False)
    print(json.dumps({'n_fine_labels': int(len(groups)), 'n_used_groups': int(table.used.sum())}, sort_keys=True))

def norm_pca(ad, genes):
    missing = [gene for gene in genes if gene not in ad.var_names]
    if missing:
        raise RuntimeError(f'missing requested genes: {len(missing)}')
    xad = ad[:, genes].copy()
    sc.pp.normalize_total(xad, target_sum=1e4)
    sc.pp.log1p(xad)
    x = xad.X.toarray() if hasattr(xad.X, 'toarray') else np.asarray(xad.X)
    pca = PCA(n_components=D, svd_solver='full', random_state=20261014)
    return pca.fit_transform(np.asarray(x, dtype=np.float32)).astype(np.float32), pca.components_.astype(np.float64)

def tangent(z):
    nbrs = NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    bases = np.empty((len(z), D, R), dtype=np.float32)
    for i in range(len(z)):
        _, u = np.linalg.eigh(np.cov(z[nbrs[i]], rowvar=False))
        bases[i] = u[:, -R:]
    return bases

def symmetrize(j): return (j + np.swapaxes(j, 1, 2)) / 2.0

def drift_config(seed):
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=seed, batch_size=256,
                             loss='esm', lr=5e-4, beta=1.0, esm_laplacian='hutchinson', device='cpu')

def fine_labels_in_k2_order():
    fine = pd.read_csv(FINE)
    result = {}
    for dataset, letter in (('A', 'a'), ('B', 'b')):
        op = load_ops(letter)
        labels = fine.loc[fine.dataset == dataset].set_index('barcode')['fine_label']
        ordered = labels.reindex(op['barcode'].astype(str))
        if ordered.isna().any():
            raise RuntimeError(f'fine label ordering mismatch for {letter}')
        result[letter] = ordered.to_numpy(dtype=str)
    return result

def coarse_tested_types():
    oa, ob = load_ops('a'), load_ops('b')
    la, lb = oa['cell_type'].astype(str), ob['cell_type'].astype(str)
    return sorted(t for t in set(la) & set(lb) if np.sum(la == t) >= 50 and np.sum(lb == t) >= 50)

def fits_stage():
    a, b = load_matched_data()
    genes = [x.strip() for x in GENES.read_text().splitlines() if x.strip()]
    if len(genes) != 2000 or genes != sorted(genes):
        raise RuntimeError('T01d genes.txt is not the specified sorted 2000-gene list')
    ads = {'a': a, 'b': b}
    symmetry_rows, repro = [], []
    pcas = {}
    for dataset in ('a', 'b'):
        z_re, w = norm_pca(ads[dataset], genes)
        z = np.load(LAT[dataset]).astype(np.float32)
        if z.shape != z_re.shape:
            raise RuntimeError(f'K2 latent shape mismatch for {dataset}')
        pcas[dataset] = (z_re, w)
        basis = tangent(z)
        op = load_ops(dataset)
        for seed in SEEDS:
            drift = EnergyDriftAdapter(D, config=drift_config(seed)).fit(z, verbose=True)
            check = drift.symmetry_check(z, n_cells=100)
            torch.save({'state_dict': drift._model.state_dict(), 'config': asdict(drift.config)}, drift_path(dataset, seed))
            m = symmetrize(drift.jacobian_autodiff(z).astype(np.float64))
            np.save(mat_path(dataset, seed), m)
            symmetry_rows.append({'dataset': dataset, 'seed': seed, 'median_asym': float(check['median_asym']), 'max_asym': float(check['max_asym'])})
            if seed == 20261014:
                jr = symmetrize(np.einsum('ndi,ndm,nmj->nij', basis.astype(np.float64), m, basis.astype(np.float64)))
                lam = np.linalg.eigvalsh(jr).astype(np.float64)
                repro.append(float(np.max(np.abs(lam - op['lam_seeds'][0].astype(np.float64)))))
    pd.DataFrame(symmetry_rows).to_csv(OUT / 'fit_symmetry.csv', index=False)
    np.savez_compressed(loadings_path(), W_A=pcas['a'][1], W_B=pcas['b'][1])
    labels = fine_labels_in_k2_order()
    table = pd.read_csv(group_table_path())
    groups = table.loc[table.used, 'fine_label'].astype(str).to_numpy(dtype=str)
    weights_unscaled = np.asarray([min(np.sum(labels['a'] == group), np.sum(labels['b'] == group)) for group in groups], dtype=np.float64)
    weights = weights_unscaled / weights_unscaled.sum()
    mean_a, mean_b = [], []
    for group in groups:
        mean_a.append(np.asarray(np.load(mat_path('a', 20261014), mmap_mode='r')[labels['a'] == group].mean(axis=0), dtype=np.float64))
        mean_b.append(np.asarray(np.load(mat_path('b', 20261014), mmap_mode='r')[labels['b'] == group].mean(axis=0), dtype=np.float64))
    rng = np.random.default_rng(20261022)
    null_h1, null_h2 = [], []
    for dataset in ('a', 'b'):
        for seed_b in SEEDS[1:]:
            h1, h2 = [], []
            for group in groups:
                idx = np.flatnonzero(labels[dataset] == group)
                perm = rng.permutation(idx)
                h = len(perm) // 2
                h1.append(np.asarray(np.load(mat_path(dataset, 20261014), mmap_mode='r')[perm[:h]].mean(axis=0), dtype=np.float64))
                h2.append(np.asarray(np.load(mat_path(dataset, seed_b), mmap_mode='r')[perm[h:2*h]].mean(axis=0), dtype=np.float64))
            null_h1.append(np.stack(h1))
            null_h2.append(np.stack(h2))
    np.savez_compressed(group_means_path(), M_bar_a=np.stack(mean_a), M_bar_b=np.stack(mean_b), null_h1=np.stack(null_h1), null_h2=np.stack(null_h2), groups=groups, w=weights)
    g = {
        'n_cells_a': int(len(a)), 'n_cells_b': int(len(b)),
        'latent_maxdiff': float(max(np.max(np.abs(pcas['a'][0] - np.load(LAT['a']))), np.max(np.abs(pcas['b'][0] - np.load(LAT['b']))))),
        'k2_repro_max_abs': float(max(repro)),
        'symmetry_max': float(max(row['median_asym'] for row in symmetry_rows)),
    }
    put_gates(g)
    print(json.dumps(g, sort_keys=True))

def load_group_means():
    with np.load(group_means_path(), allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def load_loadings():
    with np.load(loadings_path(), allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def starts(q_load, null=False):
    yield 0, np.eye(D, dtype=np.float64)
    yield 1, np.eye(D, dtype=np.float64) if null else q_load
    for i in range(1, 20):
        yield i + 1, np.asarray(ortho_group.rvs(dim=D, random_state=20261021 + i), dtype=np.float64)

def loss_np(q, ma, mb, weights):
    w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()
    dif = np.einsum('ij,tjk,lk->til', q, ma, q) - mb
    return float(np.sum(w * np.sum(dif * dif, axis=(1, 2))) / np.sum(w * np.sum(mb * mb, axis=(1, 2))))

def estimate_q(ma, mb, weights, q_load, null=False, keep_starts=False):
    ma_t = torch.as_tensor(np.asarray(ma, dtype=np.float64), dtype=torch.float64, device='cpu')
    mb_t = torch.as_tensor(np.asarray(mb, dtype=np.float64), dtype=torch.float64, device='cpu')
    weight_np = np.asarray(weights, dtype=np.float64); weight_np = weight_np / weight_np.sum()
    wt = torch.as_tensor(weight_np, dtype=torch.float64, device='cpu')
    denominator = torch.sum(wt * torch.sum(mb_t * mb_t, dim=(1, 2)))
    output = []
    for start, q0 in starts(q_load, null=null):
        q0t = torch.as_tensor(q0, dtype=torch.float64, device='cpu')
        s = torch.zeros((D, D), dtype=torch.float64, device='cpu', requires_grad=True)
        optimizer = torch.optim.LBFGS([s], lr=1, max_iter=500, tolerance_grad=1e-10, tolerance_change=1e-12, line_search_fn='strong_wolfe')
        def closure():
            optimizer.zero_grad()
            q = q0t @ torch.matrix_exp(s - s.T)
            dif = q @ ma_t @ q.T - mb_t
            loss = torch.sum(wt * torch.sum(dif * dif, dim=(1, 2))) / denominator
            loss.backward()
            return loss
        optimizer.step(closure)
        with torch.no_grad():
            q = q0t @ torch.matrix_exp(s - s.T)
            dif = q @ ma_t @ q.T - mb_t
            final_loss = float((torch.sum(wt * torch.sum(dif * dif, dim=(1, 2))) / denominator).cpu().item())
        output.append((start, final_loss, q.detach().cpu().numpy().astype(np.float64)))
    best = min(output, key=lambda row: row[1])
    if keep_starts:
        return best[2], best[1], best[0], [(row[0], row[1]) for row in output]
    return best[2], best[1]

def quotient(a, b):
    if a == 0 and b == 0: return 1.0
    return float('inf') if b == 0 else float(a / b)

def align_stage():
    tm = load_group_means()
    ma, mb, groups, weights = tm['M_bar_a'], tm['M_bar_b'], tm['groups'].astype(str), tm['w'].astype(np.float64)
    g = get_gates()
    g.update({'n_groups': int(len(groups))})
    invalid = g['n_cells_a'] != g['n_cells_b'] or g['latent_maxdiff'] > 1e-4 or g['k2_repro_max_abs'] > 1e-5 or g['symmetry_max'] > 1e-4
    if invalid:
        g['outcome'] = 'inconclusive'; g['validity_revise'] = True; put_gates(g); print(json.dumps({k:v for k,v in g.items() if k!='validity_revise'}, sort_keys=True)); return
    if g['n_groups'] < 6:
        g['outcome'] = 'stopped'; put_gates(g); print(json.dumps(g, sort_keys=True)); return
    ld = load_loadings()
    wa, wb = ld['W_A'], ld['W_B']
    u_load, singular, vt_load = np.linalg.svd(wb @ wa.T, full_matrices=False)
    q_load = (u_load @ vt_load).astype(np.float64)
    g['subspace_min_sv'] = float(np.min(singular))
    q, loss, best_start, start_rows = estimate_q(ma, mb, weights, q_load, keep_starts=True)
    np.save(OUT / 'Q.npy', q); np.save(OUT / 'Q_load.npy', q_load)
    starts_df = pd.DataFrame(start_rows, columns=['start', 'final_loss']); starts_df.to_csv(OUT / 'starts.csv', index=False)
    null_rows = []
    for i, (dataset, seed_b) in enumerate((d, s) for d in ('a', 'b') for s in SEEDS[1:]):
        q_null, loss_null = estimate_q(tm['null_h1'][i], tm['null_h2'][i], weights, q_load, null=True)
        null_rows.append({'dataset': dataset, 'seed_b': int(seed_b), 'resid': float(np.sqrt(loss_null))})
    null_df = pd.DataFrame(null_rows); null_df.to_csv(OUT / 'null.csv', index=False)
    held_rows = []
    for i, group in enumerate(groups):
        keep = np.arange(len(groups)) != i
        q_hold, _ = estimate_q(ma[keep], mb[keep], weights[keep], q_load)
        resid = np.linalg.norm(q_hold @ ma[i] @ q_hold.T - mb[i], ord='fro') / np.linalg.norm(mb[i], ord='fro')
        held_rows.append({'group': str(group), 'resid': float(resid)})
    held_df = pd.DataFrame(held_rows); held_df.to_csv(OUT / 'heldout.csv', index=False)
    labels = fine_labels_in_k2_order()
    za, zb = np.load(LAT['a']).astype(np.float64), np.load(LAT['b']).astype(np.float64)
    ca = np.stack([za[labels['a'] == group].mean(axis=0) for group in groups])
    cb = np.stack([zb[labels['b'] == group].mean(axis=0) for group in groups])
    u_pos, _, vt_pos = np.linalg.svd(ca.T @ cb, full_matrices=False)
    q_pos = (u_pos @ vt_pos).T.astype(np.float64); np.save(OUT / 'Q_pos.npy', q_pos)
    resid_null = float(np.median(null_df.resid))
    g.update({
        'loss_best': float(loss), 'best_start': int(best_start),
        'n_starts_near_best': int(np.sum(starts_df.final_loss.to_numpy() <= 1.01 * loss)),
        'resid_identity': float(np.sqrt(loss_np(np.eye(D), ma, mb, weights))),
        'resid_load': float(np.sqrt(loss_np(q_load, ma, mb, weights))),
        'resid_aligned': float(np.sqrt(loss)), 'resid_null': resid_null,
        'frame_ratio': quotient(float(np.sqrt(loss)), resid_null),
        'heldout_resid': float(held_df.resid.mean()),
        'heldout_ratio': quotient(float(held_df.resid.mean()), resid_null),
        'resid_pos': float(np.sqrt(loss_np(q_pos, ma, mb, weights))),
        'q_op_vs_pos': float(np.linalg.norm(q - q_pos, ord='fro') / np.sqrt(2 * D)),
        'q_op_vs_load': float(np.linalg.norm(q - q_load, ord='fro') / np.sqrt(2 * D)),
    })
    if g['n_starts_near_best'] < 2:
        g['outcome'] = 'stopped'
    elif g['frame_ratio'] <= 1.5 and g['heldout_ratio'] <= 1.5:
        g['outcome'] = 'frame agreement'
    elif g['frame_ratio'] > 3.0 or g['heldout_ratio'] > 3.0:
        g['outcome'] = 'no frame agreement'
    else:
        g['outcome'] = 'inconclusive'
    put_gates(g)
    print(json.dumps({k:v for k,v in g.items() if k!='validity_revise'}, sort_keys=True))

def gene_matrix(ad, genes):
    xad = ad[:, genes].copy()
    sc.pp.normalize_total(xad, target_sum=1e4)
    sc.pp.log1p(xad)
    return xad.X.toarray() if hasattr(xad.X, 'toarray') else np.asarray(xad.X)

def make_umap(rep):
    ad = sc.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32))
    ad.obsm['view'] = np.asarray(rep)
    sc.pp.neighbors(ad, n_neighbors=15, use_rep='view', random_state=20261014)
    sc.tl.umap(ad, min_dist=0.5, random_state=20261014)
    return np.asarray(ad.obsm['X_umap'], dtype=np.float64)

def mixing(rep, dataset):
    nbrs = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(dataset[nbrs] != dataset[:, None], axis=1)) / 0.5)

def plot_umap(coords, dataset, cell_type, scores):
    fig, axes = plt.subplots(3, 2, figsize=(10, 13))
    names = ('joint', 'sep_noA', 'sep_A')
    fine_types = sorted(set(cell_type.astype(str)))
    colors = plt.get_cmap('tab10')
    for row, name in enumerate(names):
        xy = coords[name]; mix, sil = scores[name]
        ax = axes[row, 0]
        for ds, color in (('A', 'C0'), ('B', 'C1')):
            mask = dataset == ds
            ax.scatter(xy[mask,0], xy[mask,1], s=4, alpha=.7, color=color, label=ds, linewidths=0)
        ax.set_title(f'{name}: dataset; mix={mix:.3f}, sil={sil:.3f}'); ax.legend(markerscale=2, fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        ax = axes[row, 1]
        for i, label in enumerate(fine_types):
            mask = cell_type == label
            ax.scatter(xy[mask,0], xy[mask,1], s=4, alpha=.7, color=colors(i % 10), label=label, linewidths=0)
        ax.set_title(f'{name}: cell type; mix={mix:.3f}, sil={sil:.3f}'); ax.legend(markerscale=2, fontsize=7, ncol=2)
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(OUT / 'fig_umap.png', dpi=180); plt.close(fig)

def plot_frames(ma, mb, q, groups):
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(1, 11); cmap = plt.get_cmap('tab20')
    for i, group in enumerate(groups):
        _, ua = np.linalg.eigh(q @ ma[i] @ q.T)
        _, ub = np.linalg.eigh(mb[i])
        _, ui = np.linalg.eigh(ma[i])
        aligned = np.abs(np.sum(ua[:, -10:] * ub[:, -10:], axis=0))
        identity = np.abs(np.sum(ui[:, -10:] * ub[:, -10:], axis=0))
        color = cmap(i % 20)
        ax.plot(x, aligned, marker='o', markersize=3, color=color, label=str(group))
        ax.plot(x, identity, linestyle='--', color=color, alpha=.6, label=f'{group} (Q=I)')
    ax.set_xlabel('Eigenvector rank among largest 10 (ascending)'); ax.set_ylabel('Absolute cosine with B eigenvector'); ax.set_ylim(-.02, 1.02)
    ax.legend(fontsize=5.5, ncol=2); fig.tight_layout(); fig.savefig(OUT / 'fig_frames.png', dpi=180); plt.close(fig)

def views_stage():
    g = get_gates()
    if g['outcome'] == 'stopped':
        print('Stopped by pre-registered alignment rule; views not run.'); return
    a, b = load_matched_data(); q = np.load(OUT / 'Q.npy').astype(np.float64)
    genes = [x.strip() for x in GENES.read_text().splitlines() if x.strip()]
    joint = PCA(n_components=D, svd_solver='full', random_state=20261014).fit_transform(np.vstack([gene_matrix(a, genes), gene_matrix(b, genes)]))
    za, zb = np.load(LAT['a']).astype(np.float64), np.load(LAT['b']).astype(np.float64)
    reps = {'joint': joint, 'sep_noA': np.vstack([za, zb]), 'sep_A': np.vstack([za @ q.T, zb])}
    dataset = np.asarray(['A'] * len(a) + ['B'] * len(b), dtype=str)
    coarse = np.concatenate([a.obs['cell_type'].astype(str).to_numpy(), b.obs['cell_type'].astype(str).to_numpy()]).astype(str)
    tested = coarse_tested_types(); mask = np.isin(coarse, tested)
    coords, scores = {}, {}
    for name, rep in reps.items():
        coords[name] = make_umap(rep)
        scores[name] = (mixing(rep, dataset), float(silhouette_score(rep[mask], coarse[mask], metric='euclidean')))
    np.savez_compressed(OUT / 'umap_views.npz', joint=coords['joint'], sep_noA=coords['sep_noA'], sep_A=coords['sep_A'], dataset=dataset, cell_type=coarse)
    plot_umap(coords, dataset, coarse, scores)
    tm = load_group_means(); plot_frames(tm['M_bar_a'], tm['M_bar_b'], q, tm['groups'].astype(str))
    g.update({'mix_joint': scores['joint'][0], 'mix_sep_noA': scores['sep_noA'][0], 'mix_sep_A': scores['sep_A'][0], 'sil_joint': scores['joint'][1], 'sil_sep_noA': scores['sep_noA'][1], 'sil_sep_A': scores['sep_A'][1]})
    put_gates(g)
    print(json.dumps({key:g[key] for key in ('mix_joint','mix_sep_noA','mix_sep_A','sil_joint','sil_sep_noA','sil_sep_A')}, sort_keys=True))

def report_stage():
    g = get_gates(); stopped = g.get('outcome') == 'stopped'
    files = ['requirements.lock', 'group_table.csv', 'group_means.npz', 'loadings.npz', 'scripts/t02b_run.py']
    if (OUT / 'Q.npy').exists(): files += ['Q.npy', 'Q_load.npy', 'Q_pos.npy', 'starts.csv', 'null.csv', 'heldout.csv']
    if not stopped: files += ['umap_views.npz', 'fig_umap.png', 'fig_frames.png']
    command = 'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T02b/scripts/t02b_run.py'
    fit_keys = ('n_cells_a','n_cells_b','latent_maxdiff','k2_repro_max_abs','symmetry_max')
    align_keys = ('n_groups','subspace_min_sv','loss_best','best_start','n_starts_near_best','resid_identity','resid_load','resid_aligned','resid_null','frame_ratio','heldout_resid','heldout_ratio','resid_pos','q_op_vs_pos','q_op_vs_load','outcome')
    regen = {path: command + ' --stage all' for path in files}
    regen.update({key: command + ' --stage fits' for key in fit_keys if key in g})
    regen.update({key: command + ' --stage align' for key in align_keys if key in g})
    regen.update({key: command + ' --stage views' for key in ('mix_joint','mix_sep_noA','mix_sep_A','sil_joint','sil_sep_noA','sil_sep_A') if key in g})
    report = {
        'task_id':'T02b', 'executor':'manus', 'branch':'task/T02b',
        'commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'], text=True).strip(),
        'seed':20261014, 'started':'2026-10-15T15:15:00-06:00', 'finished':datetime.now(timezone.utc).isoformat(),
        'status':'stopped' if stopped else 'complete', 'gated':{key:value for key,value in g.items() if key != 'validity_revise'},
        'deliverables':[{'path':path,'sha256':sha(OUT/path)} for path in files], 'regeneration':regen,
        'deviations':DEVIATIONS, 'deviations_requested':[],
        'observations':[f'Input SHA-256: {path.relative_to(ROOT)} = {sha(path)}.' for path in INPUTS],
        'known_failures':KNOWN, 'lockfile_sha256':sha(OUT/'requirements.lock'),
    }
    if stopped:
        if g['n_groups'] < 6: report['stop_rule'] = 'n_groups < 6 → inconclusive (stop) (not enough correspondences).'
        else: report['stop_rule'] = 'n_starts_near_best < 2 → inconclusive (stop).'
    (OUT/'report.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report['gated'], sort_keys=True))

def main():
    p = argparse.ArgumentParser(); p.add_argument('--stage',required=True,choices=('labels','fits','align','views','report','all')); stage=p.parse_args().stage
    if stage in ('labels','all'): labels_stage()
    if stage in ('fits','all'): fits_stage()
    if stage in ('align','all'): align_stage()
    if stage in ('views','all'): views_stage()
    if stage in ('report','all'): report_stage()
if __name__ == '__main__': main()
