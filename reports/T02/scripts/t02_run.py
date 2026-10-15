#!/usr/bin/env python
import argparse
import hashlib
import json
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

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

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'T02'
K2 = ROOT / 'reports' / 'K2'
T01D = ROOT / 'reports' / 'T01d'
DATA = ROOT / 'data' / 'pbmc_v2_v3'
LAT = {'a': K2 / 'latent_pca30_a.npy', 'b': K2 / 'latent_pca30_b.npy'}
OPS = {'a': K2 / 'ops_a.npz', 'b': K2 / 'ops_b.npz'}
V2 = DATA / 'pbmc_v2.h5ad'
V3 = DATA / 'pbmc_v3_thin.h5ad'
CELLS = DATA / 'k2_cells.csv'
GENES = T01D / 'genes.txt'
INPUTS = (LAT['a'], LAT['b'], OPS['a'], OPS['b'], CELLS, V2, V3, GENES)
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
D = 30
R = 10
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
def basis_path(dataset): return OUT / f'basis_{dataset}.npy'
def rank_path(dataset): return OUT / f'rank_{dataset}.npy'
def fit_sym_path(): return OUT / 'fit_symmetry.csv'
def types_path(): return OUT / 'type_means.npz'

def get_gates():
    return json.loads(gates_path().read_text())

def put_gates(g):
    gates_path().write_text(json.dumps(g, sort_keys=True) + '\n')

def config(seed):
    return EnergyDriftConfig(
        hidden=96, depth=3, n_epochs=1500, seed=seed, batch_size=256,
        loss='esm', lr=5e-4, beta=1.0, esm_laplacian='hutchinson', device='cpu'
    )

def symmetrize(j):
    return (j + np.swapaxes(j, 1, 2)) / 2.0

def tangent(z):
    nn = NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z, return_distance=False)[:, 1:]
    ranks = np.empty(len(z), dtype=np.int32)
    basis = np.empty((len(z), D, R), dtype=np.float32)
    for i in range(len(z)):
        w, u = np.linalg.eigh(np.cov(z[nn[i]], rowvar=False))
        ranks[i] = int(np.sum(w > 0.05 * w.max()))
        basis[i] = u[:, -R:]
    return ranks, basis

def load_k2_ops(dataset):
    with np.load(OPS[dataset], allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def load_matched_data():
    cells = pd.read_csv(CELLS)
    a0 = sc.read_h5ad(V2)
    a0 = a0[a0.obs['qc_pass'].to_numpy(dtype=bool)].copy()
    b0 = sc.read_h5ad(V3)
    bar_a = cells.loc[cells['dataset'] == 'A', 'barcode'].astype(str).to_numpy()
    bar_b = cells.loc[cells['dataset'] == 'B', 'barcode'].astype(str).to_numpy()
    a = a0[bar_a].copy()
    b = b0[bar_b].copy()
    for dataset, ad, bars in (('a', a, bar_a), ('b', b, bar_b)):
        o = load_k2_ops(dataset)
        if not np.array_equal(ad.obs_names.astype(str).to_numpy(), bars):
            raise RuntimeError(f'k2_cells ordering mismatch for {dataset}')
        if not np.array_equal(ad.obs_names.astype(str).to_numpy(), o['barcode'].astype(str)):
            raise RuntimeError(f'K2 barcode mismatch for {dataset}')
        if not np.array_equal(ad.obs['cell_type'].astype(str).to_numpy(), o['cell_type'].astype(str)):
            raise RuntimeError(f'K2 cell-type mismatch for {dataset}')
    return a, b

def tested_types(oa, ob):
    la = oa['cell_type'].astype(str)
    lb = ob['cell_type'].astype(str)
    return sorted(
        t for t in set(la) & set(lb)
        if int(np.sum(la == t)) >= 50 and int(np.sum(lb == t)) >= 50
    )

def fit_stage():
    ads = dict(zip(('a', 'b'), load_matched_data()))
    symmetry_rows = []
    repro = []
    for dataset in ('a', 'b'):
        z = np.load(LAT[dataset]).astype(np.float32)
        ad = ads[dataset]
        op = load_k2_ops(dataset)
        if z.shape != (len(ad), D):
            raise RuntimeError(f'K2 latent shape mismatch for {dataset}: {z.shape}')
        ranks, basis = tangent(z)
        np.save(rank_path(dataset), ranks)
        np.save(basis_path(dataset), basis)
        for seed in SEEDS:
            drift = EnergyDriftAdapter(D, config=config(seed)).fit(z, verbose=True)
            check = drift.symmetry_check(z, n_cells=100)
            torch.save({'state_dict': drift._model.state_dict(), 'config': asdict(drift.config)}, drift_path(dataset, seed))
            m = symmetrize(drift.jacobian_autodiff(z).astype(np.float64))
            np.save(mat_path(dataset, seed), m)
            symmetry_rows.append({
                'dataset': dataset, 'seed': seed,
                'median_asym': float(check['median_asym']), 'max_asym': float(check['max_asym'])
            })
            if seed == 20261014:
                jr = symmetrize(np.einsum('ndi,ndm,nmj->nij', basis.astype(np.float64), m, basis.astype(np.float64)))
                lam = np.linalg.eigvalsh(jr).astype(np.float64)
                repro.append(float(np.max(np.abs(lam - op['lam_seeds'][0].astype(np.float64)))))
    pd.DataFrame(symmetry_rows).to_csv(fit_sym_path(), index=False)
    oa, ob = load_k2_ops('a'), load_k2_ops('b')
    types = tested_types(oa, ob)
    mba, mbb, ns = [], [], []
    for t in types:
        ma = np.load(mat_path('a', 20261014), mmap_mode='r')[oa['cell_type'].astype(str) == t]
        mb = np.load(mat_path('b', 20261014), mmap_mode='r')[ob['cell_type'].astype(str) == t]
        if len(ma) != len(mb):
            raise RuntimeError(f'K2 count-match mismatch for {t}')
        mba.append(np.asarray(ma.mean(axis=0), dtype=np.float64))
        mbb.append(np.asarray(mb.mean(axis=0), dtype=np.float64))
        ns.append(len(ma))
    np.savez_compressed(types_path(), M_bar_a=np.stack(mba), M_bar_b=np.stack(mbb), types=np.asarray(types, dtype=str), n=np.asarray(ns, dtype=np.int64))
    g = {
        'n_cells_a': int(len(ads['a'])), 'n_cells_b': int(len(ads['b'])),
        'k2_repro_max_abs': float(max(repro)),
        'symmetry_max': float(max(row['median_asym'] for row in symmetry_rows)),
    }
    put_gates(g)
    print(json.dumps(g, sort_keys=True))

def load_type_means():
    with np.load(types_path(), allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def q_starts():
    yield 0, np.eye(D, dtype=np.float64)
    for i in range(1, 20):
        yield i, np.asarray(ortho_group.rvs(dim=D, random_state=20261021 + i), dtype=np.float64)

def loss_np(q, ma, mb, n):
    w = np.asarray(n, dtype=np.float64) / np.sum(n)
    dif = np.einsum('ij,tjk,lk->til', q, ma, q) - mb
    numerator = np.sum(w * np.sum(dif * dif, axis=(1, 2)))
    denominator = np.sum(w * np.sum(mb * mb, axis=(1, 2)))
    return float(numerator / denominator)

def estimate_q(ma, mb, n, retain_starts=False):
    ma_t = torch.as_tensor(np.asarray(ma, dtype=np.float64), dtype=torch.float64, device='cpu')
    mb_t = torch.as_tensor(np.asarray(mb, dtype=np.float64), dtype=torch.float64, device='cpu')
    w_t = torch.as_tensor(np.asarray(n, dtype=np.float64) / np.sum(n), dtype=torch.float64, device='cpu')
    denominator = torch.sum(w_t * torch.sum(mb_t * mb_t, dim=(1, 2)))
    results = []
    for start, q0 in q_starts():
        q0_t = torch.as_tensor(q0, dtype=torch.float64, device='cpu')
        s = torch.zeros((D, D), dtype=torch.float64, device='cpu', requires_grad=True)
        optimizer = torch.optim.LBFGS(
            [s], lr=1, max_iter=500, tolerance_grad=1e-10,
            tolerance_change=1e-12, line_search_fn='strong_wolfe'
        )
        def closure():
            optimizer.zero_grad()
            q = q0_t @ torch.matrix_exp(s - s.T)
            dif = q @ ma_t @ q.T - mb_t
            loss = torch.sum(w_t * torch.sum(dif * dif, dim=(1, 2))) / denominator
            loss.backward()
            return loss
        optimizer.step(closure)
        with torch.no_grad():
            q = q0_t @ torch.matrix_exp(s - s.T)
            dif = q @ ma_t @ q.T - mb_t
            final_loss = float((torch.sum(w_t * torch.sum(dif * dif, dim=(1, 2))) / denominator).cpu().item())
        results.append((start, final_loss, q.detach().cpu().numpy().astype(np.float64)))
    best = min(results, key=lambda x: x[1])
    if retain_starts:
        return best[2], best[1], [(x[0], x[1]) for x in results]
    return best[2], best[1]

def quotient(a, b):
    if a == 0 and b == 0: return 1.0
    return float('inf') if b == 0 else float(a / b)

def align_stage():
    tm = load_type_means()
    ma, mb, types, n = tm['M_bar_a'], tm['M_bar_b'], tm['types'].astype(str), tm['n'].astype(np.int64)
    oa, ob = load_k2_ops('a'), load_k2_ops('b')
    q, loss, starts = estimate_q(ma, mb, n, retain_starts=True)
    np.save(OUT / 'Q.npy', q)
    starts_df = pd.DataFrame(starts, columns=['start', 'final_loss'])
    starts_df.to_csv(OUT / 'starts.csv', index=False)
    null_rows = []
    rng = np.random.default_rng(20261022)
    mats = {dataset: {seed: np.load(mat_path(dataset, seed), mmap_mode='r') for seed in SEEDS} for dataset in ('a', 'b')}
    labels = {'a': oa['cell_type'].astype(str), 'b': ob['cell_type'].astype(str)}
    for dataset in ('a', 'b'):
        for seed_b in SEEDS[1:]:
            h1, h2, half_n = [], [], []
            for t in types:
                idx = np.flatnonzero(labels[dataset] == t)
                perm = rng.permutation(idx)
                h = len(perm) // 2
                h1.append(np.asarray(mats[dataset][20261014][perm[:h]].mean(axis=0), dtype=np.float64))
                h2.append(np.asarray(mats[dataset][seed_b][perm[h:2*h]].mean(axis=0), dtype=np.float64))
                half_n.append(h)
            q_null, _ = estimate_q(np.stack(h1), np.stack(h2), np.asarray(half_n))
            null_rows.append({'dataset': dataset, 'seed_b': int(seed_b), 'resid': float(np.sqrt(loss_np(q_null, np.stack(h1), np.stack(h2), np.asarray(half_n))))})
    null_df = pd.DataFrame(null_rows)
    null_df.to_csv(OUT / 'null.csv', index=False)
    heldout_rows = []
    for i, t in enumerate(types):
        keep = np.arange(len(types)) != i
        q_hold, _ = estimate_q(ma[keep], mb[keep], n[keep])
        numerator = np.linalg.norm(q_hold @ ma[i] @ q_hold.T - mb[i], ord='fro')
        heldout_rows.append({'type': str(t), 'resid': float(numerator / np.linalg.norm(mb[i], ord='fro'))})
    heldout_df = pd.DataFrame(heldout_rows)
    heldout_df.to_csv(OUT / 'heldout.csv', index=False)
    za, zb = np.load(LAT['a']).astype(np.float64), np.load(LAT['b']).astype(np.float64)
    ca = np.stack([za[labels['a'] == t].mean(axis=0) for t in types])
    cb = np.stack([zb[labels['b'] == t].mean(axis=0) for t in types])
    u, _, vt = np.linalg.svd(ca.T @ cb, full_matrices=False)
    q_pos = (u @ vt).T.astype(np.float64)
    np.save(OUT / 'Q_pos.npy', q_pos)
    g = get_gates()
    resid_null = float(np.median(null_df['resid']))
    g.update({
        'n_cell_types_tested': int(len(types)),
        'loss_best': float(loss),
        'n_starts_near_best': int(np.sum(starts_df['final_loss'].to_numpy() <= 1.01 * loss)),
        'resid_identity': float(np.sqrt(loss_np(np.eye(D), ma, mb, n))),
        'resid_aligned': float(np.sqrt(loss)),
        'resid_null': resid_null,
        'frame_ratio': quotient(float(np.sqrt(loss)), resid_null),
        'heldout_resid': float(heldout_df['resid'].mean()),
        'heldout_ratio': quotient(float(heldout_df['resid'].mean()), resid_null),
        'resid_pos': float(np.sqrt(loss_np(q_pos, ma, mb, n))),
        'q_op_vs_pos': float(np.linalg.norm(q - q_pos, ord='fro') / np.sqrt(2 * D)),
    })
    invalid = g['n_cells_a'] != g['n_cells_b'] or g['k2_repro_max_abs'] > 1e-5 or g['symmetry_max'] > 1e-4
    if invalid:
        g['outcome'] = 'inconclusive'
        g['validity_revise'] = True
    elif g['n_starts_near_best'] < 2:
        g['outcome'] = 'stopped'
    elif g['frame_ratio'] <= 1.5 and g['heldout_ratio'] <= 1.5:
        g['outcome'] = 'frame agreement'
    elif g['frame_ratio'] > 3.0 or g['heldout_ratio'] > 3.0:
        g['outcome'] = 'no frame agreement'
    else:
        g['outcome'] = 'inconclusive'
    put_gates(g)
    print(json.dumps({k: g[k] for k in g if k != 'validity_revise'}, sort_keys=True))

def gene_matrix(ad, genes):
    missing = [g for g in genes if g not in ad.var_names]
    if missing:
        raise RuntimeError(f'missing requested genes: {len(missing)}')
    xad = ad[:, genes].copy()
    sc.pp.normalize_total(xad, target_sum=1e4)
    sc.pp.log1p(xad)
    x = xad.X
    return x.toarray() if hasattr(x, 'toarray') else np.asarray(x)

def umap(rep):
    ad = sc.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32))
    ad.obsm['view'] = np.asarray(rep)
    sc.pp.neighbors(ad, n_neighbors=15, use_rep='view', random_state=20261014)
    sc.tl.umap(ad, min_dist=0.5, random_state=20261014)
    return np.asarray(ad.obsm['X_umap'], dtype=np.float64)

def mixing(rep, dataset):
    nn = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(dataset[nn] != dataset[:, None], axis=1)) / 0.5)

def plot_umap(views, coordinates, dataset, cell_type, types, scores):
    fig, axes = plt.subplots(3, 2, figsize=(10, 13))
    unique_types = sorted(set(cell_type.astype(str)))
    cmap = plt.get_cmap('tab10')
    type_color = {t: cmap(i % 10) for i, t in enumerate(unique_types)}
    for row, name in enumerate(('joint', 'sep_noA', 'sep_A')):
        xy = coordinates[name]
        mix, sil = scores[name]
        ax = axes[row, 0]
        for d, color in (('A', 'C0'), ('B', 'C1')):
            mask = dataset == d
            ax.scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=.7, color=color, label=d, linewidths=0)
        ax.set_title(f'{name}: dataset; mix={mix:.3f}, sil={sil:.3f}')
        ax.legend(markerscale=2, fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        ax = axes[row, 1]
        for t in unique_types:
            mask = cell_type == t
            ax.scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=.7, color=type_color[t], label=t, linewidths=0)
        ax.set_title(f'{name}: cell type; mix={mix:.3f}, sil={sil:.3f}')
        ax.legend(markerscale=2, fontsize=7, ncol=2)
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / 'fig_umap.png', dpi=180)
    plt.close(fig)

def plot_frames(ma, mb, q, types):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = np.arange(1, 11)
    colors = plt.get_cmap('tab10')
    for i, t in enumerate(types):
        _, ua = np.linalg.eigh(q @ ma[i] @ q.T)
        _, ub = np.linalg.eigh(mb[i])
        _, ui = np.linalg.eigh(ma[i])
        aligned = np.abs(np.sum(ua[:, -10:] * ub[:, -10:], axis=0))
        identity = np.abs(np.sum(ui[:, -10:] * ub[:, -10:], axis=0))
        color = colors(i % 10)
        ax.plot(x, aligned, marker='o', color=color, label=str(t))
        ax.plot(x, identity, linestyle='--', color=color, alpha=.65, label=f'{t} (Q=I)')
    ax.set_xlabel('Eigenvector rank among largest 10 (ascending)')
    ax.set_ylabel('Absolute cosine with B eigenvector')
    ax.set_ylim(-.02, 1.02)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / 'fig_frames.png', dpi=180)
    plt.close(fig)

def views_stage():
    g = get_gates()
    if g['outcome'] == 'stopped':
        print('Stopped by n_starts_near_best rule; views not run.')
        return
    a, b = load_matched_data()
    tm = load_type_means()
    types = tm['types'].astype(str)
    q = np.load(OUT / 'Q.npy').astype(np.float64)
    za, zb = np.load(LAT['a']).astype(np.float64), np.load(LAT['b']).astype(np.float64)
    genes = [x.strip() for x in GENES.read_text().splitlines() if x.strip()]
    if len(genes) != 2000:
        raise RuntimeError('T01d genes list is not 2000 genes')
    joint_input = np.vstack([gene_matrix(a, genes), gene_matrix(b, genes)])
    joint = PCA(n_components=30, svd_solver='full', random_state=20261014).fit_transform(joint_input)
    reps = {'joint': joint, 'sep_noA': np.vstack([za, zb]), 'sep_A': np.vstack([za @ q.T, zb])}
    dataset = np.asarray(['A'] * len(a) + ['B'] * len(b), dtype=str)
    cell_type = np.concatenate([a.obs['cell_type'].astype(str).to_numpy(), b.obs['cell_type'].astype(str).to_numpy()]).astype(str)
    coordinates, scores = {}, {}
    mask = np.isin(cell_type, types)
    for name, rep in reps.items():
        coordinates[name] = umap(rep)
        scores[name] = (mixing(rep, dataset), float(silhouette_score(rep[mask], cell_type[mask], metric='euclidean')))
    np.savez_compressed(OUT / 'umap_views.npz', joint=coordinates['joint'], sep_noA=coordinates['sep_noA'], sep_A=coordinates['sep_A'], dataset=dataset, cell_type=cell_type)
    plot_umap(reps, coordinates, dataset, cell_type, types, scores)
    plot_frames(tm['M_bar_a'], tm['M_bar_b'], q, types)
    g.update({
        'mix_joint': scores['joint'][0], 'mix_sep_noA': scores['sep_noA'][0], 'mix_sep_A': scores['sep_A'][0],
        'sil_joint': scores['joint'][1], 'sil_sep_noA': scores['sep_noA'][1], 'sil_sep_A': scores['sep_A'][1],
    })
    put_gates(g)
    print(json.dumps({k: g[k] for k in ('mix_joint','mix_sep_noA','mix_sep_A','sil_joint','sil_sep_noA','sil_sep_A')}, sort_keys=True))

def report_stage():
    g = get_gates()
    stopped = g.get('outcome') == 'stopped'
    files = ['requirements.lock', 'type_means.npz', 'Q.npy', 'Q_pos.npy', 'starts.csv', 'null.csv', 'heldout.csv', 'scripts/t02_run.py']
    if not stopped:
        files += ['umap_views.npz', 'fig_umap.png', 'fig_frames.png']
    command = 'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T02/scripts/t02_run.py'
    fit_keys = ('n_cells_a', 'n_cells_b', 'k2_repro_max_abs', 'symmetry_max')
    align_keys = ('n_cell_types_tested', 'loss_best', 'n_starts_near_best', 'resid_identity', 'resid_aligned', 'resid_null', 'frame_ratio', 'heldout_resid', 'heldout_ratio', 'resid_pos', 'q_op_vs_pos', 'outcome')
    regeneration = {x: command + ' --stage all' for x in files}
    regeneration.update({key: command + ' --stage fits' for key in fit_keys})
    regeneration.update({key: command + ' --stage align' for key in align_keys if key in g})
    regeneration.update({key: command + ' --stage views' for key in ('mix_joint','mix_sep_noA','mix_sep_A','sil_joint','sil_sep_noA','sil_sep_A') if key in g})
    obj = {
        'task_id': 'T02', 'executor': 'manus', 'branch': 'task/T02',
        'commit': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
        'seed': 20261014, 'started': '2026-10-15T15:02:00-06:00',
        'finished': datetime.now(timezone.utc).isoformat(),
        'status': 'stopped' if stopped else 'complete',
        'gated': {k: v for k, v in g.items() if k != 'validity_revise'},
        'deliverables': [{'path': x, 'sha256': sha(OUT / x)} for x in files],
        'regeneration': regeneration,
        'deviations': [], 'deviations_requested': [],
        'observations': [f'Input SHA-256: {p.relative_to(ROOT)} = {sha(p)}.' for p in INPUTS],
        'known_failures': KNOWN, 'lockfile_sha256': sha(OUT / 'requirements.lock'),
    }
    if stopped:
        obj['stop_rule'] = 'n_starts_near_best < 2 → inconclusive (stop) (the optimum is not found reproducibly).'
    (OUT / 'report.json').write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    print(json.dumps(obj['gated'], sort_keys=True))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', required=True, choices=('fits', 'align', 'views', 'report', 'all'))
    stage = parser.parse_args().stage
    if stage in ('fits', 'all'):
        fit_stage()
    if stage in ('align', 'all'):
        align_stage()
    if stage in ('views', 'all'):
        views_stage()
    if stage in ('report', 'all'):
        report_stage()

if __name__ == '__main__':
    main()
