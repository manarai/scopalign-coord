#!/usr/bin/env python
import argparse
import hashlib
import json
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import scipy.sparse as sp
import torch
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.compare import centroid_operator
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'T05a'
DATA = ROOT / 'data' / 'pbmc_v2_v3'
T01D = ROOT / 'reports' / 'T01d'
K2 = ROOT / 'reports' / 'K2'
T02B = ROOT / 'reports' / 'T02b'
T02C = ROOT / 'reports' / 'T02c'
CELLS = DATA / 'k2_cells.csv'
V2 = DATA / 'pbmc_v2.h5ad'
V3 = DATA / 'pbmc_v3_thin.h5ad'
GENES = T01D / 'genes.txt'
LAT_A = K2 / 'latent_pca30_a.npy'
LAT_B = K2 / 'latent_pca30_b.npy'
LOADINGS = T02B / 'loadings.npz'
T02C_REPORT = T02C / 'report.json'
INPUTS = (CELLS, V2, V3, GENES, LAT_A, LAT_B, LOADINGS, T02C_REPORT)
METHODS = ('joint', 'centre', 'harmony', 'scvi', 'ours')
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
D = 30
KNOWN = [
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility',
]

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for piece in iter(lambda: f.read(1 << 20), b''):
            h.update(piece)
    return h.hexdigest()

def gate_path(): return OUT / 'gates.json'
def get_gates(): return json.loads(gate_path().read_text())
def put_gates(gates): gate_path().write_text(json.dumps(gates, sort_keys=True) + '\n')
def emb_path(): return OUT / 'embeddings.npz'
def ops_path(): return OUT / 'centroid_ops.npz'
def drift_path(method, dataset, seed): return OUT / f'drift_{method}_{dataset}_{seed}.pt'
def matrix_path(method, dataset, seed): return OUT / f'M_{method}_{dataset}_{seed}.npy'

def load_matched():
    cells = pd.read_csv(CELLS)
    a0 = sc.read_h5ad(V2)
    a0 = a0[a0.obs['qc_pass'].to_numpy(dtype=bool)].copy()
    b0 = sc.read_h5ad(V3)
    bars_a = cells.loc[cells.dataset == 'A', 'barcode'].astype(str).to_numpy()
    bars_b = cells.loc[cells.dataset == 'B', 'barcode'].astype(str).to_numpy()
    a, b = a0[bars_a].copy(), b0[bars_b].copy()
    if not np.array_equal(a.obs_names.astype(str).to_numpy(), bars_a):
        raise RuntimeError('A barcode order does not equal k2_cells.csv')
    if not np.array_equal(b.obs_names.astype(str).to_numpy(), bars_b):
        raise RuntimeError('B barcode order does not equal k2_cells.csv')
    return a, b

def dense_log_matrix(adata, genes):
    x = adata[:, genes].copy()
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    return x.X.toarray() if sp.issparse(x.X) else np.asarray(x.X)

def raw_adata(a, b, genes):
    aa, bb = a[:, genes].copy(), b[:, genes].copy()
    if sp.issparse(aa.X) or sp.issparse(bb.X):
        x = sp.vstack([sp.csr_matrix(aa.X), sp.csr_matrix(bb.X)]).astype(np.float32)
    else:
        x = np.vstack([np.asarray(aa.X), np.asarray(bb.X)]).astype(np.float32)
    obs = pd.DataFrame({'dataset': ['A'] * aa.n_obs + ['B'] * bb.n_obs})
    var = pd.DataFrame(index=pd.Index(genes, dtype=str))
    return ad.AnnData(X=x, obs=obs, var=var)

def load_q():
    with np.load(LOADINGS, allow_pickle=False) as x:
        wa, wb = x['W_A'].astype(np.float64), x['W_B'].astype(np.float64)
    u, _, vt = np.linalg.svd(wb @ wa.T, full_matrices=False)
    return (u @ vt).astype(np.float64)

def embed_stage():
    a, b = load_matched()
    genes = [line.strip() for line in GENES.read_text().splitlines() if line.strip()]
    if len(genes) != 2000 or genes != sorted(genes):
        raise RuntimeError('T01d genes.txt must contain the sorted 2000-gene list')
    xa, xb = dense_log_matrix(a, genes), dense_log_matrix(b, genes)
    stacked = np.vstack([xa, xb])
    joint = PCA(n_components=D, svd_solver='full', random_state=20261014).fit_transform(stacked).astype(np.float64)
    centred = np.vstack([xa - xa.mean(axis=0), xb - xb.mean(axis=0)])
    centre = PCA(n_components=D, svd_solver='full', random_state=20261014).fit_transform(centred).astype(np.float64)
    dataset = np.asarray(['A'] * len(a) + ['B'] * len(b), dtype=str)
    meta = pd.DataFrame({'dataset': dataset})
    harmony = hm.run_harmony(joint, meta, 'dataset', random_state=20261014, max_iter_harmony=20).Z_corr.T.astype(np.float64)
    scvi.settings.seed = 20261014
    scvi_data = raw_adata(a, b, genes)
    scvi.model.SCVI.setup_anndata(scvi_data, batch_key='dataset')
    scvi_model = scvi.model.SCVI(scvi_data, n_latent=D)
    scvi_model.train(max_epochs=400, accelerator='cpu', devices=1)
    scvi_latent = scvi_model.get_latent_representation().astype(np.float64)
    za, zb = np.load(LAT_A).astype(np.float64), np.load(LAT_B).astype(np.float64)
    if za.shape != (len(a), D) or zb.shape != (len(b), D):
        raise RuntimeError('K2 latent array shape does not match k2_cells.csv')
    ours = np.vstack([za @ load_q().T, zb]).astype(np.float64)
    labels = np.concatenate([a.obs['cell_type'].astype(str).to_numpy(), b.obs['cell_type'].astype(str).to_numpy()]).astype(str)
    np.savez_compressed(emb_path(), joint=joint, centre=centre, harmony=harmony, scvi=scvi_latent, ours=ours, dataset=dataset, cell_type=labels)
    print(json.dumps({'n_cells_a': int(len(a)), 'n_cells_b': int(len(b)), 'n_genes': int(len(genes)), 'methods': list(METHODS)}, sort_keys=True))

def load_embeddings():
    with np.load(emb_path(), allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def tested_types(labels, dataset, threshold=50):
    a, b = labels[dataset == 'A'], labels[dataset == 'B']
    return sorted(label for label in set(a) & set(b) if np.sum(a == label) >= threshold and np.sum(b == label) >= threshold)

def mixing(embedding, dataset):
    nn = NearestNeighbors(n_neighbors=31).fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(dataset[nn] != dataset[:, None], axis=1)) / 0.5)

def score_values(embedding, dataset, labels, types):
    selected = np.isin(labels, types)
    return mixing(embedding, dataset), float(silhouette_score(embedding[selected], labels[selected], metric='euclidean'))

def scores_stage():
    e = load_embeddings()
    dataset, labels = e['dataset'].astype(str), e['cell_type'].astype(str)
    types = tested_types(labels, dataset)
    rows, gates = [], {}
    for method in METHODS:
        mix, sil = score_values(e[method].astype(np.float64), dataset, labels, types)
        rows.append({'method': method, 'mix': mix, 'sil': sil})
        gates[f'mix_{method}'] = mix
        gates[f'sil_{method}'] = sil
    score_df = pd.DataFrame(rows)
    score_df.to_csv(OUT / 'scores.csv', index=False)
    prior = json.loads(T02C_REPORT.read_text())['gated']
    ours = score_df.set_index('method').loc['ours']
    gates['ours_repro_max_abs'] = float(max(abs(float(ours['mix']) - float(prior['mix_sep_load'])), abs(float(ours['sil']) - float(prior['sil_sep_load']))))
    put_gates(gates)
    print(json.dumps(gates, sort_keys=True))

def drift_config(seed):
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=seed, batch_size=256, loss='esm', lr=5e-4, beta=1.0, esm_laplacian='hutchinson', device='cpu')

def sym(j): return (j + np.swapaxes(j, 1, 2)) / 2.0

def clusters_weights(labels, dataset):
    la, lb = labels[dataset == 'A'], labels[dataset == 'B']
    clusters = np.asarray(sorted(label for label in set(la) & set(lb) if np.sum(la == label) >= 100 and np.sum(lb == label) >= 100), dtype=str)
    mins = np.asarray([min(np.sum(la == label), np.sum(lb == label)) for label in clusters], dtype=np.float64)
    return clusters, mins / mins.sum()

def residual(a, b, w):
    weights = np.asarray(w, dtype=np.float64); weights /= weights.sum()
    delta = a - b
    return float(np.sqrt(np.sum(weights * np.sum(delta * delta, axis=(1, 2))) / np.sum(weights * np.sum(b * b, axis=(1, 2)))))

def outcome(ratio):
    if ratio <= 1.5: return 'agreement'
    if ratio > 3.0: return 'no agreement'
    return 'inconclusive'

def operators_stage():
    e = load_embeddings()
    dataset, labels = e['dataset'].astype(str), e['cell_type'].astype(str)
    clusters, weights = clusters_weights(labels, dataset)
    if len(clusters) < 3:
        gates = get_gates(); gates['n_clusters'] = int(len(clusters)); gates['status'] = 'stopped'; put_gates(gates); print(json.dumps(gates, sort_keys=True)); return
    all_arrays = {'clusters': clusters, 'w': weights}
    rows, sym_rows = [], []
    split = {'a': np.flatnonzero(dataset == 'A'), 'b': np.flatnonzero(dataset == 'B')}
    labels_by_dataset = {'a': labels[split['a']], 'b': labels[split['b']]}
    for method in METHODS:
        rng = np.random.default_rng(20261022)
        embedding = e[method].astype(np.float32)
        z = {'a': embedding[split['a']], 'b': embedding[split['b']]}
        matrices, centroids = {'a': {}, 'b': {}}, {'a': [], 'b': []}
        for dname in ('a', 'b'):
            for seed in SEEDS:
                drift = EnergyDriftAdapter(D, config=drift_config(seed)).fit(z[dname], verbose=True)
                check = drift.symmetry_check(z[dname], n_cells=100)
                torch.save({'state_dict': drift._model.state_dict(), 'config': asdict(drift.config)}, drift_path(method, dname, seed))
                m = sym(drift.jacobian_autodiff(z[dname]).astype(np.float64))
                np.save(matrix_path(method, dname, seed), m)
                matrices[dname][seed] = m
                sym_rows.append({'method': method, 'dataset': dname, 'seed': seed, 'median_asym': float(check['median_asym']), 'max_asym': float(check['max_asym'])})
                centroids[dname].append(np.stack([centroid_operator(z[dname][labels_by_dataset[dname] == cluster], m[labels_by_dataset[dname] == cluster], n_cells_centroid=30).astype(np.float64) for cluster in clusters]))
        null_h1, null_h2 = [], []
        for dname in ('a', 'b'):
            for seed_b in SEEDS[1:]:
                h1, h2 = [], []
                for cluster in clusters:
                    ix = np.flatnonzero(labels_by_dataset[dname] == cluster)
                    permutation = rng.permutation(ix)
                    half = len(permutation) // 2
                    h1.append(centroid_operator(z[dname][permutation[:half]], matrices[dname][20261014][permutation[:half]], n_cells_centroid=30).astype(np.float64))
                    h2.append(centroid_operator(z[dname][permutation[half:2*half]], matrices[dname][seed_b][permutation[half:2*half]], n_cells_centroid=30).astype(np.float64))
                null_h1.append(np.stack(h1)); null_h2.append(np.stack(h2))
        ca, cb = np.stack(centroids['a']), np.stack(centroids['b'])
        nh1, nh2 = np.stack(null_h1), np.stack(null_h2)
        frame_resid = residual(ca[0], cb[0], weights)
        null_values = [residual(nh1[i], nh2[i], weights) for i in range(8)]
        frame_null = float(np.median(null_values))
        ratio = float(frame_resid / frame_null) if frame_null else float('inf')
        rows.append({'method': method, 'frame_resid': frame_resid, 'frame_null': frame_null, 'frame_ratio': ratio, 'frame_outcome': outcome(ratio)})
        all_arrays[f'{method}_C_a'] = ca; all_arrays[f'{method}_C_b'] = cb; all_arrays[f'{method}_null_h1'] = nh1; all_arrays[f'{method}_null_h2'] = nh2
    np.savez_compressed(ops_path(), **all_arrays)
    pd.DataFrame(rows).to_csv(OUT / 'frames.csv', index=False)
    pd.DataFrame(sym_rows).to_csv(OUT / 'fit_symmetry.csv', index=False)
    gates = get_gates()
    gates['n_clusters'] = int(len(clusters))
    gates['symmetry_max'] = float(max(row['median_asym'] for row in sym_rows))
    for row in rows:
        gates[f"frame_ratio_{row['method']}"] = float(row['frame_ratio'])
        gates[f"frame_outcome_{row['method']}"] = row['frame_outcome']
    put_gates(gates)
    print(json.dumps({key: gates[key] for key in ['n_clusters', 'symmetry_max'] + [f'frame_ratio_{method}' for method in METHODS] + [f'frame_outcome_{method}' for method in METHODS]}, sort_keys=True))

def make_umap(embedding):
    a = ad.AnnData(X=np.zeros((len(embedding), 1), dtype=np.float32))
    a.obsm['view'] = np.asarray(embedding, dtype=np.float64)
    sc.pp.neighbors(a, n_neighbors=15, use_rep='view', random_state=20261014)
    sc.tl.umap(a, min_dist=0.5, random_state=20261014)
    return np.asarray(a.obsm['X_umap'], dtype=np.float64)

def views_stage():
    e = load_embeddings(); scores = pd.read_csv(OUT / 'scores.csv').set_index('method'); frames = pd.read_csv(OUT / 'frames.csv').set_index('method')
    dataset, labels = e['dataset'].astype(str), e['cell_type'].astype(str)
    fig, axes = plt.subplots(len(METHODS), 2, figsize=(10, 21)); colour_map = plt.get_cmap('tab10'); cell_types = sorted(set(labels))
    for r, method in enumerate(METHODS):
        xy = make_umap(e[method])
        mix, sil, ratio = float(scores.loc[method, 'mix']), float(scores.loc[method, 'sil']), float(frames.loc[method, 'frame_ratio'])
        ax = axes[r, 0]
        for ds, colour in (('A', 'C0'), ('B', 'C1')):
            mask = dataset == ds
            ax.scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=0.7, color=colour, label=ds, linewidths=0)
        ax.set_title(f'{method}: dataset; mix={mix:.3f}, sil={sil:.3f}, ratio={ratio:.3f}')
        ax.legend(markerscale=2, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
        ax = axes[r, 1]
        for i, label in enumerate(cell_types):
            mask = labels == label
            ax.scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=0.7, color=colour_map(i % 10), label=label, linewidths=0)
        ax.set_title(f'{method}: cell type; mix={mix:.3f}, sil={sil:.3f}, ratio={ratio:.3f}')
        ax.legend(markerscale=2, fontsize=7, ncol=2); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(OUT / 'fig_umap.png', dpi=180); plt.close(fig)
    print(json.dumps({'figure': 'fig_umap.png', 'methods': list(METHODS)}, sort_keys=True))

def report_stage():
    gates = get_gates()
    stopped = gates.get('status') == 'stopped'
    files = ['requirements.lock', 'embeddings.npz', 'scores.csv', 'scripts/t05a_run.py']
    if (OUT / 'centroid_ops.npz').exists(): files += ['centroid_ops.npz', 'frames.csv', 'fig_umap.png']
    command = 'cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T05a/scripts/t05a_run.py'
    score_keys = tuple([f'mix_{method}' for method in METHODS] + [f'sil_{method}' for method in METHODS] + ['ours_repro_max_abs'])
    operator_keys = tuple([f'frame_ratio_{method}' for method in METHODS] + [f'frame_outcome_{method}' for method in METHODS] + ['n_clusters', 'symmetry_max'])
    regeneration = {path: command + ' --stage all' for path in files}
    regeneration.update({key: command + ' --stage scores' for key in score_keys if key in gates})
    regeneration.update({key: command + ' --stage operators' for key in operator_keys if key in gates})
    report = {
        'task_id': 'T05a', 'executor': 'manus', 'branch': 'task/T05a',
        'commit': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
        'seed': 20261014, 'started': '2026-10-15T16:17:00-06:00', 'finished': datetime.now(timezone.utc).isoformat(),
        'status': 'stopped' if stopped else 'complete',
        'gated': {key: value for key, value in gates.items() if key not in {'status', 'validity_revise'}},
        'deliverables': [{'path': path, 'sha256': sha(OUT / path)} for path in files],
        'regeneration': regeneration,
        'deviations': [{'what': 'Installed harmonypy==0.0.10 in the project virtual environment.', 'why': 'Authorised explicitly by tasks/T05a.md.'}],
        'deviations_requested': [],
        'observations': [f'Input SHA-256: {path.relative_to(ROOT)} = {sha(path)}.' for path in INPUTS],
        'known_failures': KNOWN, 'lockfile_sha256': sha(OUT / 'requirements.lock'),
    }
    if stopped: report['stop_rule'] = 'n_clusters < 3 → stopped.'
    if gates.get('ours_repro_max_abs', 0.0) > 1e-6 or gates.get('symmetry_max', 0.0) > 1e-4:
        report['validity'] = 'REVISE'
    (OUT / 'report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(report['gated'], sort_keys=True))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', required=True, choices=('embed', 'scores', 'operators', 'views', 'report', 'all'))
    stage = parser.parse_args().stage
    if stage in ('embed', 'all'): embed_stage()
    if stage in ('scores', 'all'): scores_stage()
    if stage in ('operators', 'all'): operators_stage()
    if stage in ('views', 'all'): views_stage()
    if stage in ('report', 'all'): report_stage()
if __name__ == '__main__': main()
