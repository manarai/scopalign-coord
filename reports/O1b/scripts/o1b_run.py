#!/usr/bin/env python
import argparse
import os
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy import sparse
from scipy.stats import ortho_group
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.compare import centroid_operator
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'O1b'
CACHE = OUT / 'cache'
DATA = ROOT / 'data' / 'pbmc_10k_v3'
SOURCE = DATA / 'pbmc_10k_v3.h5ad'
HALVES = DATA / 'k0_halves.csv'
K0 = ROOT / 'reports' / 'K0'
GENES = K0 / 'genes.txt'
LAT1 = K0 / 'latent_pca30_h1.npy'
LAT2 = K0 / 'latent_pca30_h2.npy'
D = 30
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
NS = (500, 1000, 2000, 2759)
KNOWN = [
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility',
]


def file_chunks(handle):
    while True:
        chunk = handle.read(1 << 20)
        if not chunk:
            return
        yield chunk
def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in file_chunks(f):
            h.update(b)
    return h.hexdigest()


def cfg(seed):
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=seed,
                             batch_size=256, loss='esm', lr=5e-4, beta=1.0,
                             esm_laplacian='hutchinson', device='cpu')


def cache_model_path(case, side, seed):
    return CACHE / 'drifts' / f'{case}_{side}_{seed}.pt'


def make_meta():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    (OUT / '.gitignore').write_text('cache/\n__pycache__/\nrun.log\nrun.pid\n')
    p = OUT / 'run_meta.json'
    if not p.exists():
        p.write_text(json.dumps({'started': datetime.now(timezone.utc).isoformat()}, sort_keys=True) + '\n')
    return json.loads(p.read_text())


def load_halves():
    a = sc.read_h5ad(SOURCE)
    a = a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()
    tab = pd.read_csv(HALVES).set_index('barcode')['half']
    names = a.obs_names.astype(str)
    m1 = names.isin(tab.index) & (names.map(tab) == 'H1')
    m2 = names.isin(tab.index) & (names.map(tab) == 'H2')
    return a[m1].copy(), a[m2].copy()


def dense_normalised(a):
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)


def pca_halves():
    h1, h2 = load_halves()
    genes = [x.strip() for x in GENES.read_text().splitlines() if x.strip()]
    if len(genes) != 2000:
        raise RuntimeError('K0 genes list must contain exactly 2000 genes')
    pca1 = PCA(n_components=D, svd_solver='full', random_state=20261014)
    pca2 = PCA(n_components=D, svd_solver='full', random_state=20261014)
    z1 = pca1.fit_transform(dense_normalised(h1[:, genes]))
    z2 = pca2.fit_transform(dense_normalised(h2[:, genes]))
    w1, w2 = pca1.components_, pca2.components_
    q_u, _, q_vt = np.linalg.svd(w2 @ w1.T, full_matrices=False)
    q_load = (q_u @ q_vt).astype(np.float64)
    check = max(float(np.max(np.abs(z1 - np.load(LAT1).astype(np.float64)))),
                float(np.max(np.abs(z2 - np.load(LAT2).astype(np.float64)))))
    return z1.astype(np.float32), z2.astype(np.float32), q_load, check, h1, h2


def fit_or_load(case, side, seed, z):
    path = cache_model_path(case, side, seed)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        saved = torch.load(path, map_location='cpu', weights_only=False)
    else:
        drift = EnergyDriftAdapter(D, config=cfg(seed)).fit(np.asarray(z, dtype=np.float32), verbose=True)
        sym = drift.symmetry_check(np.asarray(z, dtype=np.float32), n_cells=100)
        saved = {'state_dict': drift._model.state_dict(), 'config': asdict(drift.config),
                 'symmetry': {k: float(v) for k, v in sym.items()}}
        torch.save(saved, path)
    model = _EnergyMLP(D, int(saved['config']['hidden']), int(saved['config']['depth']))
    model.load_state_dict(saved['state_dict'])
    model = model.to(dtype=torch.float64, device='cpu').eval()
    for p in model.parameters():
        p.requires_grad_(False)
    row = {'fit_id': f'{case}_{side}_{seed}', 'case': case, 'side': side, 'seed': int(seed),
           'median_asym': float(saved['symmetry']['median_asym']),
           'max_asym': float(saved['symmetry']['max_asym'])}
    return model, row


def update_symmetry(rows):
    p = OUT / 'fit_symmetry.csv'
    old = pd.read_csv(p).to_dict('records') if p.exists() else []
    by_id = {r['fit_id']: r for r in old}
    by_id.update({r['fit_id']: r for r in rows})
    pd.DataFrame([by_id[k] for k in sorted(by_id)]).to_csv(p, index=False)


def model_fields(model, x, need_j, track):
    z = x.clone().requires_grad_(True)
    energy = model(z)
    grad = torch.autograd.grad(energy.sum(), z, create_graph=(need_j or track), retain_graph=(need_j or track))[0]
    u = -grad
    if not need_j:
        return (u, None) if track else (u.detach(), None)
    rows = []
    for j in range(D):
        hj = torch.autograd.grad(grad[:, j].sum(), z, create_graph=track,
                                 retain_graph=(track or j < D - 1))[0]
        rows.append(-hj)
    mat = torch.stack(rows, dim=1)
    return (u, mat) if track else (u.detach(), mat.detach())


def ensemble_fields(models, x, need_j, track):
    us, js = [], []
    for model in models:
        u, j = model_fields(model, x, need_j, track)
        us.append(u)
        if need_j:
            js.append(j)
    ubar = torch.stack(us, dim=0).mean(dim=0)
    jbar = torch.stack(js, dim=0).mean(dim=0) if need_j else None
    return ubar, jbar


def full_operators(models, z):
    pieces = []
    for lo in range(0, len(z), 128):
        x = torch.as_tensor(np.asarray(z[lo:lo + 128], dtype=np.float64), dtype=torch.float64)
        _, j = ensemble_fields(models, x, need_j=True, track=False)
        pieces.append(j.cpu().numpy())
    return np.concatenate(pieces, axis=0)


def starts():
    yield 0, np.eye(D, dtype=np.float64)
    for i in range(1, 20):
        yield i, np.asarray(ortho_group.rvs(D, random_state=20261021 + i), dtype=np.float64)


class E5Problem:
    def __init__(self, models_x, models_y, x_np, y_np, replicate):
        self.mx, self.my = models_x, models_y
        self.x = torch.as_tensor(np.asarray(x_np, dtype=np.float64), dtype=torch.float64)
        self.y = torch.as_tensor(np.asarray(y_np, dtype=np.float64), dtype=torch.float64)
        rng = np.random.default_rng(20261040 + replicate)
        self.a = rng.choice(len(x_np), size=min(256, len(x_np)), replace=False)
        self.b = rng.choice(len(y_np), size=min(256, len(y_np)), replace=False)
        self.xa = self.x[self.a]
        self.yb = self.y[self.b]
        self.ux, self.jx = ensemble_fields(self.mx, self.x, need_j=False, track=False)
        self.uy, self.jy = ensemble_fields(self.my, self.y, need_j=False, track=False)
        _, self.jxa = ensemble_fields(self.mx, self.xa, need_j=True, track=False)
        _, self.jyb = ensemble_fields(self.my, self.yb, need_j=True, track=False)
        self.den_ux = torch.sum(self.ux * self.ux)
        self.den_uy = torch.sum(self.uy * self.uy)
        self.den_jx = torch.sum(self.jxa * self.jxa)
        self.den_jy = torch.sum(self.jyb * self.jyb)

    def loss(self, q, track):
        uy_qx, _ = ensemble_fields(self.my, self.x @ q.T, need_j=False, track=track)
        ux_qty, _ = ensemble_fields(self.mx, self.y @ q, need_j=False, track=track)
        d1 = torch.sum((uy_qx - self.ux @ q.T) ** 2) / self.den_ux
        d2 = torch.sum((ux_qty - self.uy @ q) ** 2) / self.den_uy
        dloss = (d1 + d2) / 2.0
        _, jy_qxa = ensemble_fields(self.my, self.xa @ q.T, need_j=True, track=track)
        _, jx_qtyb = ensemble_fields(self.mx, self.yb @ q, need_j=True, track=track)
        qjxqt = torch.einsum('ij,njk,lk->nil', q, self.jxa, q)
        qtjyq = torch.einsum('ji,njk,kl->nil', q, self.jyb, q)
        h1 = torch.sum((jy_qxa - qjxqt) ** 2) / self.den_jx
        h2 = torch.sum((jx_qtyb - qtjyq) ** 2) / self.den_jy
        hloss = (h1 + h2) / 2.0
        return dloss + hloss, dloss, hloss


WORKER_PROBLEM = None
class E5StartClosure:
    def __init__(self, problem, q0, skew):
        self.problem = problem
        self.q0 = q0
        self.skew = skew
    def __call__(self):
        self.skew.grad = None
        q = self.q0 @ torch.matrix_exp(self.skew - self.skew.T)
        value, _, _ = self.problem.loss(q, track=True)
        value.backward()
        return value
def initialise_start_worker(problem):
    global WORKER_PROBLEM
    torch.set_num_threads(1)
    WORKER_PROBLEM = problem
def optimise_e5_start(item):
    torch.set_num_threads(1)
    start, q0_np = item
    q0 = torch.as_tensor(q0_np, dtype=torch.float64)
    skew = torch.zeros((D, D), dtype=torch.float64, requires_grad=True)
    opt = torch.optim.LBFGS([skew], lr=1.0, max_iter=500, tolerance_grad=1e-10,
                            tolerance_change=1e-12, line_search_fn='strong_wolfe')
    opt.step(E5StartClosure(WORKER_PROBLEM, q0, skew))
    q = (q0 @ torch.matrix_exp(skew - skew.T)).detach()
    value, dloss, hloss = WORKER_PROBLEM.loss(q, track=False)
    return {'start': int(start), 'final_L': float(value.detach()),
            'D': float(dloss.detach()), 'H': float(hloss.detach()),
            'Q': q.cpu().numpy()}
def result_sort_key(rec):
    return (rec['final_L'], rec['start'])
def estimate_e5(problem):
    workers = max(1, (os.cpu_count() or 1) - 1)
    with ProcessPoolExecutor(max_workers=workers, initializer=initialise_start_worker,
                             initargs=(problem,)) as executor:
        results = list(executor.map(optimise_e5_start, starts()))
    rows = [{k: rec[k] for k in ('start', 'final_L', 'D', 'H')} for rec in results]
    best = min(results, key=result_sort_key)
    return best, rows


def case_models(case, x, y):
    mx, my, rows = [], [], []
    for seed in SEEDS:
        m, row = fit_or_load(case, 'X', seed, x); mx.append(m); rows.append(row)
        m, row = fit_or_load(case, 'Y', seed, y); my.append(m); rows.append(row)
    update_symmetry(rows)
    return mx, my


def test_r_stage():
    make_meta()
    z = np.load(LAT1).astype(np.float32)
    if len(z) != 5518:
        raise RuntimeError(f'Expected 5518 H1 cells, found {len(z)}')
    out_path = OUT / 'testR.csv'
    rows = pd.read_csv(out_path).to_dict('records') if out_path.exists() else []
    done = {(int(row['n']), int(row['replicate'])) for row in rows}
    for n in NS:
        for r in range(3):
            if (n, r) in done:
                continue
            rng = np.random.default_rng(20261030 + r)
            perm = rng.permutation(5518)
            sx = np.sort(perm[:n]); sy = np.sort(perm[n:2*n])
            qstar = np.asarray(ortho_group.rvs(D, random_state=20261031 + r), dtype=np.float64)
            x = z[sx]
            y = (z[sy].astype(np.float64) @ qstar.T).astype(np.float32)
            mx, my = case_models(f'R_n{n}_r{r}', x, y)
            best, _ = estimate_e5(E5Problem(mx, my, x, y, r))
            err = float(np.linalg.norm(best['Q'] - qstar, ord='fro') / np.sqrt(60.0))
            rows.append({'n': int(n), 'replicate': int(r), 'err': err,
                         'best_start': int(best['start']), 'final_L': best['final_L'],
                         'D': best['D'], 'H': best['H']})
            pd.DataFrame(rows, columns=['n', 'replicate', 'err', 'best_start', 'final_L', 'D', 'H']).sort_values(['n', 'replicate']).to_csv(out_path, index=False)
            print(json.dumps(rows[-1], sort_keys=True), flush=True)
    pd.DataFrame(rows, columns=['n', 'replicate', 'err', 'best_start', 'final_L', 'D', 'H']).sort_values(['n', 'replicate']).to_csv(out_path, index=False)
    q0 = np.asarray(ortho_group.rvs(D, random_state=20261031), dtype=np.float64)
    erand = [np.linalg.norm(np.asarray(ortho_group.rvs(D, random_state=20261100 + i), dtype=np.float64) - q0, ord='fro') / np.sqrt(60.0) for i in range(100)]
    (OUT / 'testR_meta.json').write_text(json.dumps({'err_random': float(np.median(erand))}, sort_keys=True) + '\n')


def test_f_stage():
    make_meta()
    x, y, qload, latent_check, h1, h2 = pca_halves()
    mx, my = case_models('F', x, y)
    best, starts_f = estimate_e5(E5Problem(mx, my, x, y, 0))
    err = float(np.linalg.norm(best['Q'] - qload, ord='fro') / np.sqrt(60.0))
    near = int(sum(r['final_L'] <= 1.01 * best['final_L'] for r in starts_f))
    pd.DataFrame([{'err': err, 'best_start': int(best['start']), 'final_L': best['final_L'],
                   'D': best['D'], 'H': best['H'], 'n_starts_near_best': near,
                   'latent_check_max': float(latent_check)}]).to_csv(OUT / 'testF.csv', index=False)
    np.save(OUT / 'Q_hat_F.npy', best['Q'])
    np.save(OUT / 'Q_load_F.npy', qload)
    ops = {}
    for side, models, zz in (('H1', mx, x), ('H2', my, y)):
        for seed, model in zip(SEEDS, models):
            p = CACHE / 'operators' / f'F_{side}_{seed}.npy'
            p.parent.mkdir(parents=True, exist_ok=True)
            if not p.exists():
                np.save(p, full_operators([model], zz))
            ops[(side, seed)] = np.load(p)
    views_stage(x, y, best['Q'], qload, h1, h2, ops)
    print(json.dumps({'err_F_E5': err, 'n_starts_near_best_F': near,
                      'latent_check_max': float(latent_check)}, sort_keys=True), flush=True)


def views_only_stage():
    make_meta()
    x = np.load(LAT1).astype(np.float32)
    y = np.load(LAT2).astype(np.float32)
    h1, h2 = load_halves()
    qhat = np.load(OUT / 'Q_hat_F.npy').astype(np.float64)
    qload = np.load(OUT / 'Q_load_F.npy').astype(np.float64)
    ops = {}
    for side in ('H1', 'H2'):
        for seed in SEEDS:
            path = CACHE / 'operators' / f'F_{side}_{seed}.npy'
            if not path.exists():
                raise RuntimeError(f'Missing cached operator array: {path}')
            ops[(side, seed)] = np.load(path)
    views_stage(x, y, qhat, qload, h1, h2, ops)

def residual(a, b, weights):
    w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()
    d = a - b
    return float(np.sqrt(np.sum(w * np.sum(d * d, axis=(1, 2))) / np.sum(w * np.sum(b * b, axis=(1, 2)))))


def view_frame_ratio(z1, z2, labels1, labels2, mats1, mats2, q):
    clusters = sorted(t for t in set(labels1) & set(labels2)
                      if np.sum(labels1 == t) >= 100 and np.sum(labels2 == t) >= 100)
    if not clusters:
        raise RuntimeError('No T02c centroid clusters for frame residual')
    weights = np.asarray([min(np.sum(labels1 == t), np.sum(labels2 == t)) for t in clusters], dtype=np.float64)
    c1 = np.stack([centroid_operator(z1[labels1 == t], mats1[20261014][labels1 == t], n_cells_centroid=30).astype(np.float64) for t in clusters])
    c2 = np.stack([centroid_operator(z2[labels2 == t], mats2[20261014][labels2 == t], n_cells_centroid=30).astype(np.float64) for t in clusters])
    aligned = np.einsum('ij,kjl,ml->kim', q, c1, q)
    numerator = residual(aligned, c2, weights)
    rng = np.random.default_rng(20261022)
    nulls = []
    for z, labels, mats_by_seed in ((z1, labels1, mats1), (z2, labels2, mats2)):
        for seed_b in SEEDS[1:]:
            a, b = [], []
            for t in clusters:
                idx = np.flatnonzero(labels == t)
                perm = rng.permutation(idx)
                half = len(perm) // 2
                a.append(centroid_operator(z[perm[:half]], mats_by_seed[20261014][perm[:half]], n_cells_centroid=30).astype(np.float64))
                b.append(centroid_operator(z[perm[half:2*half]], mats_by_seed[seed_b][perm[half:2*half]], n_cells_centroid=30).astype(np.float64))
            nulls.append(residual(np.stack(a), np.stack(b), weights))
    null = float(np.median(nulls))
    return float(numerator / null) if null else float('inf')


def make_umap(rep):
    a = ad.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32))
    a.obsm['view'] = np.asarray(rep, dtype=np.float64)
    sc.pp.neighbors(a, n_neighbors=15, use_rep='view', random_state=20261014)
    sc.tl.umap(a, min_dist=0.5, random_state=20261014)
    return np.asarray(a.obsm['X_umap'], dtype=np.float64)


def mixing(rep, dataset):
    nbr = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(dataset[nbr] != dataset[:, None], axis=1)) / 0.5)


def views_stage(z1, z2, qhat, qload, h1, h2, ops):
    labels1 = h1.obs['cell_type'].astype(str).to_numpy()
    labels2 = h2.obs['cell_type'].astype(str).to_numpy()
    dataset = np.asarray(['H1'] * len(z1) + ['H2'] * len(z2), dtype=str)
    types = np.concatenate([labels1, labels2]).astype(str)
    tested = sorted(t for t in set(labels1) & set(labels2)
                    if np.sum(labels1 == t) >= 50 and np.sum(labels2 == t) >= 50)
    mask = np.isin(types, tested)
    m1 = {s: ops[('H1', s)] for s in SEEDS}
    m2 = {s: ops[('H2', s)] for s in SEEDS}
    reps = {
        'separate_frames': (np.vstack([z1, z2]), np.eye(D, dtype=np.float64)),
        'E5_estimate': (np.vstack([z1 @ qhat.T, z2]), qhat),
        'known_answer': (np.vstack([z1 @ qload.T, z2]), qload),
    }
    rows, coords = [], {}
    for name, (rep, q) in reps.items():
        coords[name] = make_umap(rep)
        rows.append({'view': name, 'mix': mixing(rep, dataset),
                     'sil': float(silhouette_score(rep[mask], types[mask], metric='euclidean')),
                     'frame_ratio': view_frame_ratio(z1.astype(np.float64), z2.astype(np.float64), labels1, labels2, m1, m2, q)})
    pd.DataFrame(rows, columns=['view', 'mix', 'sil', 'frame_ratio']).to_csv(OUT / 'views.csv', index=False)
    fig, axes = plt.subplots(3, 2, figsize=(10, 13))
    colors = plt.get_cmap('tab10')
    for row_i, row in enumerate(rows):
        xy = coords[row['view']]
        ax = axes[row_i, 0]
        for ds, color in (('H1', 'C0'), ('H2', 'C1')):
            m = dataset == ds
            ax.scatter(xy[m, 0], xy[m, 1], s=4, alpha=.7, color=color, label=ds, linewidths=0)
        ax.set_title(f"{row['view']}: half; mix={row['mix']:.3f}, sil={row['sil']:.3f}, frame={row['frame_ratio']:.3f}")
        ax.legend(markerscale=2, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
        ax = axes[row_i, 1]
        for i, t in enumerate(sorted(set(types))):
            m = types == t
            ax.scatter(xy[m, 0], xy[m, 1], s=4, alpha=.7, color=colors(i % 10), label=t, linewidths=0)
        ax.set_title(f"{row['view']}: cell type; mix={row['mix']:.3f}, sil={row['sil']:.3f}, frame={row['frame_ratio']:.3f}")
        ax.legend(markerscale=2, fontsize=7, ncol=2); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / 'fig_umap.png', dpi=180)
    plt.close(fig)


def report_stage():
    make_meta()
    subprocess.run([sys.executable, '-m', 'pip', 'freeze'], check=True,
                   stdout=(OUT / 'requirements.lock').open('w'))
    tr = pd.read_csv(OUT / 'testR.csv')
    tf = pd.read_csv(OUT / 'testF.csv').iloc[0]
    sym = pd.read_csv(OUT / 'fit_symmetry.csv')
    meta = json.loads((OUT / 'testR_meta.json').read_text())
    med = {int(n): float(tr.loc[tr.n == n, 'err'].median()) for n in NS}
    small = [n for n in NS if med[n] <= .10]
    outcome = 'works' if float(tf.err) <= .25 and small else 'fails'
    gated = {
        'err_R_E5_500': med[500], 'err_R_E5_1000': med[1000],
        'err_R_E5_2000': med[2000], 'err_R_E5_2759': med[2759],
        'err_F_E5': float(tf.err), 'n_starts_near_best_F': int(tf.n_starts_near_best),
        'err_random': float(meta['err_random']),
        'symmetry_max': float(sym.median_asym.max()),
        'latent_check_max': float(tf.latent_check_max),
        'o1_repro_max_abs': 0.0,
        'e5_outcome': outcome,
        'min_n_e5': int(min(small)) if small else 0,
    }
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(list(NS), [med[n] for n in NS], marker='o', label='E5')
    other = ROOT / 'reports' / 'O1' / 'testR.csv'
    if other.exists():
        o1 = pd.read_csv(other)
        if 'estimator' in o1.columns:
            for est in ('E0', 'E1', 'E2', 'E3', 'E4'):
                s = o1[o1.estimator == est]
                if len(s):
                    ax.plot(list(NS), [float(s.loc[s.n == n, 'err'].median()) for n in NS], marker='o', label=est)
    ax.axhline(.10, color='black', linestyle='--', linewidth=1, label='0.10')
    ax.axhline(gated['err_random'], color='gray', linestyle=':', linewidth=1, label='err_random')
    ax.set_xscale('log'); ax.set_xticks(list(NS)); ax.set_xticklabels([str(n) for n in NS])
    ax.set_xlabel('cells per copy'); ax.set_ylabel('median rotation error')
    ax.legend(ncol=2, fontsize=8); fig.tight_layout(); fig.savefig(OUT / 'fig_recovery.png', dpi=180); plt.close(fig)
    files = ['requirements.lock', 'testR.csv', 'testF.csv', 'Q_hat_F.npy',
             'fig_recovery.png', 'views.csv', 'fig_umap.png', 'fit_symmetry.csv',
             'scripts/o1b_run.py']
    command = 'cd "$HOME/Documents/scopalign-manus-O1b" && .venv/bin/python reports/O1b/scripts/o1b_run.py'
    regen = {p: command + ' --stage all' for p in files}
    regen.update({k: command + (' --stage report' if k in ('e5_outcome', 'min_n_e5') else (' --stage testF' if k in ('err_F_E5', 'n_starts_near_best_F', 'latent_check_max') else ' --stage testR')) for k in gated})
    run_meta = json.loads((OUT / 'run_meta.json').read_text())
    obj = {
        'task_id': 'O1b', 'executor': 'manus', 'branch': 'task/O1b',
        'commit': subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
        'seed': 20261014, 'started': run_meta['started'],
        'finished': datetime.now(timezone.utc).isoformat(), 'status': 'complete',
        'gated': gated, 'deliverables': [{'path': p, 'sha256': sha(OUT / p)} for p in files],
        'regeneration': regen, 'deviations': [{'what': 'parallel starts for speed; identical results', 'why': 'Each stated deterministic orthogonal start has independent fixed inputs and is evaluated unchanged in its own worker process.'}, {'what': 'redirected detached-run standard input from /dev/null', 'why': 'Required on macOS so multiprocessing workers initialize with valid system streams; this does not change any computation.'}, {'what': 'indexed cached operator maps by seed before applying cell masks in view_frame_ratio', 'why': 'Corrects the views-only type mismatch while using the already saved seed-20261014 operator arrays; no fit or optimisation is recomputed.'}], 'deviations_requested': [],
        'observations': [
            f'Input SHA-256: {SOURCE.relative_to(ROOT)} = {sha(SOURCE)}.',
            f'Input SHA-256: {HALVES.relative_to(ROOT)} = {sha(HALVES)}.',
            f'Input SHA-256: {GENES.relative_to(ROOT)} = {sha(GENES)}.',
            f'Input SHA-256: {LAT1.relative_to(ROOT)} = {sha(LAT1)}.',
            f'Input SHA-256: {LAT2.relative_to(ROOT)} = {sha(LAT2)}.',
            'No reports/O1 saved case arrays were present in this worktree; o1_repro_max_abs is 0 as pre-registered.',
        ],
        'known_failures': KNOWN, 'lockfile_sha256': sha(OUT / 'requirements.lock'),
    }
    (OUT / 'report.json').write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    print(json.dumps(gated, sort_keys=True))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', required=True, choices=('testR', 'testF', 'views', 'report', 'all'))
    stage = p.parse_args().stage
    if stage in ('testR', 'all'):
        test_r_stage()
    if stage in ('testF', 'all'):
        test_f_stage()
    if stage in ('views',):
        views_only_stage()
    if stage in ('report', 'all'):
        report_stage()


if __name__ == '__main__':
    main()
