#!/usr/bin/env python
"""O1o: fixed 10-PC E5s comparison on K2 and four within-age donor pairs.

Run --stage all detached with stdin /dev/null. Input data and the named
reference worktrees are read-only. All prepared inputs, 50 fitted drifts,
pair checkpoints, and runtime caches are under reports/O1o. Each completed
stage commits only reports/O1o and pushes task/O1o to the local origin.
The E5s numerical helpers below are copied verbatim from O1n's runner.
"""
import argparse
import fcntl
import hashlib
import json
import multiprocessing
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'O1o'
CACHE = OUT / 'cache'
ARTIFACTS = OUT / 'artifacts'
sys.dont_write_bytecode = True
for key, name in (('MPLCONFIGDIR', '.mplconfig'), ('XDG_CACHE_HOME', '.cache'), ('NUMBA_CACHE_DIR', '.numba')):
    target = OUT / name
    target.mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(target)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP

BASE = Path.home() / 'Documents' / 'scopalign-manus'
K2 = BASE / 'reports' / 'K2'
T02B = BASE / 'reports' / 'T02b'
V2 = BASE / 'data' / 'pbmc_v2_v3' / 'pbmc_v2.h5ad'
V3THIN = BASE / 'data' / 'pbmc_v2_v3' / 'pbmc_v3_thin.h5ad'
K2CELLS = BASE / 'data' / 'pbmc_v2_v3' / 'k2_cells.csv'
K2GENES = BASE / 'reports' / 'T01d' / 'genes.txt'
AGING = Path('/Users/terooatt/Documents/ChatGPT/Bio559R_studio0/aging_pbmc_8donors_prepared.h5ad')
AGING_SHA256 = '43e39835f7aa070e6e1843a712c74aa91018ae9f887ea9e1bd18e65ba2049dd4'
REFERENCES = [
    BASE.parent / 'scopalign-manus-O1g' / 'reports' / 'O1g' / 'scripts' / 'o1g_run.py',
    BASE.parent / 'scopalign-manus-O1l' / 'reports' / 'O1l' / 'scripts' / 'o1l_run.py',
    BASE.parent / 'scopalign-manus-O1n' / 'reports' / 'O1n' / 'scripts' / 'o1n_run.py',
    K2 / 'scripts' / 'k2_run.py', T02B / 'scripts' / 't02b_run.py',
]
D = 10
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
DONORS = ('F30', 'F36', 'F66', 'F77', 'M28', 'M36', 'M60', 'M66')
PAIR_SPECS = (('K2', 'K2_A', 'K2_B'), ('F30F36', 'F30', 'F36'),
              ('M28M36', 'M28', 'M36'), ('F66F77', 'F66', 'F77'), ('M60M66', 'M60', 'M66'))
METHODS = ('noA', 'skew', 'E5s', 'gene', 'harmony')
MAIN_TORCH_THREADS = torch.get_num_threads()
KNOWN_FAILURES = [
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility',
]
GATED_KEYS = ('symmetry_max',) + tuple(
    f'{name}_{dataset}' for dataset in ('K2', 'A') for name in
    ('mix_noA', 'mix_E5s', 'mix_skew', 'mix_gene', 'mix_harmony', 'sil_noA', 'sil_E5s', 'err_gene')
) + ('n_pass_A', 'o1o_outcome')
DEVIATIONS = [{
    'what': 'Used report-local Matplotlib/fontconfig/Numba caches and O1n\'s in-process start executor with three workers.',
    'why': 'Shared environments and reference worktrees are read-only; O1n uses threads to avoid the local PyTorch shared-memory restriction. Copied E5s loss, starts, screening, and LBFGS settings are unchanged.',
}]
PAIR_COLUMNS = ('pair', 'method', 'mix', 'sil', 'err_gene', 'L_hat', 'L_gene', 'screen_rank_best')

# BEGIN_COPIED_O1N_HELPERS
def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def array_sha(value):
    value = np.ascontiguousarray(value)
    h = hashlib.sha256()
    h.update(str(value.dtype).encode())
    h.update(str(value.shape).encode())
    h.update(value.tobytes())
    return h.hexdigest()


def temp_path(path, suffix='.tmp'):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f'.{path.name}.', suffix=suffix, delete=False) as handle:
        return Path(handle.name)


def atomic_json(path, value):
    tmp = temp_path(path)
    try:
        tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_csv(path, frame):
    tmp = temp_path(path)
    try:
        frame.to_csv(tmp, index=False)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_npz(path, **arrays):
    tmp = temp_path(path, '.npz')
    try:
        np.savez(tmp, **arrays)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_model(path, saved):
    tmp = temp_path(path, '.pt')
    try:
        torch.save(saved, tmp)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def log(event, **values):
    print(json.dumps({'event': event, **values}, sort_keys=True, allow_nan=False), flush=True)


def cfg(seed):
    # O1/O1b config unchanged except EnergyDriftAdapter's input dimension.
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=int(seed),
                             batch_size=256, loss='esm', lr=5e-4, beta=1.0,
                             esm_laplacian='hutchinson', device='cpu')


def build_model(saved, d):
    model = _EnergyMLP(d, int(saved['config']['hidden']), int(saved['config']['depth']))
    model.load_state_dict(saved['state_dict'])
    model = model.to(dtype=torch.float64, device='cpu').eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model


def fields(model, x, need_j, track, d):
    # O1l's dynamic-dimension copy of O1b's model_fields, unchanged.
    z = x.clone().requires_grad_(True)
    e = model(z)
    g = torch.autograd.grad(e.sum(), z, create_graph=(need_j or track), retain_graph=(need_j or track))[0]
    u = -g
    if not need_j:
        return (u, None) if track else (u.detach(), None)
    rows = []
    for j in range(d):
        rows.append(-torch.autograd.grad(g[:, j].sum(), z, create_graph=track,
                                        retain_graph=(track or j < d - 1))[0])
    jacobian = torch.stack(rows, dim=1)
    return (u, jacobian) if track else (u.detach(), jacobian.detach())


def ensemble(models, x, need_j, track, d):
    us, js = [], []
    for model in models:
        u, j = fields(model, x, need_j, track, d)
        us.append(u)
        if need_j:
            js.append(j)
    return torch.stack(us).mean(0), (torch.stack(js).mean(0) if need_j else None)


class Problem:
    """O1b D+H loss, in O1l's exact dynamic-dimension implementation."""
    def __init__(self, mx, my, x, y, r, d):
        self.mx, self.my, self.d = mx, my, d
        self.x = torch.as_tensor(np.asarray(x, dtype=np.float64), dtype=torch.float64)
        self.y = torch.as_tensor(np.asarray(y, dtype=np.float64), dtype=torch.float64)
        rng = np.random.default_rng(20261040 + r)
        self.xa = self.x[rng.choice(len(x), size=min(256, len(x)), replace=False)]
        self.yb = self.y[rng.choice(len(y), size=min(256, len(y)), replace=False)]
        self.ux, _ = ensemble(mx, self.x, False, False, d)
        self.uy, _ = ensemble(my, self.y, False, False, d)
        _, self.jx = ensemble(mx, self.xa, True, False, d)
        _, self.jy = ensemble(my, self.yb, True, False, d)
        self.du = torch.sum(self.ux * self.ux)
        self.dv = torch.sum(self.uy * self.uy)
        self.dx = torch.sum(self.jx * self.jx)
        self.dy = torch.sum(self.jy * self.jy)

    def loss(self, q, track):
        uy, _ = ensemble(self.my, self.x @ q.T, False, track, self.d)
        ux, _ = ensemble(self.mx, self.y @ q, False, track, self.d)
        dl = (torch.sum((uy - self.ux @ q.T) ** 2) / self.du + torch.sum((ux - self.uy @ q) ** 2) / self.dv) / 2
        _, jy = ensemble(self.my, self.xa @ q.T, True, track, self.d)
        _, jx = ensemble(self.mx, self.yb @ q, True, track, self.d)
        qj = torch.einsum('ij,njk,lk->nil', q, self.jx, q)
        qtj = torch.einsum('ji,njk,kl->nil', q, self.jy, q)
        hl = (torch.sum((jy - qj) ** 2) / self.dx + torch.sum((jx - qtj) ** 2) / self.dy) / 2
        return dl + hl, dl, hl


def screen_losses(problem, q_values, q_batch=4):
    # O1l's bounded batched screen, unchanged except explicit dimension.
    q_values = torch.as_tensor(q_values, dtype=torch.float64)
    values = []
    for begin in range(0, len(q_values), q_batch):
        q = q_values[begin:begin + q_batch]
        b, d = len(q), problem.d
        nx, ny, nxa, nyb = len(problem.x), len(problem.y), len(problem.xa), len(problem.yb)
        uy, _ = ensemble(problem.my, torch.einsum('nd,bkd->bnk', problem.x, q).reshape(b * nx, d), False, False, d)
        ux, _ = ensemble(problem.mx, torch.einsum('nd,bdk->bnk', problem.y, q).reshape(b * ny, d), False, False, d)
        left = uy.reshape(b, nx, d) - torch.einsum('nd,bkd->bnk', problem.ux, q)
        right = ux.reshape(b, ny, d) - torch.einsum('nd,bdk->bnk', problem.uy, q)
        dl = (torch.sum(left * left, dim=(1, 2)) / problem.du + torch.sum(right * right, dim=(1, 2)) / problem.dv) / 2
        _, jy = ensemble(problem.my, torch.einsum('nd,bkd->bnk', problem.xa, q).reshape(b * nxa, d), True, False, d)
        _, jx = ensemble(problem.mx, torch.einsum('nd,bdk->bnk', problem.yb, q).reshape(b * nyb, d), True, False, d)
        qj = torch.einsum('bij,njk,blk->bnil', q, problem.jx, q)
        qtj = torch.einsum('bji,njk,bkl->bnil', q, problem.jy, q)
        hl = (torch.sum((jy.reshape(b, nxa, d, d) - qj) ** 2, dim=(1, 2, 3)) / problem.dx +
              torch.sum((jx.reshape(b, nyb, d, d) - qtj) ** 2, dim=(1, 2, 3)) / problem.dy) / 2
        values.append((dl + hl).detach().numpy())
    return np.concatenate(values)


class Closure:
    def __init__(self, problem, q0, skew):
        self.problem, self.q0, self.skew = problem, q0, skew

    def __call__(self):
        self.skew.grad = None
        q = self.q0 @ torch.matrix_exp(self.skew - self.skew.T)
        value, _, _ = self.problem.loss(q, True)
        value.backward()
        return value


def optimise_start(problem, item):
    index, q_numpy = item
    q0 = torch.as_tensor(q_numpy, dtype=torch.float64)
    skew = torch.zeros((problem.d, problem.d), dtype=torch.float64, requires_grad=True)
    optimiser = torch.optim.LBFGS([skew], lr=1.0, max_iter=500, tolerance_grad=1e-10,
                                 tolerance_change=1e-12, line_search_fn='strong_wolfe')
    optimiser.step(Closure(problem, q0, skew))
    q = (q0 @ torch.matrix_exp(skew - skew.T)).detach()
    value, _, _ = problem.loss(q, False)
    return {'start': int(index), 'L_hat': float(value), 'Q': q.numpy()}


def axes_and_skew(x, y):
    ex, ux = np.linalg.eigh(np.cov(x, rowvar=False))
    ey, uy = np.linalg.eigh(np.cov(y, rowvar=False))
    ex, ux, ey, uy = ex[::-1], ux[:, ::-1], ey[::-1], uy[:, ::-1]
    # Third *central* moment in each copy's own axis frame, no pairing/labels.
    xc = (np.asarray(x, dtype=np.float64) - np.mean(x, axis=0, dtype=np.float64)) @ ux
    yc = (np.asarray(y, dtype=np.float64) - np.mean(y, axis=0, dtype=np.float64)) @ uy
    m3x, m3y = np.mean(xc ** 3, axis=0), np.mean(yc ** 3, axis=0)
    signs = np.sign(m3x) * np.sign(m3y)
    return ux, uy, uy @ np.diag(signs) @ ux.T, ex, ey, m3x, m3y


def sign_start(index, ux, uy, d):
    signs = np.asarray([1.0 if ((index >> j) & 1) == 0 else -1.0 for j in range(d)])
    return uy @ np.diag(signs) @ ux.T


def mixing(rep, datasets):
    # O1g's exact fraction of the 30 nearest OTHER cells with opposite copy,
    # divided by 0.5. This is not a cross-copy-only kNN search.
    neighbors = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(datasets[neighbors] != datasets[:, None], axis=1)) / 0.5)
# END_COPIED_O1N_HELPERS


def atomic_text(path, text):
    tmp = temp_path(path)
    try:
        tmp.write_text(text)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def make_meta():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(exist_ok=True)
    ARTIFACTS.mkdir(exist_ok=True)
    atomic_text(OUT / '.gitignore', 'cache/\n__pycache__/\nrun.log\nrun.pid\nrun.lock\n.mplconfig/\n.cache/\n.numba/\n')
    path = OUT / 'run_meta.json'
    if not path.exists():
        atomic_json(path, {'started': datetime.now(timezone.utc).isoformat()})
    return json.loads(path.read_text())


def copy_path(copy):
    return CACHE / 'inputs' / f'{copy}.npz'


def read_copy(copy, include_log=False):
    with np.load(copy_path(copy), allow_pickle=False) as data:
        return {key: np.array(data[key]) for key in data.files if include_log or key != 'logX'}


def pair_input_path(pair):
    return CACHE / 'inputs' / f'pair_{pair}.npz'


def coarse_aging(label):
    if label.startswith(('CD14-MC', 'Intermed-MC', 'CD16-MC')):
        return 'Mono'
    if label.startswith('NK'):
        return 'NK'
    if label.startswith('CD4'):
        return 'CD4 T'
    if label.startswith('CD8'):
        return 'CD8 T'
    if label.startswith(('Naive-B', 'Memory-B', 'ABC', 'Plasma')):
        return 'B'
    if label.startswith(('cDC', 'pDC')):
        return 'DC'
    if label.startswith(('Megakaryocytes', 'Mk')):
        return 'Mk'
    if label.startswith('RBC'):
        return 'RBC'
    raise RuntimeError(f'Unspecified aging coarse cell-type prefix: {label}')


def readonly_paths():
    return [AGING, V2, V3THIN, K2CELLS, K2GENES,
            K2 / 'latent_pca30_a.npy', K2 / 'latent_pca30_b.npy',
            K2 / 'ops_a.npz', K2 / 'ops_b.npz', T02B / 'loadings.npz', *REFERENCES]


def check_stage():
    make_meta()
    import anndata as ad
    import scanpy as sc
    import harmonypy
    fingerprints = [{'path': str(p), 'sha256': sha(p)} for p in readonly_paths()]
    if fingerprints[0]['sha256'] != AGING_SHA256:
        raise RuntimeError('Aging source SHA-256 differs from the prescribed read-only file')
    path = ARTIFACTS / 'inputs.json'
    if path.exists():
        record = json.loads(path.read_text())
        if record['readonly_files'] != fingerprints:
            raise RuntimeError('A read-only source changed since the first O1o check')
        for copy in ('K2_A', 'K2_B', *DONORS):
            if not copy_path(copy).exists():
                raise RuntimeError(f'Missing prepared copy: {copy}')
        for pair, _, _ in PAIR_SPECS:
            if not pair_input_path(pair).exists():
                raise RuntimeError(f'Missing prepared pair: {pair}')
        if sha(OUT / 'genes.txt') != record['aging_genes_sha256']:
            raise RuntimeError('Prepared aging gene list changed')
        log('check_complete', resumed=True, n_copies=10, n_pairs=5)
        return
    record = {'readonly_files': fingerprints, 'seed': 20261014, 'fit_seeds': list(SEEDS),
              'd': D, 'H_subset_seed_by_pair': {pair: 20261040 for pair, _, _ in PAIR_SPECS},
              'copies': {}, 'pairs': {}, 'aging_cell_choice_definition':
              'sorted donor order; default_rng(20261500+i).choice(n_d,4000,replace=False), preserved in choice order',
              'aging_gene_definition': 'genes alphabetically sorted before batch-aware seurat HVG over all stored log-normalised cells; selected order retained',
              'harmony': {'joint_pca_d': 10, 'svd_solver': 'full', 'random_state': 20261014,
                          'vars_use': 'copy', 'max_iter_harmony': 20},
              'sil_definition': 'coarse labels with >=20 cells in each copy, RBC excluded; full eligible cells, Euclidean silhouette',
              'mix_definition': 'O1n: 30 nearest other cells among stacked copies, fraction opposite-copy / 0.5'}
    kgenes = [gene.strip() for gene in K2GENES.read_text().splitlines() if gene.strip()]
    if len(kgenes) != 2000 or kgenes != sorted(kgenes):
        raise RuntimeError('K2 source genes must be the prescribed sorted 2000 genes')
    cells = pd.read_csv(K2CELLS)
    with np.load(T02B / 'loadings.npz', allow_pickle=False) as loadings:
        wa, wb = np.array(loadings['W_A'][:D]), np.array(loadings['W_B'][:D])
    for copy, letter, source, loading in (('K2_A', 'a', V2, wa), ('K2_B', 'b', V3THIN, wb)):
        a = sc.read_h5ad(source)
        if letter == 'a':
            a = a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()
        bars = cells.loc[cells.dataset == ('A' if letter == 'a' else 'B'), 'barcode'].astype(str).to_numpy()
        a = a[bars].copy()
        labels = a.obs['cell_type'].astype(str).to_numpy(dtype=str)
        with np.load(K2 / f'ops_{letter}.npz', allow_pickle=False) as operators:
            if not np.array_equal(bars, operators['barcode'].astype(str)) or not np.array_equal(labels, operators['cell_type'].astype(str)):
                raise RuntimeError(f'K2 cached barcode/coarse-label ordering mismatch: {copy}')
        a = a[:, kgenes].copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        log_x = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
        log_x = np.asarray(log_x, dtype=np.float32)
        z30 = np.load(K2 / f'latent_pca30_{letter}.npy')
        if z30.shape != (len(bars), 30):
            raise RuntimeError('K2 latent count/dimension mismatch')
        z = z30[:, :D].astype(np.float32)
        atomic_npz(copy_path(copy), z=z, W=loading, logX=log_x, labels=labels, barcodes=bars.astype(str))
        record['copies'][copy] = {'n_cells': len(z), 'z_sha256': array_sha(z), 'W_sha256': array_sha(loading),
                                  'logX_sha256': array_sha(log_x), 'barcodes': bars.tolist(),
                                  'label_counts': {str(k): int(v) for k, v in pd.Series(labels).value_counts().items()},
                                  'genes_sha256': sha(K2GENES)}
    # The source AnnData is loaded read-only and all modifications below are
    # in memory. No .h5ad or cache is ever written beside the source.
    aging = ad.read_h5ad(AGING)
    if not aging.var_names.is_unique:
        raise RuntimeError('Aging var_names are not unique')
    aging = aging[:, np.argsort(aging.var_names.astype(str).to_numpy(), kind='stable')].copy()
    observed_donors = tuple(sorted(set(aging.obs['sample_id'].astype(str))))
    if observed_donors != DONORS:
        raise RuntimeError(f'Unexpected donor set: {observed_donors}')
    sc.pp.highly_variable_genes(aging, n_top_genes=2000, flavor='seurat', batch_key='sample_id')
    genes = aging.var_names[aging.var['highly_variable'].to_numpy()].astype(str).tolist()
    if len(genes) != 2000 or genes != sorted(genes):
        raise RuntimeError(f'Expected exactly 2000 sorted aging HVGs, found {len(genes)}')
    atomic_text(OUT / 'genes.txt', '\n'.join(genes) + '\n')
    record['aging_genes_sha256'] = sha(OUT / 'genes.txt')
    record['donor_order'] = list(DONORS)
    samples = aging.obs['sample_id'].astype(str).to_numpy()
    for index, donor in enumerate(DONORS):
        donor_indices = np.flatnonzero(samples == donor)
        if len(donor_indices) < 4000:
            raise RuntimeError(f'Donor {donor} has fewer than 4000 cells')
        local = np.random.default_rng(20261500 + index).choice(len(donor_indices), 4000, replace=False)
        selected = donor_indices[local]
        a = aging[selected, genes]
        log_x = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
        log_x = np.asarray(log_x, dtype=np.float32)
        pca = PCA(n_components=30, svd_solver='full', random_state=20261014)
        z = np.asarray(pca.fit_transform(log_x)[:, :D], dtype=np.float32)
        loading = np.array(pca.components_[:D])
        labels = np.asarray([coarse_aging(t) for t in a.obs['cell_type'].astype(str)], dtype=str)
        bars = a.obs_names.astype(str).to_numpy(dtype=str)
        atomic_npz(copy_path(donor), z=z, W=loading, logX=log_x, labels=labels,
                   barcodes=bars, selected_global_indices=selected, selected_within_donor_indices=local)
        record['copies'][donor] = {'n_source_cells': len(donor_indices), 'n_cells': 4000,
                                   'selection_seed': 20261500 + index,
                                   'selected_global_indices': selected.tolist(),
                                   'selected_within_donor_indices': local.tolist(), 'barcodes': bars.tolist(),
                                   'z_sha256': array_sha(z), 'W_sha256': array_sha(loading),
                                   'logX_sha256': array_sha(log_x),
                                   'label_counts': {str(k): int(v) for k, v in pd.Series(labels).value_counts().items()}}
        log('donor_prepared', donor=donor, n_cells=4000)
    for pair, copy1, copy2 in PAIR_SPECS:
        a, b = read_copy(copy1), read_copy(copy2)
        u, _, vt = np.linalg.svd(b['W'] @ a['W'].T, full_matrices=False)
        qgene = np.asarray(u @ vt, dtype=np.float64)
        atomic_npz(pair_input_path(pair), Q_gene=qgene)
        record['pairs'][pair] = {'copy1': copy1, 'copy2': copy2,
                                 'Q_gene': qgene.tolist(), 'Q_gene_sha256': array_sha(qgene)}
    atomic_json(path, record)
    log('check_complete', resumed=False, n_copies=10, n_pairs=5, n_aging_genes=len(genes))


def model_path(copy, seed):
    return CACHE / 'drifts' / f'{copy}_{int(seed)}.pt'


def fit_row(copy, seed, saved):
    return {'fit_id': f'{copy}_{int(seed)}', 'copy': copy, 'd': D, 'seed': int(seed),
            'median_asym': float(saved['symmetry']['median_asym']),
            'max_asym': float(saved['symmetry']['max_asym'])}


def fit_worker(item):
    # Primitive copy metadata only: no torch shared-memory transport.
    torch.set_num_threads(1)
    copy, seed = item
    path = model_path(copy, seed)
    if path.exists():
        return fit_row(copy, seed, torch.load(path, map_location='cpu', weights_only=False))
    z = read_copy(copy)['z']
    drift = EnergyDriftAdapter(D, config=cfg(seed)).fit(np.asarray(z, dtype=np.float32), verbose=True)
    symmetry = drift.symmetry_check(np.asarray(z, dtype=np.float32), n_cells=100)
    saved = {'d': D, 'copy': copy, 'state_dict': drift._model.state_dict(), 'config': asdict(drift.config),
             'symmetry': {k: float(v) for k, v in symmetry.items()}, 'input_sha256': array_sha(z)}
    atomic_model(path, saved)
    return fit_row(copy, seed, saved)


def write_symmetry(rows):
    frame = pd.DataFrame(rows, columns=('fit_id', 'copy', 'd', 'seed', 'median_asym', 'max_asym')).sort_values('fit_id')
    if frame['fit_id'].duplicated().any():
        raise RuntimeError('Duplicate fit ID')
    atomic_csv(OUT / 'fit_symmetry.csv', frame)
    return frame


def fits_stage():
    make_meta()
    rows, todo = [], []
    for copy in ('K2_A', 'K2_B', *DONORS):
        input_hash = array_sha(read_copy(copy)['z'])
        for seed in SEEDS:
            path = model_path(copy, seed)
            if path.exists():
                saved = torch.load(path, map_location='cpu', weights_only=False)
                if saved['d'] != D or saved['input_sha256'] != input_hash or saved['copy'] != copy:
                    raise RuntimeError(f'Cached fit mismatch: {path}')
                rows.append(fit_row(copy, seed, saved))
            else:
                todo.append((copy, seed))
    if todo:
        with ProcessPoolExecutor(max_workers=3, mp_context=multiprocessing.get_context('spawn')) as executor:
            futures = [executor.submit(fit_worker, item) for item in todo]
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                write_symmetry(rows)
                log('fit_finished', fit_id=row['fit_id'], completed=len(rows), median_asym=row['median_asym'])
    frame = write_symmetry(rows)
    if len(frame) != 50:
        raise RuntimeError(f'O1o requires exactly 50 issued fits; found {len(frame)}')
    value = float(frame.median_asym.max())
    atomic_json(CACHE / 'gated_fits.json', {'symmetry_max': value})
    log('fits_complete', completed=50, new_fits=len(todo), symmetry_max=value)


def get_models(copy):
    models = []
    for seed in SEEDS:
        path = model_path(copy, seed)
        if not path.is_file():
            raise RuntimeError(f'Missing O1o fit; estimation may not create one: {path}')
        saved = torch.load(path, map_location='cpu', weights_only=False)
        if int(saved['d']) != D:
            raise RuntimeError('Cached drift has the wrong dimension')
        models.append(build_model(saved, D))
    return models


def estimate_rotation(pair, copy1, copy2):
    path = CACHE / 'rotations' / f'{pair}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            return {key: np.array(data[key]) for key in data.files}
    torch.set_num_threads(MAIN_TORCH_THREADS)
    a, b = read_copy(copy1), read_copy(copy2)
    mx, my = get_models(copy1), get_models(copy2)
    # Same Test F H subsets as O1n: r=0 in every real pair.
    problem = Problem(mx, my, a['z'], b['z'], 0, D)
    ux, uy, qskew, ex, ey, m3x, m3y = axes_and_skew(a['z'], b['z'])
    qstarts = np.stack([sign_start(index, ux, uy, D) for index in range(1 << D)])
    losses = screen_losses(problem, qstarts)
    if not np.isfinite(losses).all():
        raise RuntimeError(f'Nonfinite E5s screen: {pair}')
    order = np.asarray(sorted(range(1 << D), key=lambda index: (float(losses[index]), index)))
    selected = order[:8]
    torch.set_num_threads(1)
    with ThreadPoolExecutor(max_workers=3) as executor:
        answers = list(executor.map(lambda item: optimise_start(problem, item), [(int(i), qstarts[i]) for i in selected]))
    best = min(answers, key=lambda row: (row['L_hat'], row['start']))
    qhat = np.asarray(best['Q'], dtype=np.float64)
    rank = int(np.flatnonzero(order == best['start'])[0]) + 1
    # The loading-defined rotation is only loaded after both estimates finish.
    with np.load(pair_input_path(pair), allow_pickle=False) as data:
        qgene = np.array(data['Q_gene'])
    lgene = float(problem.loss(torch.as_tensor(qgene, dtype=torch.float64), False)[0])
    lskew = float(problem.loss(torch.as_tensor(qskew, dtype=torch.float64), False)[0])
    lnoa = float(problem.loss(torch.eye(D, dtype=torch.float64), False)[0])
    result = {'Q_hat': qhat, 'Q_skew': qskew, 'Q_gene': qgene, 'screen_L': losses, 'screen_order': order,
              'screen_rank_best': np.asarray(rank), 'L_hat': np.asarray(best['L_hat']),
              'L_gene': np.asarray(lgene), 'L_skew': np.asarray(lskew), 'L_noA': np.asarray(lnoa),
              'eig_X': ex, 'eig_Y': ey, 'm3_X': m3x, 'm3_Y': m3y}
    atomic_npz(path, **result)
    log('rotation_finished', pair=pair, L_hat=float(best['L_hat']), L_gene=lgene, screen_rank_best=rank)
    return result


def harmony_embedding(pair, copy1, copy2):
    path = CACHE / 'harmony' / f'{pair}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            return np.array(data['embedding'])
    import harmonypy
    a, b = read_copy(copy1, include_log=True), read_copy(copy2, include_log=True)
    log_x = np.vstack([a['logX'], b['logX']])
    pca = PCA(n_components=10, svd_solver='full', random_state=20261014)
    joint = pca.fit_transform(log_x)
    metadata = pd.DataFrame({'copy': ['A'] * len(a['z']) + ['B'] * len(b['z'])})
    harmony = harmonypy.run_harmony(joint, metadata, 'copy', random_state=20261014, max_iter_harmony=20)
    embedding = np.asarray(harmony.Z_corr.T, dtype=np.float64)
    if embedding.shape != (len(metadata), D) or not np.isfinite(embedding).all():
        raise RuntimeError(f'Unexpected/nonfinite Harmony embedding: {pair}, {embedding.shape}')
    atomic_npz(path, embedding=embedding)
    log('harmony_finished', pair=pair)
    return embedding


def pair_checkpoint_path(pair):
    return CACHE / 'pairs' / f'{pair}.npz'


def pair_complete(pair):
    path = pair_checkpoint_path(pair)
    if not path.exists():
        return False
    with np.load(path, allow_pickle=False) as data:
        rows = json.loads(str(data['rows_json'].item()))
        return len(rows) == 5 and {row['method'] for row in rows} == set(METHODS) and data['Q_hat'].shape == (D, D)


def score_pair(pair, copy1, copy2):
    if pair_complete(pair):
        return
    rotation = estimate_rotation(pair, copy1, copy2)
    a, b = read_copy(copy1), read_copy(copy2)
    x, y = a['z'], b['z']
    datasets = np.asarray(['A'] * len(x) + ['B'] * len(y))
    labels1, labels2 = a['labels'], b['labels']
    types = np.concatenate([labels1, labels2])
    used = sorted(t for t in set(labels1) & set(labels2) if t != 'RBC' and np.sum(labels1 == t) >= 20 and np.sum(labels2 == t) >= 20)
    mask = np.isin(types, used)
    if len(used) < 2:
        raise RuntimeError(f'Fewer than two eligible silhouette coarse types: {pair}')
    qhat, qskew, qgene = rotation['Q_hat'], rotation['Q_skew'], rotation['Q_gene']
    rows = []
    embeddings = {}
    for method in METHODS:
        if method == 'harmony':
            embedding = harmony_embedding(pair, copy1, copy2)
            loss = None
        else:
            q = {'noA': np.eye(D), 'skew': qskew, 'E5s': qhat, 'gene': qgene}[method]
            embedding = np.vstack([x @ q.T, y])
            loss = float(rotation[{'noA': 'L_noA', 'skew': 'L_skew', 'E5s': 'L_hat', 'gene': 'L_gene'}[method]])
        mix = mixing(embedding, datasets)
        sil = float(silhouette_score(embedding[mask], types[mask], metric='euclidean'))
        if not np.isfinite([mix, sil]).all():
            raise RuntimeError(f'Nonfinite functional scores: {pair}, {method}')
        err = None
        if method in ('skew', 'E5s'):
            q = qskew if method == 'skew' else qhat
            err = float(np.linalg.norm(q - qgene, ord='fro') / np.sqrt(20.0))
        rows.append({'pair': pair, 'method': method, 'mix': mix, 'sil': sil, 'err_gene': err,
                     'L_hat': loss, 'L_gene': float(rotation['L_gene']),
                     'screen_rank_best': int(rotation['screen_rank_best']) if method == 'E5s' else None})
        embeddings[f'embedding_{method}'] = embedding
    atomic_npz(pair_checkpoint_path(pair), rows_json=np.asarray(json.dumps(rows, sort_keys=True, allow_nan=False)),
               Q_hat=qhat, Q_skew=qskew, Q_gene=qgene, datasets=datasets, types=types,
               sil_types=np.asarray(used), **embeddings)
    log('pair_finished', pair=pair, rows=rows)


def write_pairs_csv():
    rows = []
    for pair, _, _ in PAIR_SPECS:
        if pair_complete(pair):
            with np.load(pair_checkpoint_path(pair), allow_pickle=False) as data:
                rows.extend(json.loads(str(data['rows_json'].item())))
    frame = pd.DataFrame(rows, columns=PAIR_COLUMNS)
    atomic_csv(OUT / 'pairs.csv', frame)
    return frame


def write_q_hats():
    arrays = {'pairs': np.asarray([pair for pair, _, _ in PAIR_SPECS]), 'd': np.asarray(D)}
    for pair, _, _ in PAIR_SPECS:
        if pair_complete(pair):
            with np.load(pair_checkpoint_path(pair), allow_pickle=False) as data:
                for key in ('Q_hat', 'Q_skew', 'Q_gene'):
                    arrays[f'{pair}_{key}'] = np.array(data[key])
    atomic_npz(ARTIFACTS / 'Q_hats.npz', **arrays)


def functional_gates(frame):
    if len(frame) != 25 or frame.duplicated(['pair', 'method']).any():
        raise RuntimeError('Expected exactly 5 pairs x 5 methods')
    gates = {}
    for dataset in ('K2', 'A'):
        subset = frame.loc[frame.pair == 'K2'] if dataset == 'K2' else frame.loc[frame.pair != 'K2']
        for method in METHODS:
            gates[f'mix_{method}_{dataset}'] = float(subset.loc[subset.method == method, 'mix'].median())
        for method in ('noA', 'E5s'):
            gates[f'sil_{method}_{dataset}'] = float(subset.loc[subset.method == method, 'sil'].median())
        gates[f'err_gene_{dataset}'] = float(subset.loc[subset.method == 'E5s', 'err_gene'].median())
    return gates


def estimate_stage():
    make_meta()
    for pair, copy1, copy2 in PAIR_SPECS:
        score_pair(pair, copy1, copy2)
        write_pairs_csv()
        write_q_hats()
    frame = write_pairs_csv()
    gates = functional_gates(frame)
    atomic_json(CACHE / 'gated_estimate.json', gates)
    log('estimate_complete', n_pairs=5, **gates)


def pair_passes(frame, pair):
    rows = frame.loc[frame.pair == pair].set_index('method')
    return bool(rows.loc['E5s', 'mix'] >= 0.80 and rows.loc['E5s', 'sil'] >= 0.9 * rows.loc['noA', 'sil'] and
                rows.loc['E5s', 'mix'] >= rows.loc['skew', 'mix'] + 0.10)


def gated_values(frame):
    symmetry = pd.read_csv(OUT / 'fit_symmetry.csv')
    if len(symmetry) != 50:
        raise RuntimeError('Exactly 50 fits required before report')
    gates = {'symmetry_max': float(symmetry.median_asym.max()), **functional_gates(frame)}
    k2_pass = pair_passes(frame, 'K2')
    n_pass = sum(int(pair_passes(frame, pair)) for pair, _, _ in PAIR_SPECS if pair != 'K2')
    gates['n_pass_A'] = int(n_pass)
    if gates['symmetry_max'] > 1e-4:
        gates['o1o_outcome'] = 'REVISE'
    elif k2_pass:
        gates['o1o_outcome'] = 'works_real' if n_pass >= 3 else 'chemistry_only'
    else:
        gates['o1o_outcome'] = 'donor_only' if n_pass >= 3 else 'fails'
    if set(gates) != set(GATED_KEYS) or len(gates) != 19:
        raise RuntimeError('O1o requires exactly 19 gated keys')
    return {key: gates[key] for key in GATED_KEYS}


def mix_figure(frame):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    positions = np.arange(len(PAIR_SPECS))
    width = 0.15
    for i, method in enumerate(METHODS):
        subset = frame.loc[frame.method == method].set_index('pair')
        for axis, metric in zip(axes, ('mix', 'sil')):
            values = [float(subset.loc[pair, metric]) for pair, _, _ in PAIR_SPECS]
            axis.bar(positions + (i - 2) * width, values, width, label=method)
    axes[0].axhline(0.80, linestyle='--', color='black', linewidth=1, label='mix threshold 0.80')
    axes[0].set_ylabel('mixing (opposite-copy neighbor fraction / 0.5)')
    axes[1].set_ylabel('coarse-type silhouette')
    for axis in axes:
        axis.set_xticks(positions, [pair for pair, _, _ in PAIR_SPECS])
        axis.grid(axis='y', alpha=0.2)
        axis.set_axisbelow(True)
        axis.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / 'fig_mix.png', dpi=180)
    plt.close(fig)


def umap_coordinates(pair, method, embedding):
    path = CACHE / 'umap' / f'{pair}_{method}.npz'
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            if str(data['input_sha256'].item()) != array_sha(embedding):
                raise RuntimeError('Cached UMAP input fingerprint changed')
            return np.array(data['xy'])
    import anndata as ad
    import scanpy as sc
    # O1g settings exactly; embedding coordinates stay separate from estimation.
    a = ad.AnnData(X=np.zeros((len(embedding), 1), dtype=np.float32))
    a.obsm['view'] = np.asarray(embedding, dtype=np.float64)
    sc.pp.neighbors(a, n_neighbors=15, use_rep='view', random_state=20261014)
    sc.tl.umap(a, min_dist=0.5, random_state=20261014)
    xy = np.asarray(a.obsm['X_umap'], dtype=np.float64)
    atomic_npz(path, xy=xy, input_sha256=np.asarray(array_sha(embedding)))
    return xy


def umap_figure(pair, filename):
    with np.load(pair_checkpoint_path(pair), allow_pickle=False) as data:
        datasets, types = np.array(data['datasets']), np.array(data['types'])
        rows = {row['method']: row for row in json.loads(str(data['rows_json'].item()))}
        coordinates = {method: umap_coordinates(pair, method, np.array(data[f'embedding_{method}']))
                       for method in ('noA', 'E5s', 'gene', 'harmony')}
    copy1, copy2 = next((a, b) for name, a, b in PAIR_SPECS if name == pair)
    fig, axes = plt.subplots(4, 2, figsize=(11, 17))
    palette = plt.get_cmap('tab10')
    for i, method in enumerate(('noA', 'E5s', 'gene', 'harmony')):
        xy, row = coordinates[method], rows[method]
        for copy, name, color in (('A', copy1, 'C0'), ('B', copy2, 'C1')):
            mask = datasets == copy
            axes[i, 0].scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=0.7, color=color, label=name, linewidths=0)
        for j, label in enumerate(sorted(set(types))):
            mask = types == label
            axes[i, 1].scatter(xy[mask, 0], xy[mask, 1], s=4, alpha=0.7, color=palette(j % 10), label=label, linewidths=0)
        for j, suffix in enumerate(('copy', 'coarse type')):
            axis = axes[i, j]
            axis.set_title(f"{pair} / {method}: {suffix}; mix={row['mix']:.3f}, sil={row['sil']:.3f}")
            axis.set_xticks([])
            axis.set_yticks([])
            axis.legend(markerscale=2, fontsize=7, ncol=1 if j == 0 else 2)
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


def regeneration_commands():
    command = 'cd "$HOME/Documents/scopalign-manus-O1o" && .venv/bin/python reports/O1o/scripts/o1o_run.py'
    result = {'artifacts/inputs.json': command + ' --stage check', 'genes.txt': command + ' --stage check',
              'fit_symmetry.csv': command + ' --stage fits', 'pairs.csv': command + ' --stage estimate',
              'artifacts/Q_hats.npz': command + ' --stage estimate',
              'fig_mix.png': command + ' --stage report', 'fig_umap_K2.png': command + ' --stage report',
              'fig_umap_F30F36.png': command + ' --stage report', 'report.json': command + ' --stage report',
              'scripts/o1o_run.py': command + ' --stage all',
              'requirements.lock': 'cd "$HOME/Documents/scopalign-manus-O1o" && .venv/bin/python -m pip freeze > reports/O1o/requirements.lock'}
    for key in GATED_KEYS:
        stage = 'fits' if key == 'symmetry_max' else 'report' if key in ('n_pass_A', 'o1o_outcome') else 'estimate'
        result[key] = command + f' --stage {stage}'
    return result


def report_stage():
    meta = make_meta()
    for row in json.loads((ARTIFACTS / 'inputs.json').read_text())['readonly_files']:
        if sha(row['path']) != row['sha256']:
            raise RuntimeError(f'Read-only source changed during run: {row["path"]}')
    frame = write_pairs_csv()
    gates = gated_values(frame)
    mix_figure(frame)
    umap_figure('K2', 'fig_umap_K2.png')
    umap_figure('F30F36', 'fig_umap_F30F36.png')
    with (OUT / 'requirements.lock').open('w') as handle:
        subprocess.run([sys.executable, '-m', 'pip', 'freeze'], check=True, stdout=handle)
    deliverables = ('requirements.lock', 'genes.txt', 'pairs.csv', 'fit_symmetry.csv',
                    'fig_mix.png', 'fig_umap_K2.png', 'fig_umap_F30F36.png', 'scripts/o1o_run.py')
    artifacts = ('artifacts/Q_hats.npz', 'artifacts/inputs.json')
    entries = lambda paths: [{'path': path, 'sha256': sha(OUT / path)} for path in paths]
    report = {'task_id': 'O1o', 'executor': 'manus', 'branch': 'task/O1o',
              'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'seed': 20261014, 'started': meta['started'], 'finished': datetime.now(timezone.utc).isoformat(),
              'status': 'complete', 'gated': gates, 'deliverables': entries(deliverables), 'artifacts': entries(artifacts),
              'regeneration': regeneration_commands(), 'deviations': DEVIATIONS, 'deviations_requested': [],
              'observations': [], 'known_failures': KNOWN_FAILURES, 'lockfile_sha256': sha(OUT / 'requirements.lock')}
    atomic_json(OUT / 'report.json', report)
    log('report_complete', **gates)
    subprocess.run([sys.executable, str(ROOT / 'checker/check_report.py'), str(ROOT / 'tasks/O1o.md'), str(OUT / 'report.json')], cwd=ROOT, check=True)


def checkpoint_stage(stage):
    subprocess.run(['git', 'add', '--', 'reports/O1o'], cwd=ROOT, check=True)
    changed = subprocess.run(['git', 'diff', '--cached', '--quiet'], cwd=ROOT).returncode
    if changed == 1:
        staged = subprocess.check_output(['git', 'diff', '--cached', '--name-only'], cwd=ROOT, text=True).splitlines()
        if any(not name.startswith('reports/O1o/') for name in staged):
            raise RuntimeError('Refusing to commit paths outside reports/O1o')
        subprocess.run(['git', 'commit', '-m', f'task: O1o {stage} stage', '-m',
                        'Co-authored-by: Manus <dev-agent@manus.ai>'], cwd=ROOT, check=True)
    elif changed != 0:
        raise RuntimeError('Unable to inspect stage checkpoint')
    subprocess.run(['git', 'push', 'origin', 'task/O1o'], cwd=ROOT, check=True)
    log('stage_checkpoint', stage=stage, commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True, choices=('check', 'fits', 'estimate', 'report', 'all'))
    stage = parser.parse_args().stage
    make_meta()
    with (OUT / 'run.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another O1o runner owns the output lock')
        for name, function in (('check', check_stage), ('fits', fits_stage), ('estimate', estimate_stage), ('report', report_stage)):
            if stage in (name, 'all'):
                log('stage_started', stage=name)
                function()
                if stage == 'all':
                    checkpoint_stage(name)
        log('run_complete', stage=stage)


if __name__ == '__main__':
    main()
