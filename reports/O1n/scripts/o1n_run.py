#!/usr/bin/env python
"""O1n: pre-registered 2/3/5/10-PC rotation sweep against the skew baseline.

All source datasets and O1/O1b/O1g/O1l references are read-only. Only the
160 issued low-dimensional models are fitted. Run --stage all detached;
each fit and case is atomically checkpointed, and each stage is committed
and pushed. Runtime caches and prepared cell arrays stay local and ignored.
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
OUT = ROOT / 'reports' / 'O1n'
CACHE = OUT / 'cache'
ARTIFACTS = OUT / 'artifacts'
# Configure caches before any matplotlib/scanpy import, including worker imports.
for key, name in (('MPLCONFIGDIR', '.mplconfig'), ('XDG_CACHE_HOME', '.cache'),
                  ('NUMBA_CACHE_DIR', '.numba')):
    target = OUT / name
    target.mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(target)
sys.dont_write_bytecode = True

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.stats import gaussian_kde, ortho_group, spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP

BASE = Path.home() / 'Documents' / 'scopalign-manus'
K0 = BASE / 'reports' / 'K0'
LAT1, LAT2 = K0 / 'latent_pca30_h1.npy', K0 / 'latent_pca30_h2.npy'
GENES = K0 / 'genes.txt'
SOURCE = BASE / 'data' / 'pbmc_10k_v3' / 'pbmc_10k_v3.h5ad'
HALVES = BASE / 'data' / 'pbmc_10k_v3' / 'k0_halves.csv'
REFERENCE_SCRIPTS = [
    BASE / 'reports' / 'O1' / 'scripts' / 'o1_run.py',
    BASE.parent / 'scopalign-manus-O1b' / 'reports' / 'O1b' / 'scripts' / 'o1b_run.py',
    BASE.parent / 'scopalign-manus-O1g' / 'reports' / 'O1g' / 'scripts' / 'o1g_run.py',
    BASE.parent / 'scopalign-manus-O1l' / 'reports' / 'O1l' / 'scripts' / 'o1l_run.py',
]
DS = (2, 3, 5, 10)
REPS = (0, 1, 2)
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
N_R = 2759
MAIN_TORCH_THREADS = torch.get_num_threads()
KNOWN_FAILURES = [
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]',
    'tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility',
]
GATED_KEYS = ('symmetry_max', 'latent_check_rel') + tuple(
    f'{prefix}_d{d}' for prefix in
    ('err_R', 'err_R_skew', 'mix_F', 'mix_F_skew', 'mix_F_noA') for d in DS
) + ('rho_kde_d2', 'rho_kde_d3', 'd_max_op', 'o1n_outcome')
DEVIATIONS = [{
    'what': 'Used report-local Matplotlib/fontconfig/Numba caches and O1l\'s in-process refinement executor, capped at three workers.',
    'why': 'Shared environment and input worktrees are read-only. O1l\'s issued runner uses an in-process thread executor to avoid the local PyTorch spawned shared-memory restriction; all loss and LBFGS settings are unchanged.',
}]


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


def make_meta():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(exist_ok=True)
    ARTIFACTS.mkdir(exist_ok=True)
    (OUT / '.gitignore').write_text('cache/\n__pycache__/\nrun.log\nrun.pid\nrun.lock\n.mplconfig/\n.cache/\n.numba/\n')
    path = OUT / 'run_meta.json'
    if not path.exists():
        atomic_json(path, {'started': datetime.now(timezone.utc).isoformat()})
    return json.loads(path.read_text())


def case_specs():
    # Fixed outer d, then Test R r=0,1,2, then Test F.
    for d in DS:
        for r in REPS:
            yield f'R_d{d}_r{r}', 'R', d, r
        yield f'F_d{d}', 'F', d, 0


def input_path(case):
    return CACHE / 'inputs' / f'{case}.npz'


def read_case(case):
    with np.load(input_path(case), allow_pickle=False) as data:
        return {key: np.array(data[key]) for key in data.files}


def load_halves():
    # O1b load_halves exactly: QC filter, original order, and K0 membership.
    import scanpy as sc
    a = sc.read_h5ad(SOURCE)
    a = a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()
    tab = pd.read_csv(HALVES).set_index('barcode')['half']
    names = a.obs_names.astype(str)
    m1 = names.isin(tab.index) & (names.map(tab) == 'H1')
    m2 = names.isin(tab.index) & (names.map(tab) == 'H2')
    return a[m1].copy(), a[m2].copy()


def dense_normalised(a):
    import scanpy as sc
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)


def check_stage():
    make_meta()
    readonly_paths = [LAT1, LAT2, GENES, SOURCE, HALVES, *REFERENCE_SCRIPTS]
    fingerprints = [{'path': str(p), 'sha256': sha(p)} for p in readonly_paths]
    previous_path = ARTIFACTS / 'inputs.json'
    if previous_path.exists():
        previous = json.loads(previous_path.read_text())
        if previous['readonly_files'] != fingerprints:
            raise RuntimeError('A prescribed read-only input changed since the initial check')
        for case, _, _, _ in case_specs():
            if not input_path(case).exists():
                raise RuntimeError(f'Missing prepared case after check: {case}')
        log('check_complete', latent_check_rel=previous['latent_check_rel'], resumed=True)
        return
    z1 = np.load(LAT1).astype(np.float32)
    z2 = np.load(LAT2).astype(np.float32)
    if z1.shape != (5518, 30) or z2.shape[1] != 30:
        raise RuntimeError(f'Unexpected K0 latent shapes: {z1.shape}, {z2.shape}')
    h1, h2 = load_halves()
    if len(h1) != len(z1) or len(h2) != len(z2):
        raise RuntimeError('Source half order/size does not match K0 latents')
    genes = [s.strip() for s in GENES.read_text().splitlines() if s.strip()]
    if len(genes) != 2000:
        raise RuntimeError('K0 genes list must contain exactly 2000 genes')
    derived, weights, rel_fro, rel_max, max_abs = [], [], [], [], []
    for half, cached in ((h1, z1), (h2, z2)):
        pca = PCA(n_components=30, svd_solver='full', random_state=20261014)
        latent = pca.fit_transform(dense_normalised(half[:, genes]))
        delta = latent.astype(np.float64) - cached.astype(np.float64)
        rel_fro.append(float(np.linalg.norm(delta) / np.linalg.norm(cached.astype(np.float64))))
        rel_max.append(float(np.max(np.abs(delta)) / np.max(np.abs(cached))))
        max_abs.append(float(np.max(np.abs(delta))))
        derived.append(latent)
        weights.append(pca.components_)
    # Relative counterpart of O1's maximum absolute-difference check: infinity
    # norm of latent difference divided by infinity norm of that half's cache.
    latent_check_rel = float(max(rel_max))
    labels1 = h1.obs['cell_type'].astype(str).to_numpy(dtype=str)
    labels2 = h2.obs['cell_type'].astype(str).to_numpy(dtype=str)
    atomic_npz(CACHE / 'prepared.npz', W1=weights[0], W2=weights[1], labels1=labels1, labels2=labels2)
    cases = []
    targets = {}
    for case, test, d, r in case_specs():
        if test == 'R':
            rng = np.random.default_rng(20261030 + r)
            perm = rng.permutation(5518)
            sx, sy = np.sort(perm[:N_R]), np.sort(perm[N_R:2 * N_R])
            target = np.asarray(ortho_group.rvs(dim=d, random_state=20261031 + r), dtype=np.float64)
            x = z1[sx, :d]
            y = (z1[sy, :d].astype(np.float64) @ target.T).astype(np.float32)
            arrays = {'X': x, 'Y': y, 'Q_target': target, 'S_X': sx, 'S_Y': sy}
        else:
            u, _, vt = np.linalg.svd(weights[1][:d] @ weights[0][:d].T, full_matrices=False)
            target = np.asarray(u @ vt, dtype=np.float64)
            x, y = z1[:, :d], z2[:, :d]
            arrays = {'X': x, 'Y': y, 'Q_target': target}
        atomic_npz(input_path(case), **arrays)
        targets[case] = target.tolist()
        cases.append({'case': case, 'test': test, 'd': d, 'r': r, 'n_X': len(x), 'n_Y': len(y),
                      'X_sha256': array_sha(x), 'Y_sha256': array_sha(y), 'Q_target_sha256': array_sha(target)})
    record = {'readonly_files': fingerprints, 'latent_shapes': [list(z1.shape), list(z2.shape)],
              'latent_check_rel': latent_check_rel,
              'latent_check_rel_definition': 'max over halves of max(abs(rederived-K0))/max(abs(K0)); PCA-30 re-derived before truncation',
              'relative_frobenius_by_half': rel_fro, 'relative_max_by_half': rel_max,
              'max_abs_by_half': max_abs, 'cases': cases, 'Q_targets': targets,
              'fit_seeds': list(SEEDS), 'dimensions': list(DS),
              'H_subset_seeds': {case: 20261040 + r for case, _, _, r in case_specs()}}
    atomic_json(ARTIFACTS / 'inputs.json', record)
    log('check_complete', latent_check_rel=latent_check_rel, resumed=False)


def cfg(seed):
    # O1/O1b config unchanged except EnergyDriftAdapter's input dimension.
    return EnergyDriftConfig(hidden=96, depth=3, n_epochs=1500, seed=int(seed),
                             batch_size=256, loss='esm', lr=5e-4, beta=1.0,
                             esm_laplacian='hutchinson', device='cpu')


def model_path(case, side, seed):
    return CACHE / 'drifts' / f'{case}_{side}_{int(seed)}.pt'


def fit_row(case, side, seed, saved):
    return {'fit_id': f'{case}_{side}_{seed}', 'case': case, 'd': int(saved['d']),
            'side': side, 'seed': int(seed), 'median_asym': float(saved['symmetry']['median_asym']),
            'max_asym': float(saved['symmetry']['max_asym'])}


def fit_worker(item):
    # Only primitive metadata is passed to this spawn worker: no shared torch
    # storage and no torch model/Problem is transported between processes.
    torch.set_num_threads(1)
    case, d, side, seed = item
    path = model_path(case, side, seed)
    if path.exists():
        return fit_row(case, side, seed, torch.load(path, map_location='cpu', weights_only=False))
    z = read_case(case)[side]
    drift = EnergyDriftAdapter(d, config=cfg(seed)).fit(np.asarray(z, dtype=np.float32), verbose=True)
    symmetry = drift.symmetry_check(np.asarray(z, dtype=np.float32), n_cells=100)
    saved = {'d': int(d), 'state_dict': drift._model.state_dict(), 'config': asdict(drift.config),
             'symmetry': {k: float(v) for k, v in symmetry.items()}, 'input_sha256': array_sha(z)}
    atomic_model(path, saved)
    return fit_row(case, side, seed, saved)


def write_symmetry(rows):
    columns = ('fit_id', 'case', 'd', 'side', 'seed', 'median_asym', 'max_asym')
    frame = pd.DataFrame(rows, columns=columns).sort_values('fit_id', kind='stable')
    if frame['fit_id'].duplicated().any():
        raise RuntimeError('Duplicate fit ID')
    atomic_csv(OUT / 'fit_symmetry.csv', frame)
    return frame


def fits_stage():
    make_meta()
    rows, todo = [], []
    for case, _, d, _ in case_specs():
        for side in ('X', 'Y'):
            for seed in SEEDS:
                path = model_path(case, side, seed)
                if path.exists():
                    saved = torch.load(path, map_location='cpu', weights_only=False)
                    if saved['d'] != d or saved['input_sha256'] != array_sha(read_case(case)[side]):
                        raise RuntimeError(f'Cached low-d fit mismatch: {path}')
                    rows.append(fit_row(case, side, seed, saved))
                else:
                    todo.append((case, d, side, seed))
    if todo:
        with ProcessPoolExecutor(max_workers=3, mp_context=multiprocessing.get_context('spawn')) as executor:
            futures = {executor.submit(fit_worker, item): item for item in todo}
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                write_symmetry(rows)
                log('fit_finished', fit_id=row['fit_id'], completed=len(rows), median_asym=row['median_asym'])
    frame = write_symmetry(rows)
    if len(frame) != 160:
        raise RuntimeError(f'Expected exactly 160 low-dimensional fits, got {len(frame)}')
    value = float(frame['median_asym'].max())
    atomic_json(CACHE / 'gated_fits.json', {'symmetry_max': value})
    log('fits_complete', new_fits=len(todo), completed=160, symmetry_max=value)


def build_model(saved, d):
    model = _EnergyMLP(d, int(saved['config']['hidden']), int(saved['config']['depth']))
    model.load_state_dict(saved['state_dict'])
    model = model.to(dtype=torch.float64, device='cpu').eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model


def get_models(case, d):
    mx, my = [], []
    for seed in SEEDS:
        for side, bucket in (('X', mx), ('Y', my)):
            path = model_path(case, side, seed)
            if not path.is_file():
                raise RuntimeError(f'Missing issued O1n fit: {path}; no fit is created during estimation')
            saved = torch.load(path, map_location='cpu', weights_only=False)
            if int(saved['d']) != d:
                raise RuntimeError('Wrong input dimension in cached fit')
            bucket.append(build_model(saved, d))
    return mx, my


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


def estimate(case, d, r):
    torch.set_num_threads(MAIN_TORCH_THREADS)
    data = read_case(case)
    x, y = data['X'], data['Y']
    mx, my = get_models(case, d)
    problem = Problem(mx, my, x, y, r, d)
    ux, uy, qskew, ex, ey, m3x, m3y = axes_and_skew(x, y)
    qstarts = np.stack([sign_start(i, ux, uy, d) for i in range(1 << d)])
    losses = screen_losses(problem, qstarts)
    if not np.isfinite(losses).all():
        raise RuntimeError(f'Nonfinite screen loss: {case}')
    order = np.asarray(sorted(range(1 << d), key=lambda i: (float(losses[i]), i)))
    selected = order[:min(8, 1 << d)]
    torch.set_num_threads(1)
    with ThreadPoolExecutor(max_workers=3) as executor:
        answers = list(executor.map(lambda item: optimise_start(problem, item),
                                    [(int(i), qstarts[i]) for i in selected]))
    best = min(answers, key=lambda row: (row['L_hat'], row['start']))
    qhat = best['Q']
    # Known targets enter evaluation only after screening, optimization, and
    # positions-only estimation have completed.
    target = data['Q_target']
    result = {'err': float(np.linalg.norm(qhat - target, 'fro') / np.sqrt(2 * d)),
              'err_skew': float(np.linalg.norm(qskew - target, 'fro') / np.sqrt(2 * d)),
              'L_hat': best['L_hat'], 'screen_rank_best': int(np.flatnonzero(order == best['start'])[0]) + 1}
    if case.startswith('R'):
        result['L_star'] = float(problem.loss(torch.as_tensor(target, dtype=torch.float64), False)[0])
    return result, {'Q_hat': qhat, 'Q_skew': qskew, 'Q_target': target,
                    'screen_L': losses, 'screen_order': order, 'eig_X': ex, 'eig_Y': ey,
                    'm3_X': m3x, 'm3_Y': m3y}


def checkpoint_path(case):
    return CACHE / 'cases' / f'{case}.npz'


def case_complete(case):
    path = checkpoint_path(case)
    if not path.exists():
        return False
    with np.load(path, allow_pickle=False) as data:
        return {'result_json', 'Q_hat', 'Q_skew', 'Q_target', 'screen_L'}.issubset(data.files)


def read_result(case):
    with np.load(checkpoint_path(case), allow_pickle=False) as data:
        return json.loads(str(data['result_json'].item()))


def save_result(case, result, arrays):
    atomic_npz(checkpoint_path(case), result_json=np.asarray(json.dumps(result, sort_keys=True, allow_nan=False)), **arrays)


def write_results(test):
    rows = []
    for case, t, _, _ in case_specs():
        if t == test and case_complete(case):
            rows.append(read_result(case))
    if test == 'R':
        columns = ('d', 'r', 'err', 'err_skew', 'L_hat', 'L_star', 'screen_rank_best')
        frame = pd.DataFrame(rows, columns=columns).sort_values(['d', 'r'])
    else:
        columns = ('d', 'err', 'err_skew', 'L_hat', 'screen_rank_best') + tuple(
            f'{metric}_{view}' for view in ('noA', 'skew', 'E5s', 'Q_load') for metric in ('mix', 'sil'))
        frame = pd.DataFrame(rows, columns=columns).sort_values('d')
    atomic_csv(OUT / f'test{test}.csv', frame)
    return frame


def write_q_hats():
    arrays = {'dimensions': np.asarray(DS), 'replicates': np.asarray(REPS)}
    for case, _, _, _ in case_specs():
        if case_complete(case):
            with np.load(checkpoint_path(case), allow_pickle=False) as data:
                for key in ('Q_hat', 'Q_skew', 'Q_target'):
                    arrays[f'{case}_{key}'] = np.array(data[key])
    atomic_npz(ARTIFACTS / 'Q_hats.npz', **arrays)


def rho_stage():
    rows = []
    for d in (2, 3):
        case = f'R_d{d}_r0'
        path = CACHE / f'rho_d{d}.json'
        if path.exists():
            rows.extend(json.loads(path.read_text()))
            continue
        x = read_case(case)['X'].astype(np.float64)
        density = gaussian_kde(x.T)(x.T)
        if np.any(density <= 0):
            raise RuntimeError('Default Gaussian KDE returned nonpositive density at a cell')
        negative_log_density = -np.log(density)
        new = []
        for seed in SEEDS:
            saved = torch.load(model_path(case, 'X', seed), map_location='cpu', weights_only=False)
            model = build_model(saved, d)
            with torch.no_grad():
                energy = model(torch.as_tensor(x, dtype=torch.float64)).numpy()
            rho = float(spearmanr(energy, negative_log_density).statistic)
            if not np.isfinite(rho):
                raise RuntimeError('Nonfinite fitted-energy/KDE Spearman correlation')
            new.append({'d': d, 'r': 0, 'side': 'X', 'seed': seed, 'rho': rho})
        atomic_json(path, new)
        rows.extend(new)
    frame = pd.DataFrame(rows).sort_values(['d', 'seed'])
    atomic_csv(OUT / 'rho_kde.csv', frame)
    return frame


def test_r_stage():
    make_meta()
    for case, test, d, r in case_specs():
        if test != 'R':
            continue
        if not case_complete(case):
            result, arrays = estimate(case, d, r)
            result.update({'d': d, 'r': r})
            save_result(case, result, arrays)
            log('case_finished', case=case, **result)
        write_results('R')
        write_q_hats()
    frame = write_results('R')
    if len(frame) != 12:
        raise RuntimeError('Test R requires 12 complete cases')
    rho = rho_stage()
    gates = {f'err_R_d{d}': float(frame.loc[frame.d == d, 'err'].median()) for d in DS}
    gates.update({f'err_R_skew_d{d}': float(frame.loc[frame.d == d, 'err_skew'].median()) for d in DS})
    gates.update({f'rho_kde_d{d}': float(rho.loc[rho.d == d, 'rho'].median()) for d in (2, 3)})
    atomic_json(CACHE / 'gated_testR.json', gates)
    log('testR_complete', **gates)


def mixing(rep, datasets):
    # O1g's exact fraction of the 30 nearest OTHER cells with opposite copy,
    # divided by 0.5. This is not a cross-copy-only kNN search.
    neighbors = NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep, return_distance=False)[:, 1:]
    return float(np.mean(np.mean(datasets[neighbors] != datasets[:, None], axis=1)) / 0.5)


def functional_views(d, qhat, qskew, target):
    data = read_case(f'F_d{d}')
    x, y = data['X'], data['Y']
    with np.load(CACHE / 'prepared.npz', allow_pickle=False) as prepared:
        l1, l2 = prepared['labels1'], prepared['labels2']
    datasets = np.asarray(['H1'] * len(x) + ['H2'] * len(y))
    types = np.concatenate([l1, l2])
    tested = sorted(t for t in set(l1) & set(l2) if np.sum(l1 == t) >= 50 and np.sum(l2 == t) >= 50)
    mask = np.isin(types, tested)
    rows, coordinates = {}, {}
    for name, q in (('noA', np.eye(d)), ('skew', qskew), ('E5s', qhat), ('Q_load', target)):
        rep = np.vstack([x @ q.T, y])
        rows[f'mix_{name}'] = mixing(rep, datasets)
        rows[f'sil_{name}'] = float(silhouette_score(rep[mask], types[mask], metric='euclidean'))
        coordinates[name] = rep
    return rows, coordinates, datasets, types


def test_f_stage():
    make_meta()
    for d in DS:
        case = f'F_d{d}'
        if not case_complete(case):
            result, arrays = estimate(case, d, 0)
            views, _, _, _ = functional_views(d, arrays['Q_hat'], arrays['Q_skew'], arrays['Q_target'])
            result.update({'d': d, **views})
            save_result(case, result, arrays)
            log('case_finished', case=case, **result)
        write_results('F')
        write_q_hats()
    frame = write_results('F')
    if len(frame) != 4:
        raise RuntimeError('Test F requires four complete cases')
    gates = {}
    for d in DS:
        row = frame.loc[frame.d == d].iloc[0]
        gates.update({f'mix_F_d{d}': float(row['mix_E5s']), f'mix_F_skew_d{d}': float(row['mix_skew']),
                      f'mix_F_noA_d{d}': float(row['mix_noA'])})
    atomic_json(CACHE / 'gated_testF.json', gates)
    log('testF_complete', **gates)


def gated_values():
    gates = {}
    gates['symmetry_max'] = float(pd.read_csv(OUT / 'fit_symmetry.csv')['median_asym'].max())
    gates['latent_check_rel'] = float(json.loads((ARTIFACTS / 'inputs.json').read_text())['latent_check_rel'])
    gates.update(json.loads((CACHE / 'gated_testR.json').read_text()))
    gates.update(json.loads((CACHE / 'gated_testF.json').read_text()))
    passes, thresholds = [], []
    for d in DS:
        threshold_ok = gates[f'err_R_d{d}'] <= 0.10 and gates[f'mix_F_d{d}'] >= 0.80
        beats = gates[f'err_R_d{d}'] <= gates[f'err_R_skew_d{d}'] - 0.05 or gates[f'mix_F_d{d}'] >= gates[f'mix_F_skew_d{d}'] + 0.10
        if threshold_ok:
            thresholds.append(d)
        if threshold_ok and beats:
            passes.append(d)
    gates['d_max_op'] = int(max(passes) if passes else 0)
    if gates['symmetry_max'] > 1e-4 or gates['latent_check_rel'] > 1e-4:
        gates['o1n_outcome'] = 'REVISE'
    else:
        gates['o1n_outcome'] = 'op_works_lowd' if passes else 'positions_suffice' if thresholds else 'fails'
    if set(gates) != set(GATED_KEYS) or len(gates) != 26:
        raise RuntimeError('O1n must contain exactly the 26 prescribed gates')
    return {key: gates[key] for key in GATED_KEYS}


def figures(gates):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    axes[0].plot(DS, [gates[f'err_R_d{d}'] for d in DS], marker='o', label='E5s')
    axes[0].plot(DS, [gates[f'err_R_skew_d{d}'] for d in DS], marker='s', label='skew baseline')
    axes[0].axhline(0.10, color='black', linestyle='--', linewidth=1, label='0.10')
    axes[0].set_ylabel('median Test R rotation error')
    axes[1].plot(DS, [gates[f'mix_F_d{d}'] for d in DS], marker='o', label='E5s')
    axes[1].plot(DS, [gates[f'mix_F_skew_d{d}'] for d in DS], marker='s', label='skew baseline')
    axes[1].axhline(0.80, color='black', linestyle='--', linewidth=1, label='0.80')
    axes[1].set_ylabel('Test F mixing')
    for axis in axes:
        axis.set_xlabel('number of PCs (d)')
        axis.set_xticks(DS)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'fig_sweep.png', dpi=180)
    plt.close(fig)
    with np.load(checkpoint_path('F_d2'), allow_pickle=False) as data:
        metrics, coordinates, datasets, types = functional_views(2, data['Q_hat'], data['Q_skew'], data['Q_target'])
    fig, axes = plt.subplots(4, 2, figsize=(11, 18))
    palette = plt.get_cmap('tab10')
    for index, name in enumerate(('noA', 'skew', 'E5s', 'Q_load')):
        rep = coordinates[name]
        copy_axis, type_axis = axes[index]
        for label, color in (('H1', 'C0'), ('H2', 'C1')):
            mask = datasets == label
            copy_axis.scatter(rep[mask, 0], rep[mask, 1], s=4, alpha=0.7, color=color, label=label, linewidths=0)
        for j, label in enumerate(sorted(set(types))):
            mask = types == label
            type_axis.scatter(rep[mask, 0], rep[mask, 1], s=4, alpha=0.7, color=palette(j % 10), label=label, linewidths=0)
        for axis, suffix in ((copy_axis, 'copy'), (type_axis, 'coarse type')):
            axis.set_title(f"{name}: {suffix}; mix={metrics[f'mix_{name}']:.3f}, sil={metrics[f'sil_{name}']:.3f}")
            axis.set_xlabel('PC 1 in target frame')
            axis.set_ylabel('PC 2 in target frame')
            axis.legend(markerscale=2, fontsize=7, ncol=2 if suffix == 'coarse type' else 1)
    fig.tight_layout()
    fig.savefig(OUT / 'fig_views_d2.png', dpi=180)
    plt.close(fig)


def regeneration_commands():
    command = 'cd "$HOME/Documents/scopalign-manus-O1n" && .venv/bin/python reports/O1n/scripts/o1n_run.py'
    result = {'artifacts/inputs.json': command + ' --stage check', 'fit_symmetry.csv': command + ' --stage fits',
              'testR.csv': command + ' --stage testR', 'rho_kde.csv': command + ' --stage testR',
              'artifacts/Q_hats.npz': command + ' --stage testF', 'testF.csv': command + ' --stage testF',
              'fig_sweep.png': command + ' --stage report', 'fig_views_d2.png': command + ' --stage report',
              'report.json': command + ' --stage report', 'scripts/o1n_run.py': command + ' --stage all',
              'requirements.lock': 'cd "$HOME/Documents/scopalign-manus-O1n" && .venv/bin/python -m pip freeze > reports/O1n/requirements.lock'}
    for key in GATED_KEYS:
        stage = ('fits' if key == 'symmetry_max' else 'check' if key == 'latent_check_rel' else
                 'testR' if key.startswith(('err_R', 'rho_kde')) else
                 'testF' if key.startswith('mix_F') else 'report')
        result[key] = command + f' --stage {stage}'
    return result


def report_stage():
    meta = make_meta()
    if len(pd.read_csv(OUT / 'fit_symmetry.csv')) != 160 or len(pd.read_csv(OUT / 'testR.csv')) != 12 or len(pd.read_csv(OUT / 'testF.csv')) != 4 or len(pd.read_csv(OUT / 'rho_kde.csv')) != 10:
        raise RuntimeError('Report requires all 160 fits, 12 Test R cases, 4 Test F cases and 10 KDE correlations')
    inputs = json.loads((ARTIFACTS / 'inputs.json').read_text())
    for row in inputs['readonly_files']:
        if sha(row['path']) != row['sha256']:
            raise RuntimeError(f'Read-only source changed during execution: {row["path"]}')
    gates = gated_values()
    figures(gates)
    with (OUT / 'requirements.lock').open('w') as handle:
        subprocess.run([sys.executable, '-m', 'pip', 'freeze'], check=True, stdout=handle)
    deliverables = ('requirements.lock', 'testR.csv', 'testF.csv', 'rho_kde.csv', 'fit_symmetry.csv',
                    'fig_sweep.png', 'fig_views_d2.png', 'scripts/o1n_run.py')
    artifacts = ('artifacts/Q_hats.npz', 'artifacts/inputs.json')
    entries = lambda paths: [{'path': str(p), 'sha256': sha(OUT / p)} for p in paths]
    report = {'task_id': 'O1n', 'executor': 'manus', 'branch': 'task/O1n',
              'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'seed': 20261014, 'started': meta['started'], 'finished': datetime.now(timezone.utc).isoformat(),
              'status': 'complete', 'gated': gates, 'deliverables': entries(deliverables),
              'artifacts': entries(artifacts), 'regeneration': regeneration_commands(),
              'deviations': DEVIATIONS, 'deviations_requested': [], 'observations': [],
              'known_failures': KNOWN_FAILURES, 'lockfile_sha256': sha(OUT / 'requirements.lock')}
    atomic_json(OUT / 'report.json', report)
    log('report_complete', **gates)
    subprocess.run([sys.executable, str(ROOT / 'checker' / 'check_report.py'),
                    str(ROOT / 'tasks' / 'O1n.md'), str(OUT / 'report.json')], cwd=ROOT, check=True)


def checkpoint_stage(stage):
    # Never stage shared input/environment links or any other task's outputs.
    subprocess.run(['git', 'add', '--', 'reports/O1n'], cwd=ROOT, check=True)
    changed = subprocess.run(['git', 'diff', '--cached', '--quiet'], cwd=ROOT).returncode
    if changed == 1:
        staged = subprocess.check_output(['git', 'diff', '--cached', '--name-only'], cwd=ROOT, text=True).splitlines()
        if any(not p.startswith('reports/O1n/') for p in staged):
            raise RuntimeError('Refusing a stage checkpoint with a path outside reports/O1n')
        subprocess.run(['git', 'commit', '-m', f'task: O1n {stage} stage',
                        '-m', 'Co-authored-by: Manus <dev-agent@manus.ai>'], cwd=ROOT, check=True)
    elif changed != 0:
        raise RuntimeError(f'Cannot inspect staged files after {stage}')
    # Push even if a resumed stage had no new changes, including after a prior
    # failed push. Never merge or update main.
    subprocess.run(['git', 'push', 'origin', 'task/O1n'], cwd=ROOT, check=True)
    log('stage_checkpoint', stage=stage,
        commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True, choices=('check', 'fits', 'testR', 'testF', 'report', 'all'))
    stage = parser.parse_args().stage
    make_meta()
    # Prevent duplicate detached runners from ever writing the same case/fit.
    with (OUT / 'run.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('An O1n runner already owns the output lock')
        stages = (('check', check_stage), ('fits', fits_stage), ('testR', test_r_stage),
                  ('testF', test_f_stage), ('report', report_stage))
        for name, function in stages:
            if stage in (name, 'all'):
                log('stage_started', stage=name)
                function()
                if stage == 'all':
                    checkpoint_stage(name)
        log('run_complete', stage=stage)


if __name__ == '__main__':
    main()
