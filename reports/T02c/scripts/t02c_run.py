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
from scjdo2.compare import centroid_operator
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'T02c'
K2 = ROOT / 'reports' / 'K2'
T02B = ROOT / 'reports' / 'T02b'
T01D = ROOT / 'reports' / 'T01d'
DATA = ROOT / 'data' / 'pbmc_v2_v3'
LAT = {'a': K2 / 'latent_pca30_a.npy', 'b': K2 / 'latent_pca30_b.npy'}
OPS = {'a': K2 / 'ops_a.npz', 'b': K2 / 'ops_b.npz'}
LOADINGS = T02B / 'loadings.npz'
V2 = DATA / 'pbmc_v2.h5ad'
V3_THIN = DATA / 'pbmc_v3_thin.h5ad'
CELLS = DATA / 'k2_cells.csv'
GENES = T01D / 'genes.txt'
INPUTS = (LAT['a'], LAT['b'], OPS['a'], OPS['b'], LOADINGS, CELLS, V2, V3_THIN, GENES)
SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)
D, R = 30, 10
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
def centroid_path(): return OUT / 'centroid_ops.npz'
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

def clusters_and_weights(oa, ob):
    la, lb = oa['cell_type'].astype(str), ob['cell_type'].astype(str)
    clusters = np.asarray(sorted(t for t in set(la) & set(lb) if np.sum(la == t) >= 100 and np.sum(lb == t) >= 100), dtype=str)
    mins = np.asarray([min(np.sum(la == t), np.sum(lb == t)) for t in clusters], dtype=np.float64)
    return clusters, mins / mins.sum()

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

def fits_stage():
    a, b = load_matched_data()
    ads = {'a': a, 'b': b}
    oa, ob = load_ops('a'), load_ops('b')
    clusters, weights = clusters_and_weights(oa, ob)
    labels = {'a': oa['cell_type'].astype(str), 'b': ob['cell_type'].astype(str)}
    z = {dataset: np.load(LAT[dataset]).astype(np.float32) for dataset in ('a', 'b')}
    if z['a'].shape != (len(a), D) or z['b'].shape != (len(b), D):
        raise RuntimeError('K2 latent shape does not match K2 cell ordering')
    bases = {dataset: tangent(z[dataset]) for dataset in ('a', 'b')}
    all_mats, all_centroids, repro, symmetry_rows = {'a': {}, 'b': {}}, {'a': [], 'b': []}, [], []
    for dataset in ('a', 'b'):
        op = load_ops(dataset)
        for seed in SEEDS:
            drift = EnergyDriftAdapter(D, config=drift_config(seed)).fit(z[dataset], verbose=True)
            check = drift.symmetry_check(z[dataset], n_cells=100)
            torch.save({'state_dict': drift._model.state_dict(), 'config': asdict(drift.config)}, drift_path(dataset, seed))
            m = symmetrize(drift.jacobian_autodiff(z[dataset]).astype(np.float64))
            np.save(mat_path(dataset, seed), m)
            all_mats[dataset][seed] = m
            symmetry_rows.append({'dataset': dataset, 'seed': seed, 'median_asym': float(check['median_asym']), 'max_asym': float(check['max_asym'])})
            all_centroids[dataset].append(np.stack([centroid_operator(z[dataset][labels[dataset] == cluster], m[labels[dataset] == cluster], n_cells_centroid=30).astype(np.float64) for cluster in clusters]))
            if seed == 20261014:
                jr = symmetrize(np.einsum('ndi,ndm,nmj->nij', bases[dataset].astype(np.float64), m, bases[dataset].astype(np.float64)))
                lam = np.linalg.eigvalsh(jr).astype(np.float64)
                repro.append(float(np.max(np.abs(lam - op['lam_seeds'][0].astype(np.float64)))))
    pd.DataFrame(symmetry_rows).to_csv(OUT / 'fit_symmetry.csv', index=False)
    rng_null = np.random.default_rng(20261022)
    null_h1, null_h2 = [], []
    for dataset in ('a', 'b'):
        for seed_b in SEEDS[1:]:
            h1, h2 = [], []
            for cluster in clusters:
                idx = np.flatnonzero(labels[dataset] == cluster)
                perm = rng_null.permutation(idx)
                h = len(perm) // 2
                h1.append(centroid_operator(z[dataset][perm[:h]], all_mats[dataset][20261014][perm[:h]], n_cells_centroid=30).astype(np.float64))
                h2.append(centroid_operator(z[dataset][perm[h:2*h]], all_mats[dataset][seed_b][perm[h:2*h]], n_cells_centroid=30).astype(np.float64))
            null_h1.append(np.stack(h1)); null_h2.append(np.stack(h2))
    rng_split = np.random.default_rng(20261023)
    split = {'a': [[], []], 'b': [[], []]}
    for dataset in ('a', 'b'):
        for cluster in clusters:
            idx = np.flatnonzero(labels[dataset] == cluster)
            perm = rng_split.permutation(idx)
            h = len(perm) // 2
            split[dataset][0].append(centroid_operator(z[dataset][perm[:h]], all_mats[dataset][20261014][perm[:h]], n_cells_centroid=30).astype(np.float64))
            split[dataset][1].append(centroid_operator(z[dataset][perm[h:2*h]], all_mats[dataset][20261014][perm[h:2*h]], n_cells_centroid=30).astype(np.float64))
    np.savez_compressed(centroid_path(), C_a=np.stack(all_centroids['a']), C_b=np.stack(all_centroids['b']), null_h1=np.stack(null_h1), null_h2=np.stack(null_h2), split_a=np.stack([np.stack(split['a'][0]), np.stack(split['a'][1])]), split_b=np.stack([np.stack(split['b'][0]), np.stack(split['b'][1])]), clusters=clusters, w=weights)
    g = {'n_cells_a': int(len(a)), 'n_cells_b': int(len(b)), 'k2_repro_max_abs': float(max(repro)), 'symmetry_max': float(max(row['median_asym'] for row in symmetry_rows))}
    put_gates(g)
    print(json.dumps(g, sort_keys=True))

def load_centroids():
    with np.load(centroid_path(), allow_pickle=False) as x:
        return {key: x[key] for key in x.files}

def load_qs():
    with np.load(LOADINGS, allow_pickle=False) as x:
        wa, wb = x['W_A'].astype(np.float64), x['W_B'].astype(np.float64)
    g = wb @ wa.T
    u, _, vt = np.linalg.svd(g, full_matrices=False)
    q_load = (u @ vt).astype(np.float64)
    u10, s10, vt10 = np.linalg.svd(g[:R, :R], full_matrices=False)
    q_load10 = (u10 @ vt10).astype(np.float64)
    return q_load, q_load10, s10

def residual(a, b, weights):
    w = np.asarray(weights, dtype=np.float64); w = w / w.sum()
    dif = a - b
    return float(np.sqrt(np.sum(w * np.sum(dif * dif, axis=(1,2))) / np.sum(w * np.sum(b * b, axis=(1,2)))))

def part1_stage():
    c = load_centroids(); q_load, q_load10, s10 = load_qs()
    ca, cb, weights, clusters = c['C_a'], c['C_b'], c['w'].astype(np.float64), c['clusters'].astype(str)
    g = get_gates(); g['n_clusters'] = int(len(clusters))
    invalid = g['n_cells_a'] != g['n_cells_b'] or g['k2_repro_max_abs'] > 1e-5 or g['symmetry_max'] > 1e-4
    if invalid:
        g['outcome'] = 'inconclusive'; g['validity_revise'] = True; put_gates(g); print(json.dumps({k:v for k,v in g.items() if k!='validity_revise'}, sort_keys=True)); return
    if g['n_clusters'] < 3:
        g['outcome'] = 'stopped'; put_gates(g); print(json.dumps(g, sort_keys=True)); return
    load_ca = np.einsum('ij,kjl,ml->kim', q_load, ca[0], q_load)
    resid_identity = residual(ca[0, :, :R, :R], cb[0, :, :R, :R], weights)
    resid_load = residual(load_ca[:, :R, :R], cb[0, :, :R, :R], weights)
    null_rows = []
    for i, (dataset, seed_b) in enumerate((d,s) for d in ('a','b') for s in SEEDS[1:]):
        r = residual(c['null_h1'][i, :, :R, :R], c['null_h2'][i, :, :R, :R], weights)
        null_rows.append({'dataset': dataset, 'seed_b': int(seed_b), 'resid': r})
    null_df = pd.DataFrame(null_rows); null_df.to_csv(OUT / 'null.csv', index=False)
    resid_null = float(np.median(null_df.resid))
    g.update({'subspace_min_sv_10': float(np.min(s10)), 'resid_identity_10': resid_identity, 'resid_load_10': resid_load, 'resid_null_10': resid_null, 'frame_ratio_load': float(resid_load / resid_null) if resid_null else float('inf')})
    if g['frame_ratio_load'] <= 1.5: g['outcome'] = 'frame agreement'
    elif g['frame_ratio_load'] > 3.0: g['outcome'] = 'no frame agreement'
    else: g['outcome'] = 'inconclusive'
    put_gates(g)
    part1 = pd.DataFrame([{'quantity': key, 'value': g[key]} for key in ('n_clusters','subspace_min_sv_10','resid_identity_10','resid_load_10','resid_null_10','frame_ratio_load')])
    part1.to_csv(OUT / 'part1.csv', index=False)
    print(json.dumps({key:g[key] for key in ('n_clusters','subspace_min_sv_10','resid_identity_10','resid_load_10','resid_null_10','frame_ratio_load','outcome')}, sort_keys=True))

def r_starts(q_load10):
    yield 0, np.eye(R, dtype=np.float64)
    yield 1, q_load10
    for i in range(1, 19):
        yield i + 1, np.asarray(ortho_group.rvs(dim=R, random_state=20261021 + i), dtype=np.float64)

def loss_r(r, ca, cb, weights):
    w = np.asarray(weights, dtype=np.float64); w = w / w.sum()
    d = np.einsum('ij,kjl,ml->kim', r, ca, r) - cb
    return float(np.sum(w * np.sum(d*d, axis=(1,2))) / np.sum(w * np.sum(cb*cb, axis=(1,2))))

def estimate_r(ca, cb, weights, q_load10, retain=False):
    ca_t = torch.as_tensor(np.asarray(ca, dtype=np.float64), dtype=torch.float64, device='cpu')
    cb_t = torch.as_tensor(np.asarray(cb, dtype=np.float64), dtype=torch.float64, device='cpu')
    wnp = np.asarray(weights, dtype=np.float64); wnp /= wnp.sum()
    wt = torch.as_tensor(wnp, dtype=torch.float64, device='cpu')
    denom = torch.sum(wt * torch.sum(cb_t*cb_t, dim=(1,2)))
    out = []
    for start, r0 in r_starts(q_load10):
        r0t = torch.as_tensor(r0, dtype=torch.float64, device='cpu')
        s = torch.zeros((R,R), dtype=torch.float64, device='cpu', requires_grad=True)
        optimizer = torch.optim.LBFGS([s], lr=1, max_iter=500, tolerance_grad=1e-10, tolerance_change=1e-12, line_search_fn='strong_wolfe')
        def closure():
            optimizer.zero_grad()
            r = r0t @ torch.matrix_exp(s-s.T)
            d = r @ ca_t @ r.T - cb_t
            loss = torch.sum(wt * torch.sum(d*d, dim=(1,2))) / denom
            loss.backward()
            return loss
        optimizer.step(closure)
        with torch.no_grad():
            r = r0t @ torch.matrix_exp(s-s.T)
            d = r @ ca_t @ r.T - cb_t
            loss = float((torch.sum(wt * torch.sum(d*d, dim=(1,2))) / denom).cpu().item())
        out.append((start, loss, r.detach().cpu().numpy().astype(np.float64)))
    best = min(out, key=lambda x:x[1])
    return (best[2], best[1], [(a,b) for a,b,_ in out]) if retain else (best[2], best[1])

def part2_stage():
    g = get_gates()
    if g.get('outcome') == 'stopped': print('Stopped by Part 1 pre-registered rule; Part 2 not run.'); return
    if g.get('validity_revise'): print('Stopped by validity revision condition; Part 2 not run.'); return
    c = load_centroids(); _, q_load10, _ = load_qs()
    ca, cb, weights = c['C_a'][:, :, :R, :R], c['C_b'][:, :, :R, :R], c['w'].astype(np.float64)
    rhat, loss, start_rows = estimate_r(ca[0], cb[0], weights, q_load10, retain=True)
    np.save(OUT/'R.npy', rhat); pd.DataFrame(start_rows, columns=['start','final_loss']).to_csv(OUT/'starts_R.csv', index=False)
    r1, _ = estimate_r(c['split_a'][0, :, :R, :R], c['split_b'][0, :, :R, :R], weights, q_load10)
    r2, _ = estimate_r(c['split_a'][1, :, :R, :R], c['split_b'][1, :, :R, :R], weights, q_load10)
    np.save(OUT/'R_split.npy', np.stack([r1,r2]))
    rs, _ = estimate_r(ca[1], cb[1], weights, q_load10); np.save(OUT/'R_seed.npy', rs)
    starts_df = pd.read_csv(OUT/'starts_R.csv')
    g.update({'loss_best_R': float(loss), 'n_starts_near_best_R': int(np.sum(starts_df.final_loss.to_numpy() <= 1.01*loss)), 'resid_opfit_10': float(np.sqrt(loss)), 'q_op_vs_load_10': float(np.linalg.norm(rhat-q_load10,ord='fro')/np.sqrt(2*R)), 'r_split': float(np.linalg.norm(r1-r2,ord='fro')/np.sqrt(2*R)), 'r_seed': float(np.linalg.norm(rhat-rs,ord='fro')/np.sqrt(2*R))})
    put_gates(g)
    print(json.dumps({key:g[key] for key in ('loss_best_R','n_starts_near_best_R','resid_opfit_10','q_op_vs_load_10','r_split','r_seed')}, sort_keys=True))

def gene_matrix(ad, genes):
    xad=ad[:,genes].copy(); sc.pp.normalize_total(xad,target_sum=1e4); sc.pp.log1p(xad)
    return xad.X.toarray() if hasattr(xad.X,'toarray') else np.asarray(xad.X)

def make_umap(rep):
    ad=sc.AnnData(X=np.zeros((len(rep),1),dtype=np.float32)); ad.obsm['view']=np.asarray(rep)
    sc.pp.neighbors(ad,n_neighbors=15,use_rep='view',random_state=20261014); sc.tl.umap(ad,min_dist=0.5,random_state=20261014)
    return np.asarray(ad.obsm['X_umap'],dtype=np.float64)

def mixing(rep,dataset):
    n=NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep,return_distance=False)[:,1:]
    return float(np.mean(np.mean(dataset[n]!=dataset[:,None],axis=1))/0.5)

def plot_umap(coords,dataset,types,scores):
    fig,axes=plt.subplots(3,2,figsize=(10,13)); cmap=plt.get_cmap('tab10'); labels=sorted(set(types.astype(str)))
    for row,name in enumerate(('joint','sep_noA','sep_load')):
        xy=coords[name]; mix,sil=scores[name]; ax=axes[row,0]
        for ds,color in (('A','C0'),('B','C1')):
            m=dataset==ds;ax.scatter(xy[m,0],xy[m,1],s=4,alpha=.7,color=color,label=ds,linewidths=0)
        ax.set_title(f'{name}: dataset; mix={mix:.3f}, sil={sil:.3f}');ax.legend(markerscale=2,fontsize=8);ax.set_xticks([]);ax.set_yticks([])
        ax=axes[row,1]
        for i,label in enumerate(labels):
            m=types==label;ax.scatter(xy[m,0],xy[m,1],s=4,alpha=.7,color=cmap(i%10),label=label,linewidths=0)
        ax.set_title(f'{name}: cell type; mix={mix:.3f}, sil={sil:.3f}');ax.legend(markerscale=2,fontsize=7,ncol=2);ax.set_xticks([]);ax.set_yticks([])
    fig.tight_layout();fig.savefig(OUT/'fig_umap.png',dpi=180);plt.close(fig)

def plot_frames(ca,cb,qload,clusters):
    fig,ax=plt.subplots(figsize=(8,5));x=np.arange(1,11);cmap=plt.get_cmap('tab20')
    for i,cluster in enumerate(clusters):
        aligned=(qload@ca[i]@qload.T)[:R,:R]; plain=ca[i,:R,:R]; target=cb[i,:R,:R]
        _,ua=np.linalg.eigh(aligned);_,ub=np.linalg.eigh(target);_,ui=np.linalg.eigh(plain)
        y=np.abs(np.sum(ua[:,-10:]*ub[:,-10:],axis=0));y0=np.abs(np.sum(ui[:,-10:]*ub[:,-10:],axis=0));color=cmap(i%20)
        ax.plot(x,y,marker='o',markersize=3,color=color,label=str(cluster));ax.plot(x,y0,linestyle='--',color=color,alpha=.6,label=f'{cluster} (Q=I)')
    ax.set_xlabel('Eigenvector rank among largest 10 (ascending)');ax.set_ylabel('Absolute cosine with B eigenvector');ax.set_ylim(-.02,1.02);ax.legend(fontsize=6,ncol=2);fig.tight_layout();fig.savefig(OUT/'fig_frames.png',dpi=180);plt.close(fig)

def views_stage():
    g=get_gates()
    if g.get('outcome')=='stopped' or g.get('validity_revise'): print('Views not run after pre-registered stop or validity revision.'); return
    a,b=load_matched_data(); qload,_,_=load_qs(); genes=[x.strip() for x in GENES.read_text().splitlines() if x.strip()]
    if len(genes)!=2000:raise RuntimeError('T01d genes list must have 2000 genes')
    joint=PCA(n_components=D,svd_solver='full',random_state=20261014).fit_transform(np.vstack([gene_matrix(a,genes),gene_matrix(b,genes)]))
    za,zb=np.load(LAT['a']).astype(np.float64),np.load(LAT['b']).astype(np.float64)
    reps={'joint':joint,'sep_noA':np.vstack([za,zb]),'sep_load':np.vstack([za@qload.T,zb])}
    dataset=np.asarray(['A']*len(a)+['B']*len(b),dtype=str); types=np.concatenate([a.obs['cell_type'].astype(str).to_numpy(),b.obs['cell_type'].astype(str).to_numpy()]).astype(str)
    oa,ob=load_ops('a'),load_ops('b'); tested=sorted(t for t in set(oa['cell_type'].astype(str))&set(ob['cell_type'].astype(str)) if np.sum(oa['cell_type'].astype(str)==t)>=50 and np.sum(ob['cell_type'].astype(str)==t)>=50);mask=np.isin(types,tested)
    coords={};scores={}
    for name,rep in reps.items():coords[name]=make_umap(rep);scores[name]=(mixing(rep,dataset),float(silhouette_score(rep[mask],types[mask],metric='euclidean')))
    np.savez_compressed(OUT/'umap_views.npz',joint=coords['joint'],sep_noA=coords['sep_noA'],sep_load=coords['sep_load'],dataset=dataset,cell_type=types)
    plot_umap(coords,dataset,types,scores);c=load_centroids();plot_frames(c['C_a'][0],c['C_b'][0],qload,c['clusters'].astype(str))
    g.update({'mix_joint':scores['joint'][0],'mix_sep_noA':scores['sep_noA'][0],'mix_sep_load':scores['sep_load'][0],'sil_joint':scores['joint'][1],'sil_sep_noA':scores['sep_noA'][1],'sil_sep_load':scores['sep_load'][1]});put_gates(g)
    print(json.dumps({key:g[key] for key in ('mix_joint','mix_sep_noA','mix_sep_load','sil_joint','sil_sep_noA','sil_sep_load')},sort_keys=True))

def report_stage():
    g=get_gates();stopped=g.get('outcome')=='stopped'
    files=['requirements.lock','centroid_ops.npz','scripts/t02c_run.py']
    if (OUT/'part1.csv').exists():files+=['part1.csv','null.csv']
    if (OUT/'R.npy').exists():files+=['R.npy','R_split.npy','R_seed.npy','starts_R.csv']
    if not stopped and (OUT/'umap_views.npz').exists():files+=['umap_views.npz','fig_umap.png','fig_frames.png']
    command='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T02c/scripts/t02c_run.py'
    fits=('n_cells_a','n_cells_b','k2_repro_max_abs','symmetry_max');part1=('n_clusters','subspace_min_sv_10','resid_identity_10','resid_load_10','resid_null_10','frame_ratio_load','outcome');part2=('loss_best_R','n_starts_near_best_R','resid_opfit_10','q_op_vs_load_10','r_split','r_seed');views=('mix_joint','mix_sep_noA','mix_sep_load','sil_joint','sil_sep_noA','sil_sep_load')
    regen={path:command+' --stage all' for path in files};regen.update({key:command+' --stage fits' for key in fits if key in g});regen.update({key:command+' --stage part1' for key in part1 if key in g});regen.update({key:command+' --stage part2' for key in part2 if key in g});regen.update({key:command+' --stage views' for key in views if key in g})
    obj={'task_id':'T02c','executor':'manus','branch':'task/T02c','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':'2026-10-15T15:38:00-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'stopped' if stopped else 'complete','gated':{key:value for key,value in g.items() if key!='validity_revise'},'deliverables':[{'path':path,'sha256':sha(OUT/path)} for path in files],'regeneration':regen,'deviations':[],'deviations_requested':[],'observations':[f'Input SHA-256: {path.relative_to(ROOT)} = {sha(path)}.' for path in INPUTS],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')}
    if stopped:obj['stop_rule']='n_clusters < 3 → inconclusive (stop).'
    (OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n');print(json.dumps(obj['gated'],sort_keys=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('fits','part1','part2','views','report','all'));stage=p.parse_args().stage
    if stage in ('fits','all'):fits_stage()
    if stage in ('part1','all'):part1_stage()
    if stage in ('part2','all'):part2_stage()
    if stage in ('views','all'):views_stage()
    if stage in ('report','all'):report_stage()
if __name__=='__main__':main()
