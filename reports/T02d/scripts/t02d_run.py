#!/usr/bin/env python
import argparse
import hashlib
import json
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'reports' / 'T02d'
K2 = ROOT / 'reports' / 'K2'
T02B = ROOT / 'reports' / 'T02b'
T02C = ROOT / 'reports' / 'T02c'
T01D = ROOT / 'reports' / 'T01d'
DATA = ROOT / 'data' / 'pbmc_v2_v3'
LAT = {'a': K2 / 'latent_pca30_a.npy', 'b': K2 / 'latent_pca30_b.npy'}
OPS = {'a': K2 / 'ops_a.npz', 'b': K2 / 'ops_b.npz'}
LOADINGS = T02B / 'loadings.npz'
T02C_REPORT = T02C / 'report.json'
V2, V3_THIN, CELLS, GENES = DATA / 'pbmc_v2.h5ad', DATA / 'pbmc_v3_thin.h5ad', DATA / 'k2_cells.csv', T01D / 'genes.txt'
INPUTS = (LAT['a'], LAT['b'], OPS['a'], OPS['b'], LOADINGS, CELLS, V2, V3_THIN, GENES, T02C_REPORT)
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
        for chunk in iter(lambda: f.read(1 << 20), b''): h.update(chunk)
    return h.hexdigest()
def gates_path(): return OUT / 'gates.json'
def get_gates(): return json.loads(gates_path().read_text())
def put_gates(g): gates_path().write_text(json.dumps(g, sort_keys=True) + '\n')
def mat_path(dataset, seed): return OUT / f'M_{dataset}_{seed}.npy'
def drift_path(dataset, seed): return OUT / f'drift_{dataset}_{seed}.pt'
def floor_path(): return OUT / 'floor_ops.npz'

def load_ops(dataset):
    with np.load(OPS[dataset], allow_pickle=False) as x: return {key:x[key] for key in x.files}

def load_matched_data():
    cells = pd.read_csv(CELLS)
    a0 = sc.read_h5ad(V2); a0 = a0[a0.obs['qc_pass'].to_numpy(dtype=bool)].copy()
    b0 = sc.read_h5ad(V3_THIN)
    bars_a = cells.loc[cells.dataset == 'A', 'barcode'].astype(str).to_numpy(); bars_b = cells.loc[cells.dataset == 'B', 'barcode'].astype(str).to_numpy()
    a,b = a0[bars_a].copy(), b0[bars_b].copy()
    for dataset,ad,bars in (('a',a,bars_a),('b',b,bars_b)):
        op=load_ops(dataset)
        if not np.array_equal(ad.obs_names.astype(str).to_numpy(),bars): raise RuntimeError(f'cell order mismatch {dataset}')
        if not np.array_equal(ad.obs_names.astype(str).to_numpy(),op['barcode'].astype(str)): raise RuntimeError(f'K2 barcodes mismatch {dataset}')
        if not np.array_equal(ad.obs['cell_type'].astype(str).to_numpy(),op['cell_type'].astype(str)): raise RuntimeError(f'K2 labels mismatch {dataset}')
    return a,b

def clusters_and_weights(oa,ob):
    la,lb=oa['cell_type'].astype(str),ob['cell_type'].astype(str)
    clusters=np.asarray(sorted(t for t in set(la)&set(lb) if np.sum(la==t)>=100 and np.sum(lb==t)>=100),dtype=str)
    mins=np.asarray([min(np.sum(la==t),np.sum(lb==t)) for t in clusters],dtype=np.float64)
    return clusters,mins/mins.sum()

def tangent(z):
    nbr=NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z,return_distance=False)[:,1:]
    b=np.empty((len(z),D,R),dtype=np.float32)
    for i in range(len(z)):
        _,u=np.linalg.eigh(np.cov(z[nbr[i]],rowvar=False));b[i]=u[:,-R:]
    return b

def config(seed): return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=seed,batch_size=256,loss='esm',lr=5e-4,beta=1.0,esm_laplacian='hutchinson',device='cpu')
def symmetric(j): return (j+np.swapaxes(j,1,2))/2.0

def floor_operator(indices, energies, matrices):
    chosen=indices[np.argsort(energies[indices],kind='stable')[:30]]
    return np.asarray(matrices[chosen].mean(axis=0),dtype=np.float64),chosen

def centroid_indices(z,indices):
    sub=z[indices];mu=sub.mean(axis=0);dist=np.linalg.norm(sub-mu,axis=1)
    return indices[np.argsort(dist,kind='stable')[:30]]

def fits_stage():
    a,b=load_matched_data();oa,ob=load_ops('a'),load_ops('b');clusters,w=clusters_and_weights(oa,ob)
    labels={'a':oa['cell_type'].astype(str),'b':ob['cell_type'].astype(str)};z={'a':np.load(LAT['a']).astype(np.float32),'b':np.load(LAT['b']).astype(np.float32)}
    if z['a'].shape!=(len(a),D) or z['b'].shape!=(len(b),D):raise RuntimeError('K2 latent shape mismatch')
    bases={dataset:tangent(z[dataset]) for dataset in ('a','b')}
    mats={'a':{},'b':{}};energies={'a':{},'b':{}};floors={'a':[],'b':[]};repro=[];sym=[];floor_idx={'a':None,'b':None}
    for dataset in ('a','b'):
        op=load_ops(dataset)
        cluster_indices=[np.flatnonzero(labels[dataset]==cluster) for cluster in clusters]
        for seed in SEEDS:
            drift=EnergyDriftAdapter(D,config=config(seed)).fit(z[dataset],verbose=True)
            check=drift.symmetry_check(z[dataset],n_cells=100)
            torch.save({'state_dict':drift._model.state_dict(),'config':asdict(drift.config)},drift_path(dataset,seed))
            m=symmetric(drift.jacobian_autodiff(z[dataset]).astype(np.float64));e=drift.energy(z[dataset]).astype(np.float64)
            np.save(mat_path(dataset,seed),m);mats[dataset][seed]=m;energies[dataset][seed]=e
            sym.append({'dataset':dataset,'seed':seed,'median_asym':float(check['median_asym']),'max_asym':float(check['max_asym'])})
            these=[];idx_rows=[]
            for indices in cluster_indices:
                f,idx=floor_operator(indices,e,m);these.append(f);idx_rows.append(idx)
            floors[dataset].append(np.stack(these))
            if seed==20261014:floor_idx[dataset]=np.stack(idx_rows)
            if seed==20261014:
                jr=symmetric(np.einsum('ndi,ndm,nmj->nij',bases[dataset].astype(np.float64),m,bases[dataset].astype(np.float64)))
                repro.append(float(np.max(np.abs(np.linalg.eigvalsh(jr).astype(np.float64)-op['lam_seeds'][0].astype(np.float64)))))
    pd.DataFrame(sym).to_csv(OUT/'fit_symmetry.csv',index=False)
    rng=np.random.default_rng(20261022);null_h1=[];null_h2=[]
    for dataset in ('a','b'):
        for seed_b in SEEDS[1:]:
            h1=[];h2=[]
            for cluster in clusters:
                indices=np.flatnonzero(labels[dataset]==cluster);perm=rng.permutation(indices);h=len(perm)//2
                f1,_=floor_operator(perm[:h],energies[dataset][20261014],mats[dataset][20261014]);f2,_=floor_operator(perm[h:2*h],energies[dataset][seed_b],mats[dataset][seed_b]);h1.append(f1);h2.append(f2)
            null_h1.append(np.stack(h1));null_h2.append(np.stack(h2))
    overlaps=[]
    for dataset in ('a','b'):
        for i,cluster in enumerate(clusters):
            indices=np.flatnonzero(labels[dataset]==cluster);center=centroid_indices(z[dataset],indices);inter=np.intersect1d(floor_idx[dataset][i],center,assume_unique=False)
            overlaps.append(len(inter)/30.0)
    np.savez_compressed(floor_path(),F_a=np.stack(floors['a']),F_b=np.stack(floors['b']),null_h1=np.stack(null_h1),null_h2=np.stack(null_h2),floor_idx_a=floor_idx['a'].astype(np.int64),floor_idx_b=floor_idx['b'].astype(np.int64),clusters=clusters,w=w)
    g={'n_cells_a':int(len(a)),'n_cells_b':int(len(b)),'k2_repro_max_abs':float(max(repro)),'symmetry_max':float(max(row['median_asym'] for row in sym)),'floor_centroid_overlap':float(np.mean(overlaps))}
    put_gates(g);print(json.dumps(g,sort_keys=True))

def load_floor():
    with np.load(floor_path(),allow_pickle=False) as x:return {key:x[key] for key in x.files}

def q_load():
    with np.load(LOADINGS,allow_pickle=False) as x:wa,wb=x['W_A'].astype(np.float64),x['W_B'].astype(np.float64)
    u,_,vt=np.linalg.svd(wb@wa.T,full_matrices=False);return (u@vt).astype(np.float64)

def residual(a,b,w):
    weight=np.asarray(w,dtype=np.float64);weight/=weight.sum();d=a-b
    return float(np.sqrt(np.sum(weight*np.sum(d*d,axis=(1,2)))/np.sum(weight*np.sum(b*b,axis=(1,2)))))

def part1_stage():
    f=load_floor();g=get_gates();g['n_clusters']=int(len(f['clusters']))
    invalid=g['n_cells_a']!=g['n_cells_b'] or g['k2_repro_max_abs']>1e-5 or g['symmetry_max']>1e-4
    if invalid:
        g['outcome']='inconclusive';g['validity_revise']=True;put_gates(g);print(json.dumps({k:v for k,v in g.items() if k!='validity_revise'},sort_keys=True));return
    if g['n_clusters']<3:
        g['outcome']='stopped';put_gates(g);print(json.dumps(g,sort_keys=True));return
    q=q_load();fa,fb,w=f['F_a'][0],f['F_b'][0],f['w'].astype(np.float64)
    rotated=np.einsum('ij,kjl,ml->kim',q,fa,q)
    rid=residual(fa[:,:R,:R],fb[:,:R,:R],w);rload=residual(rotated[:,:R,:R],fb[:,:R,:R],w)
    rows=[]
    for i,(dataset,seed_b) in enumerate((d,s) for d in ('a','b') for s in SEEDS[1:]):
        r=residual(f['null_h1'][i,:,:R,:R],f['null_h2'][i,:,:R,:R],w);rows.append({'dataset':dataset,'seed_b':int(seed_b),'resid':r})
    ndf=pd.DataFrame(rows);ndf.to_csv(OUT/'null.csv',index=False);rnull=float(np.median(ndf.resid))
    prior=json.loads(T02C_REPORT.read_text());prior_ratio=prior.get('gated',{}).get('frame_ratio_load')
    g.update({'resid_identity_10':rid,'resid_load_10':rload,'resid_null_10':rnull,'frame_ratio_load':float(rload/rnull) if rnull else float('inf')})
    if prior_ratio is not None:g['floor_vs_centroid']=float(g['frame_ratio_load']/prior_ratio) if prior_ratio else float('inf')
    if g['frame_ratio_load']<=1.5:g['outcome']='frame agreement'
    elif g['frame_ratio_load']>3.0:g['outcome']='no frame agreement'
    else:g['outcome']='inconclusive'
    put_gates(g)
    keys=['n_clusters','resid_identity_10','resid_load_10','resid_null_10','frame_ratio_load']+(['floor_vs_centroid'] if 'floor_vs_centroid' in g else [])+['floor_centroid_overlap']
    pd.DataFrame([{'quantity':key,'value':g[key]} for key in keys]).to_csv(OUT/'part1.csv',index=False)
    print(json.dumps({key:g[key] for key in keys+['outcome']},sort_keys=True))

def report_stage():
    g=get_gates();stopped=g.get('outcome')=='stopped';files=['requirements.lock','floor_ops.npz','scripts/t02d_run.py']
    if (OUT/'part1.csv').exists():files+=['part1.csv','null.csv']
    command='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T02d/scripts/t02d_run.py'
    fits=('n_cells_a','n_cells_b','k2_repro_max_abs','symmetry_max');part1=('n_clusters','resid_identity_10','resid_load_10','resid_null_10','frame_ratio_load','floor_vs_centroid','floor_centroid_overlap','outcome')
    regen={path:command+' --stage all' for path in files};regen.update({key:command+' --stage fits' for key in fits if key in g});regen.update({key:command+' --stage part1' for key in part1 if key in g})
    observations=[f'Input SHA-256: {path.relative_to(ROOT)} = {sha(path)}.' for path in INPUTS]
    if 'floor_vs_centroid' not in g:observations.append('T02c did not provide frame_ratio_load, so floor_vs_centroid was omitted.')
    obj={'task_id':'T02d','executor':'manus','branch':'task/T02d','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':'2026-10-15T15:43:00-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'stopped' if stopped else 'complete','gated':{key:value for key,value in g.items() if key!='validity_revise'},'deliverables':[{'path':path,'sha256':sha(OUT/path)} for path in files],'regeneration':regen,'deviations':[],'deviations_requested':[],'observations':observations,'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')}
    if stopped:obj['stop_rule']='n_clusters < 3 → inconclusive (stop).'
    (OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n');print(json.dumps(obj['gated'],sort_keys=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('fits','part1','report','all'));stage=p.parse_args().stage
    if stage in ('fits','all'):fits_stage()
    if stage in ('part1','all'):part1_stage()
    if stage in ('report','all'):report_stage()
if __name__=='__main__':main()
