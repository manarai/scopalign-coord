#!/usr/bin/env python
import argparse, hashlib, itertools, json, subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import scanpy as sc
import torch
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'reports'/'K2'; DATA=ROOT/'data'/'pbmc_v2_v3'; T01D=ROOT/'reports'/'T01d'
V2=DATA/'pbmc_v2.h5ad'; V3=DATA/'pbmc_v3_thin.h5ad'; GENES=T01D/'genes.txt'; CELLS=DATA/'k2_cells.csv'; SEEDS=(20261014,20261016,20261017,20261018,20261019); MATCH_SEED=20261014; MAIN_SEED=20261015; REFIT_SEED=20261020; R=10
INPUTS=(V2,V3,GENES); CLASSES=('', 'holder','fork','drift'); K1REF={'spectral_ratio_refit':1.5294658061679813,'basin_ratio_refit':0.8718013824357383,'class_tv_ratio_refit':1.0941906957865921}
KNOWN=["tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]","tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility"]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def lat(d):return OUT/f'latent_pca30_{d}.npy'
def bas(d):return OUT/f'basis_{d}.npy'
def rank(d):return OUT/f'rank_{d}.npy'
def ops(d):return OUT/f'ops_{d}.npz'
def drift_path(d,s):return OUT/f'drift_{d}_{s}.pt'
def indices_path():return OUT/'match_indices.npz'
def gates_path():return OUT/'gates.json'
def depth(a):
 x=a.X.sum(axis=1);return np.asarray(x).reshape(-1).astype(np.float64)
def source_data():
 a=sc.read_h5ad(V2);a=a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy();b=sc.read_h5ad(V3);return a,b
def matched_data():
 a,b=source_data();z=np.load(indices_path(),allow_pickle=False);return a[z['a']],b[z['b']]
def match_stage():
 a,b=source_data();la=a.obs['cell_type'].astype(str).to_numpy();lb=b.obs['cell_type'].astype(str).to_numpy();rng=np.random.default_rng(MATCH_SEED);ia=[];ib=[];rows=[]
 for t in sorted(set(la)&set(lb)):
  xa=np.flatnonzero(la==t);xb=np.flatnonzero(lb==t);n=min(len(xa),len(xb));sa=rng.choice(xa,n,replace=False);sb=rng.choice(xb,n,replace=False);ia.extend(sa.tolist());ib.extend(sb.tolist())
 ia=np.asarray(sorted(ia),dtype=np.int64);ib=np.asarray(sorted(ib),dtype=np.int64);ma,mb=a[ia],b[ib];np.savez_compressed(indices_path(),a=ia,b=ib)
 pd.DataFrame({'barcode':np.concatenate([ma.obs_names.astype(str).to_numpy(),mb.obs_names.astype(str).to_numpy()]),'dataset':['A']*len(ma)+['B']*len(mb)}).to_csv(CELLS,index=False)
 da,db=depth(ma),depth(mb);mla=ma.obs['cell_type'].astype(str).to_numpy();mlb=mb.obs['cell_type'].astype(str).to_numpy()
 for t in sorted(set(la)&set(lb)):
  x0=la==t;x1=lb==t;y0=mla==t;y1=mlb==t;rows.append({'cell_type':t,'n_v2_qc':int(x0.sum()),'n_v3thin':int(x1.sum()),'n_matched':int(y0.sum()),'median_depth_a':float(np.median(da[y0])),'median_depth_b':float(np.median(db[y1]))})
 tab=pd.DataFrame(rows);tab.to_csv(OUT/'match_table.csv',index=False);tested=tab[(tab.n_matched>=50)];rel=np.abs(tested.median_depth_b-tested.median_depth_a)/tested.median_depth_a;g={'n_cells_a':int(len(ma)),'n_cells_b':int(len(mb)),'count_match_max_abs':int(max((abs(int((mla==t).sum())-int((mlb==t).sum())) for t in tab.cell_type),default=0)),'depth_match_max_rel':float(rel.max()) if len(rel) else float('inf')};gates_path().write_text(json.dumps(g,sort_keys=True)+'\n');print(json.dumps(g,sort_keys=True))
def norm_pca(a,genes):
 missing=[g for g in genes if g not in a.var_names]
 if missing:raise RuntimeError(f'missing requested genes: {len(missing)}')
 a=a[:,genes].copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);x=a.X
 if hasattr(x,'toarray'):x=x.toarray()
 return a,PCA(n_components=30,svd_solver='full',random_state=20261014).fit_transform(np.asarray(x,dtype=np.float32)).astype(np.float32)
def tangent(z):
 nbr=NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z,return_distance=False)[:,1:];v=np.empty((len(z),30,R),dtype=np.float32);r=np.empty(len(z),dtype=np.int32)
 for i in range(len(z)):
  w,u=np.linalg.eigh(np.cov(z[nbr[i]],rowvar=False));r[i]=int(np.sum(w>.05*w.max()));v[i]=u[:,-R:]
 return r,v
def latent_stage():
 genes=[x.strip() for x in GENES.read_text().splitlines() if x.strip()]
 if len(genes)!=2000 or genes!=sorted(genes):raise RuntimeError('genes.txt is not the specified sorted 2000-gene list')
 a,b=matched_data();g=json.loads(gates_path().read_text());rs=[]
 for d,ad in (('a',a),('b',b)):
  _,z=norm_pca(ad,genes);r,v=tangent(z);np.save(lat(d),z);np.save(rank(d),r);np.save(bas(d),v);g[f'rank_pca_{d}']=float(np.median(r));rs.append(r)
 g['frac_rank_ge10']=float(np.mean(np.concatenate(rs)>=10));gates_path().write_text(json.dumps(g,sort_keys=True)+'\n');print(json.dumps({k:g[k] for k in ('rank_pca_a','rank_pca_b','frac_rank_ge10')},sort_keys=True))
def cfg(seed):return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=seed,batch_size=256,loss='esm',lr=5e-4,beta=1.,esm_laplacian='hutchinson',device='cpu')
def sym(j):return (j+np.swapaxes(j,1,2))/2
def classify(lam):
 tol=.1*float(np.median(np.abs(lam)));top=lam.max(axis=1);return (-top).astype(np.float64),np.where(top>tol,'fork',np.where(top<-tol,'holder','drift')).astype(str),tol
def fit_stage():
 a,b=matched_data();symrows=[]
 for d,ad in (('a',a),('b',b)):
  z=np.load(lat(d)).astype(np.float32);v=np.load(bas(d)).astype(np.float64);lams=[];bis=[];clss=[];tols=[]
  for seed in SEEDS:
   dr=EnergyDriftAdapter(30,config=cfg(seed)).fit(z,verbose=True);check=dr.symmetry_check(z,n_cells=100);torch.save({'state_dict':dr._model.state_dict(),'config':asdict(dr.config)},drift_path(d,seed));j=sym(dr.jacobian_autodiff(z).astype(np.float64));jr=sym(np.einsum('ndi,ndm,nmj->nij',v,j,v));lam=np.linalg.eigvalsh(jr).astype(np.float64);bi,cls,tol=classify(lam);lams.append(lam);bis.append(bi);clss.append(cls);tols.append(tol);symrows.append({'dataset':d,'seed':seed,'median_asym':float(check['median_asym']),'max_asym':float(check['max_asym'])})
  np.savez_compressed(ops(d),lam_seeds=np.stack(lams).astype(np.float64),basin_index_seeds=np.stack(bis).astype(np.float64),cls_seeds=np.asarray(clss,dtype=str),cell_type=np.asarray(ad.obs['cell_type'].astype(str).to_numpy(),dtype=str),barcode=np.asarray(ad.obs_names.astype(str).to_numpy(),dtype=str),tol_seeds=np.asarray(tols,dtype=np.float64),symmetry_median=np.asarray([x['median_asym'] for x in symrows if x['dataset']==d],dtype=np.float64),seeds=np.asarray(SEEDS,dtype=np.int64))
 pd.DataFrame(symrows).to_csv(OUT/'fit_symmetry.csv',index=False);g=json.loads(gates_path().read_text());g['symmetry_max']=float(max(x['median_asym'] for x in symrows));gates_path().write_text(json.dumps(g,sort_keys=True)+'\n');print(json.dumps({'symmetry_max':g['symmetry_max']},sort_keys=True))
def readops(d):
 with np.load(ops(d),allow_pickle=False) as z:return {k:z[k] for k in z.files}
def w2(x,y):return float(np.sqrt(ot.wasserstein_1d(x,y,p=2)))
def w1(x,y):return float(ot.wasserstein_1d(x,y,p=1))
def tv(x,y):return float(.5*sum(abs(np.mean(x==q)-np.mean(y==q)) for q in CLASSES))
def quot(a,b):return 1. if a==0 and b==0 else (float('inf') if b==0 else float(a/b))
def majority(x):
 c=pd.Series(x).value_counts();w=c[c==c.max()].index.astype(str).tolist();return w[0] if len(w)==1 else ''
def types_m(o):
 types=sorted(t for t in set(o['a']['cell_type'].astype(str))&set(o['b']['cell_type'].astype(str)) if (o['a']['cell_type'].astype(str)==t).sum()>=50 and (o['b']['cell_type'].astype(str)==t).sum()>=50);idx={d:{t:np.flatnonzero(o[d]['cell_type'].astype(str)==t) for t in types} for d in ('a','b')};m={t:min(len(idx['a'][t]),len(idx['b'][t]))//2 for t in types};return types,idx,m
def mean_metrics(l1,b1,c1,l2,b2,c2,pairs):
 sp=[];bi=[];cl=[]
 for x,y in pairs:sp.append(w2(l1[x].ravel(),l2[y].ravel()));bi.append(w1(b1[x],b2[y]));cl.append(tv(c1[x],c2[y]))
 return float(np.mean(sp)),float(np.mean(bi)),float(np.mean(cl))
def compare_stage():
 o={d:readops(d) for d in ('a','b')};types,idx,m=types_m(o);rng=np.random.default_rng(MAIN_SEED);rows=[]
 la,ba,ca=o['a']['lam_seeds'][0],o['a']['basin_index_seeds'][0],o['a']['cls_seeds'][0];lb,bb,cb=o['b']['lam_seeds'][0],o['b']['basin_index_seeds'][0],o['b']['cls_seeds'][0]
 for d in range(200):
  ia={t:rng.choice(idx['a'][t],m[t],replace=False) for t in types};ib={t:rng.choice(idx['b'][t],m[t],replace=False) for t in types};sp,bi,cl=mean_metrics(la,ba,ca,lb,bb,cb,[(ia[t],ib[t]) for t in types]);rows.append({'kind':'pooled','draw':d,'seed_b':20261014,'spectral_w2':sp,'basin_w1':bi,'class_tv':cl})
 for which,l,b,c in (('a',la,ba,ca),('b',lb,bb,cb)):
  for d in range(100):
   q={t:rng.choice(idx[which][t],2*m[t],replace=False) for t in types};sp,bi,cl=mean_metrics(l,b,c,l,b,c,[(q[t][:m[t]],q[t][m[t]:]) for t in types]);rows.append({'kind':f'null_{which}','draw':d,'seed_b':20261014,'spectral_w2':sp,'basin_w1':bi,'class_tv':cl})
 for which,l,b,c in (('a',la,ba,ca),('b',lb,bb,cb)):
  for d in range(100):
   pairs=[]
   for u,v in itertools.combinations(types,2):
    n=min(m[u],m[v]);pairs.append((rng.choice(idx[which][u],n,replace=False),rng.choice(idx[which][v],n,replace=False)))
   sp,bi,cl=mean_metrics(l,b,c,l,b,c,pairs);rows.append({'kind':f'sep_{which}','draw':d,'seed_b':20261014,'spectral_w2':sp,'basin_w1':bi,'class_tv':cl})
 rng=np.random.default_rng(REFIT_SEED)
 for which in ('a','b'):
  for d in range(100):
   si=1+d%4;q={t:rng.choice(idx[which][t],2*m[t],replace=False) for t in types};l0,b0,c0=o[which]['lam_seeds'][0],o[which]['basin_index_seeds'][0],o[which]['cls_seeds'][0];l1,b1,c1=o[which]['lam_seeds'][si],o[which]['basin_index_seeds'][si],o[which]['cls_seeds'][si];sp,bi,cl=mean_metrics(l0,b0,c0,l1,b1,c1,[(q[t][:m[t]],q[t][m[t]:]) for t in types]);rows.append({'kind':f'refit_{which}','draw':d,'seed_b':SEEDS[si],'spectral_w2':sp,'basin_w1':bi,'class_tv':cl})
 df=pd.DataFrame(rows);df.to_csv(OUT/'draws.csv',index=False);med=lambda mask,col:float(np.median(df.loc[mask,col]));po={c:med(df.kind=='pooled',c) for c in ('spectral_w2','basin_w1','class_tv')};rf={c:med(df.kind.str.startswith('refit_'),c) for c in po};nu={c:med(df.kind.str.startswith('null_'),c) for c in po};g=json.loads(gates_path().read_text());g.update({'n_cell_types_tested':len(types),'spectral_sep_a':quot(med(df.kind=='sep_a','spectral_w2'),med(df.kind=='null_a','spectral_w2')),'spectral_sep_b':quot(med(df.kind=='sep_b','spectral_w2'),med(df.kind=='null_b','spectral_w2')),'basin_sep_a':quot(med(df.kind=='sep_a','basin_w1'),med(df.kind=='null_a','basin_w1')),'basin_sep_b':quot(med(df.kind=='sep_b','basin_w1'),med(df.kind=='null_b','basin_w1')),'spectral_w2_pooled':po['spectral_w2'],'basin_w1_pooled':po['basin_w1'],'class_tv_pooled':po['class_tv'],'spectral_w2_refit':rf['spectral_w2'],'basin_w1_refit':rf['basin_w1'],'class_tv_refit':rf['class_tv'],'spectral_ratio_refit':quot(po['spectral_w2'],rf['spectral_w2']),'basin_ratio_refit':quot(po['basin_w1'],rf['basin_w1']),'class_tv_ratio_refit':quot(po['class_tv'],rf['class_tv']),'spectral_ratio_withinfit':quot(po['spectral_w2'],nu['spectral_w2']),'basin_ratio_withinfit':quot(po['basin_w1'],nu['basin_w1']),'class_tv_ratio_withinfit':quot(po['class_tv'],nu['class_tv'])})
 g['spectral_vs_k1']=quot(g['spectral_ratio_refit'],K1REF['spectral_ratio_refit']);g['basin_vs_k1']=quot(g['basin_ratio_refit'],K1REF['basin_ratio_refit']);g['class_tv_vs_k1']=quot(g['class_tv_ratio_refit'],K1REF['class_tv_ratio_refit'])
 if g['n_cell_types_tested']<3:g['outcome']='inconclusive'
 elif min(g['spectral_sep_a'],g['spectral_sep_b'],g['basin_sep_a'],g['basin_sep_b'])<2:g['outcome']='inconclusive'
 elif g['spectral_ratio_refit']>3 or g['basin_ratio_refit']>3:g['outcome']='fails'
 elif g['spectral_ratio_refit']<=1.5 and g['basin_ratio_refit']<=1.5 and g['class_tv_ratio_refit']<=1.5:g['outcome']='survives'
 else:g['outcome']='inconclusive'
 tab=[]
 for t in sorted(set(o['a']['cell_type'].astype(str))|set(o['b']['cell_type'].astype(str))):
  ma=o['a']['cell_type'].astype(str)==t;mb=o['b']['cell_type'].astype(str)==t;tab.append({'cell_type':t,'n_a':int(ma.sum()),'n_b':int(mb.sum()),'median_bi_a':float(np.median(ba[ma])) if ma.any() else np.nan,'median_bi_b':float(np.median(bb[mb])) if mb.any() else np.nan,'class_a':majority(ca[ma]) if ma.any() else '', 'class_b':majority(cb[mb]) if mb.any() else ''})
 table=pd.DataFrame(tab);table.to_csv(OUT/'basin_table.csv',index=False)
 fig,axs=plt.subplots(1,len(types),figsize=(3*len(types),3.2),squeeze=False)
 for ax,t in zip(axs[0],types):
  for d,l,col in (('a',la,'C0'),('b',lb,'C1')):ax.hist(l[o[d]['cell_type'].astype(str)==t].ravel(),bins=36,density=True,histtype='step',color=col,label=d.upper())
  ax.set_title(t,fontsize=8);ax.set_xlabel('restricted eigenvalue')
 axs[0][0].set_ylabel('density');axs[0][-1].legend(fontsize=8);ins=axs[0][-1].inset_axes([.48,.46,.5,.48]);ins.hist(df.loc[df.kind.str.startswith('refit_'),'spectral_w2'],bins=18,density=True,histtype='step',label='refit');ins.hist(df.loc[df.kind=='pooled','spectral_w2'],bins=18,density=True,histtype='step',label='pooled');ins.axvline(po['spectral_w2'],color='C1',ls='--');ins.axvline(rf['spectral_w2'],color='C0',ls='--');ins.legend(fontsize=5);fig.tight_layout();fig.savefig(OUT/'fig_spectra.png',dpi=180);plt.close(fig)
 x=table[table.cell_type.isin(types)].set_index('cell_type').loc[types];fig,ax=plt.subplots(figsize=(5,4));ax.scatter(x.median_bi_a,x.median_bi_b)
 for t,u,v in zip(types,x.median_bi_a,x.median_bi_b):ax.annotate(t,(u,v),fontsize=8)
 lo,hi=min(x.median_bi_a.min(),x.median_bi_b.min()),max(x.median_bi_a.max(),x.median_bi_b.max());ax.plot([lo,hi],[lo,hi],color='.3',ls='--');ax.set_xlabel('median basin_index A');ax.set_ylabel('median basin_index B');fig.tight_layout();fig.savefig(OUT/'fig_basin.png',dpi=180);plt.close(fig)
 gates_path().write_text(json.dumps(g,sort_keys=True)+'\n');print(json.dumps(g,sort_keys=True))
def report_stage():
 g=json.loads(gates_path().read_text());base='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/K2/scripts/k2_run.py';files=['requirements.lock','match_table.csv','latent_pca30_a.npy','latent_pca30_b.npy','ops_a.npz','ops_b.npz','draws.csv','basin_table.csv','fig_spectra.png','fig_basin.png','scripts/k2_run.py'];regen={**{x:base+' --stage all' for x in files},**{x:base+' --stage match' for x in ('n_cells_a','n_cells_b','count_match_max_abs','depth_match_max_rel')},**{x:base+' --stage latent' for x in ('rank_pca_a','rank_pca_b','frac_rank_ge10')},'symmetry_max':base+' --stage fits',**{x:base+' --stage compare' for x in g if x not in ('n_cells_a','n_cells_b','count_match_max_abs','depth_match_max_rel','rank_pca_a','rank_pca_b','frac_rank_ge10','symmetry_max')}};obj={'task_id':'K2','executor':'manus','branch':'task/K2','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':'2026-10-15T14:33:00-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':g,'deliverables':[{'path':x,'sha256':sha(OUT/x)} for x in files],'regeneration':regen,'deviations':[],'deviations_requested':[],'observations':[f'Input SHA-256: {p.relative_to(ROOT)} = {sha(p)}.' for p in INPUTS],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')};(OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n');print(json.dumps(g,sort_keys=True))
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('match','latent','fits','compare','report','all'));s=p.parse_args().stage
 if s in ('match','all'):match_stage()
 if s in ('latent','all'):latent_stage()
 if s in ('fits','all'):fit_stage()
 if s in ('compare','all'):compare_stage()
 if s in ('report','all'):report_stage()
if __name__=='__main__':main()
