#!/usr/bin/env python
import argparse, hashlib, itertools, json, subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import scanpy as sc
import torch
from scipy import sparse
from sklearn.decomposition import FactorAnalysis, PCA
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'data'/'pbmc_10k_v3';OUT=ROOT/'reports'/'K0';SOURCE=DATA/'pbmc_10k_v3.h5ad';HALVES=DATA/'k0_halves.csv'
SEED=20261014;DRAW_SEED=20261015;K=40;R=10;FLOOR=.05;CLASSES=('', 'holder','fork','drift');CONFIGS={'pca30':30,'fa30':30,'fa10':10}
KNOWN=["tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]","tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility"]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def nz(a):return np.asarray((a.X>0).sum(axis=0)).ravel()
def lat(c,h):return OUT/f'latent_{c}_{h}.npy'
def rankp(c,h):return OUT/f'rank_{c}_{h}.npy'
def basp(c,h):return OUT/f'basis_{c}_{h}.npy'
def drp(c,h):return OUT/f'drift_{c}_{h}.pt'
def opp(c,h):return OUT/f'ops_{c}_{h}.npz'
def raw_qc():
 a=sc.read_h5ad(SOURCE);return a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()
def half_data():
 a=raw_qc();csv=pd.read_csv(HALVES);assigned=csv.set_index('barcode')['half'];m1=a.obs_names.astype(str).isin(assigned.index)&(a.obs_names.astype(str).map(assigned)=='H1');m2=a.obs_names.astype(str).isin(assigned.index)&(a.obs_names.astype(str).map(assigned)=='H2');return a[m1].copy(),a[m2].copy()
def split_metrics():
 a=raw_qc();h1,h2=half_data();t=pd.read_csv(OUT/'split_table.csv');return {'n_cells_qc':int(a.n_obs),'n_cells_h1':int(h1.n_obs),'n_cells_h2':int(h2.n_obs),'count_match_max_abs':int(np.abs(t.n_h1-t.n_h2).max())}
def split():
 a=raw_qc();labels=a.obs['cell_type'].astype(str).to_numpy();rng=np.random.default_rng(SEED);halves=np.full(len(a),'',dtype=object);rows=[]
 for t in sorted(set(labels)):
  ids=np.flatnonzero(labels==t);perm=rng.permutation(len(ids));h=len(ids)//2;halves[ids[perm[:h]]]='H1';halves[ids[perm[h:2*h]]]='H2';rows.append({'cell_type':t,'n_qc':len(ids),'n_h1':h,'n_h2':h})
 keep=halves!='';pd.DataFrame({'barcode':a.obs_names.astype(str)[keep],'half':halves[keep]}).to_csv(HALVES,index=False);pd.DataFrame(rows,columns=['cell_type','n_qc','n_h1','n_h2']).to_csv(OUT/'split_table.csv',index=False);print(json.dumps(split_metrics(),sort_keys=True))
def select(write=False):
 h1,h2=half_data();detected=sorted(set(h1.var_names[nz(h1)>=3].astype(str))&set(h2.var_names[nz(h2)>=3].astype(str)))
 if len(detected)<2000:raise RuntimeError('Too few detected genes')
 h1,h2=h1[:,detected].copy(),h2[:,detected].copy();both=ad.concat([h1,h2],axis=0,join='inner',merge='same',index_unique=None)
 if list(both.var_names.astype(str))!=sorted(both.var_names.astype(str)):raise RuntimeError('Genes not sorted alphabetically before HVG selection')
 sc.pp.highly_variable_genes(both,flavor='seurat_v3',n_top_genes=2000,batch_key=None);genes=sorted(both.var_names[both.var['highly_variable'].to_numpy()].astype(str))
 if len(genes)!=2000:raise RuntimeError('Expected 2000 HVGs')
 if write:(OUT/'genes.txt').write_text('\n'.join(genes)+'\n')
 return {'h1':h1,'h2':h2},detected,genes
def norm(a):
 a=a.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
def tangent(z,need_basis):
 nbr=NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z,return_distance=False)[:,1:];rank=np.zeros(len(z),dtype=np.int32);basis=np.empty((len(z),z.shape[1],R),dtype=np.float32) if need_basis else None
 for i in range(len(z)):
  w,u=np.linalg.eigh(np.cov(z[nbr[i]],rowvar=False));rank[i]=int((w>FLOOR*float(w.max())).sum())
  if need_basis:basis[i]=u[:,-R:]
 return rank,basis
def fit(c,x):
 if c=='pca30':return PCA(n_components=30,svd_solver='full',random_state=SEED).fit_transform(x)
 return FactorAnalysis(CONFIGS[c],random_state=SEED,max_iter=2000).fit_transform(x)
def latent_metrics():
 genes=(OUT/'genes.txt').read_text().splitlines();h1,h2=half_data();detected=sorted(set(h1.var_names[nz(h1)>=3].astype(str))&set(h2.var_names[nz(h2)>=3].astype(str)));out={'n_genes_detected_both':len(detected),'n_hvg_zero_h1':int((nz(h1[:,genes])==0).sum()),'n_hvg_zero_h2':int((nz(h2[:,genes])==0).sum())}
 for c in CONFIGS:
  r1,r2=np.load(rankp(c,'h1')),np.load(rankp(c,'h2'));out[f'rank_h1_{c}']=float(np.median(r1));out[f'rank_h2_{c}']=float(np.median(r2));out[f'frac_rank_ge10_{c}']=float((np.concatenate([r1,r2])>=10).mean())
 return out
def latent():
 data,_,genes=select(True)
 for c in CONFIGS:
  for h in ('h1','h2'):
   z=np.asarray(fit(c,norm(data[h][:,genes].copy())),dtype=np.float32);np.save(lat(c,h),z);r,b=tangent(z,c!='fa10');np.save(rankp(c,h),r)
   if b is not None:np.save(basp(c,h),b)
 print(json.dumps(latent_metrics(),sort_keys=True))
def cfg(dim):return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=SEED,batch_size=256,loss='esm',lr=5e-4,beta=1.,esm_laplacian='hutchinson',device='cpu')
def save_dr(c,h,dr):torch.save({'state_dict':dr._model.state_dict(),'dim':CONFIGS[c],'config':asdict(dr.config)},drp(c,h))
def load_dr(c,h):
 x=torch.load(drp(c,h),map_location='cpu',weights_only=False);dr=EnergyDriftAdapter(int(x['dim']),config=EnergyDriftConfig(**x['config']));m=_EnergyMLP(int(x['dim']),dr.config.hidden,dr.config.depth).to('cpu');m.load_state_dict(x['state_dict']);m.eval();dr._model,dr._device,dr._fitted=m,'cpu',True;return dr
def drift():
 for c,d in CONFIGS.items():
  for h in ('h1','h2'):save_dr(c,h,EnergyDriftAdapter(d,config=cfg(d)).fit(np.load(lat(c,h)).astype(np.float32),verbose=True))
def labels(h):
 a={'h1':half_data()[0],'h2':half_data()[1]}[h];return a.obs['cell_type'].astype(str).to_numpy(),a.obs_names.astype(str).to_numpy()
def readops(c,h):
 with np.load(opp(c,h),allow_pickle=False) as x:return {k:x[k] for k in x.files}
def operators():
 gates={}
 for c in CONFIGS:
  for h in ('h1','h2'):
   z=np.load(lat(c,h)).astype(np.float32);dr=load_dr(c,h);J=dr.jacobian_autodiff(z).astype(np.float32);sy=dr.symmetry_check(z,n_cells=100);sym=(J+np.swapaxes(J,1,2))/2
   if c=='fa10':lam=np.linalg.eigvalsh(sym).astype(np.float32);full=np.empty((len(z),0),dtype=np.float32)
   else:
    V=np.load(basp(c,h));jr=np.einsum('ndi,ndm,nmj->nij',V,sym,V);jr=(jr+np.swapaxes(jr,1,2))/2;lam=np.linalg.eigvalsh(jr).astype(np.float32);full=np.linalg.eigvalsh(sym).astype(np.float32)
   top=lam.max(axis=1);tol=.1*float(np.median(np.abs(lam)));cls=np.where(top>tol,'fork',np.where(top<-tol,'holder','drift')).astype(str);lab,bar=labels(h);lab=np.asarray(lab,dtype=str);bar=np.asarray(bar,dtype=str);np.savez_compressed(opp(c,h),lam=lam,lam_full=full,basin_index=(-top).astype(np.float32),cls=cls,cell_type=lab,barcode=bar,symmetry_median=np.asarray(sy['median_asym']),symmetry_max=np.asarray(sy['max_asym']),tol=np.asarray(tol));gates[f'symmetry_{h}_{c}']=float(sy['median_asym'])
 print(json.dumps(gates,sort_keys=True))
def w2(a,b):return float(np.sqrt(ot.wasserstein_1d(a,b,p=2)))
def w1(a,b):return float(ot.wasserstein_1d(a,b,p=1))
def tv(a,b):return float(.5*sum(abs(np.mean(a==q)-np.mean(b==q)) for q in CLASSES))
def quotient(a,b):return 1. if b==0 and a==0 else (float('inf') if b==0 else a/b)
def majority(a):
 x=pd.Series(a).value_counts();w=x[x==x.max()].index.astype(str).tolist();return w[0] if len(w)==1 else ''
def types_indices():
 a,b=readops('pca30','h1'),readops('pca30','h2');types=sorted([t for t in set(a['cell_type'].astype(str))&set(b['cell_type'].astype(str)) if (a['cell_type'].astype(str)==t).sum()>=50 and (b['cell_type'].astype(str)==t).sum()>=50]);idx={'h1':{t:np.flatnonzero(a['cell_type'].astype(str)==t) for t in types},'h2':{t:np.flatnonzero(b['cell_type'].astype(str)==t) for t in types}};m={t:min(len(idx['h1'][t]),len(idx['h2'][t]))//2 for t in types};return types,idx,m
def draws_once():
 types,idx,m=types_indices();rng=np.random.default_rng(DRAW_SEED);pack={'types':np.array(types,dtype=str),'m':np.array([m[t] for t in types],dtype=np.int32)};plan=[]
 for d in range(200):
  rec={'kind':'pooled','draw':d,'h1':{},'h2':{}}
  for i,t in enumerate(types):rec['h1'][t]=rng.choice(idx['h1'][t],m[t],replace=False);pack[f'p_h1_{d}_{i}']=rec['h1'][t]
  for i,t in enumerate(types):rec['h2'][t]=rng.choice(idx['h2'][t],m[t],replace=False);pack[f'p_h2_{d}_{i}']=rec['h2'][t]
  plan.append(rec)
 for h in ('h1','h2'):
  for d in range(100):
   rec={'kind':f'null_{h}','draw':d,h:{}}
   for i,t in enumerate(types):rec[h][t]=rng.choice(idx[h][t],2*m[t],replace=False);pack[f'n_{h}_{d}_{i}']=rec[h][t]
   plan.append(rec)
 for h in ('h1','h2'):
  for d in range(100):
   rec={'kind':f'sep_{h}','draw':d,h:{}}
   for i,(u,v) in enumerate(itertools.combinations(types,2)):
    size=min(m[u],m[v]);rec[h][(u,v)]=(rng.choice(idx[h][u],size,replace=False),rng.choice(idx[h][v],size,replace=False));pack[f's_{h}_{d}_{i}_a']=rec[h][(u,v)][0];pack[f's_{h}_{d}_{i}_b']=rec[h][(u,v)][1]
   plan.append(rec)
 np.savez_compressed(OUT/'draw_indices.npz',**pack);return types,plan
def metrics_config(c,types,plan):
 o={h:readops(c,h) for h in ('h1','h2')};rows=[]
 for rec in plan:
  kind=rec['kind'];sp=[];bi=[];cl=[];sf=[]
  if kind=='pooled':
   for t in types:
    a,b=rec['h1'][t],rec['h2'][t];sp.append(w2(o['h1']['lam'][a].ravel(),o['h2']['lam'][b].ravel()));bi.append(w1(o['h1']['basin_index'][a],o['h2']['basin_index'][b]));cl.append(tv(o['h1']['cls'][a],o['h2']['cls'][b]));
    if c!='fa10':sf.append(w2(o['h1']['lam_full'][a].ravel(),o['h2']['lam_full'][b].ravel()))
  elif kind.startswith('null_'):
   h=kind[-2:]
   for t in types:
    q=rec[h][t];a,b=q[:len(q)//2],q[len(q)//2:];sp.append(w2(o[h]['lam'][a].ravel(),o[h]['lam'][b].ravel()));bi.append(w1(o[h]['basin_index'][a],o[h]['basin_index'][b]));cl.append(tv(o[h]['cls'][a],o[h]['cls'][b]));
    if c!='fa10':sf.append(w2(o[h]['lam_full'][a].ravel(),o[h]['lam_full'][b].ravel()))
  else:
   h=kind[-2:]
   for u,v in itertools.combinations(types,2):
    a,b=rec[h][(u,v)];sp.append(w2(o[h]['lam'][a].ravel(),o[h]['lam'][b].ravel()));bi.append(w1(o[h]['basin_index'][a],o[h]['basin_index'][b]));
    if c!='fa10':sf.append(w2(o[h]['lam_full'][a].ravel(),o[h]['lam_full'][b].ravel()))
  rows.append({'config':c,'kind':kind,'draw':rec['draw'],'spectral_w2':float(np.mean(sp)),'basin_w1':float(np.mean(bi)),'class_tv':float(np.mean(cl)) if cl else np.nan,'spectral_w2_full':float(np.mean(sf)) if sf else np.nan})
 df=pd.DataFrame(rows);med=lambda mask,col:float(np.median(df.loc[mask,col]));sp,spn=med(df.kind=='pooled','spectral_w2'),med(df.kind.str.startswith('null_'),'spectral_w2');bi,bin_=med(df.kind=='pooled','basin_w1'),med(df.kind.str.startswith('null_'),'basin_w1');cl,cln=med(df.kind=='pooled','class_tv'),med(df.kind.str.startswith('null_'),'class_tv');g={f'spectral_ratio_{c}':quotient(sp,spn),f'basin_ratio_{c}':quotient(bi,bin_),f'class_tv_ratio_{c}':quotient(cl,cln),f'spectral_sep_h1_{c}':quotient(med(df.kind=='sep_h1','spectral_w2'),med(df.kind=='null_h1','spectral_w2')),f'spectral_sep_h2_{c}':quotient(med(df.kind=='sep_h2','spectral_w2'),med(df.kind=='null_h2','spectral_w2')),f'basin_sep_h1_{c}':quotient(med(df.kind=='sep_h1','basin_w1'),med(df.kind=='null_h1','basin_w1')),f'basin_sep_h2_{c}':quotient(med(df.kind=='sep_h2','basin_w1'),med(df.kind=='null_h2','basin_w1'))}
 if c!='fa10':g[f'spectral_ratio_full_{c}']=quotient(med(df.kind=='pooled','spectral_w2_full'),med(df.kind.str.startswith('null_'),'spectral_w2_full'))
 return df,g
def compare(write=True):
 types,plan=draws_once();allrows=[];gates={'n_cell_types_tested':len(types)};gates.update(latent_metrics());bas=[]
 for c in CONFIGS:
  df,g=metrics_config(c,types,plan);allrows.append(df);gates.update(g);o={h:readops(c,h) for h in ('h1','h2')};gates[f'symmetry_h1_{c}']=float(o['h1']['symmetry_median']);gates[f'symmetry_h2_{c}']=float(o['h2']['symmetry_median'])
  for t in sorted(set(o['h1']['cell_type'].astype(str))|set(o['h2']['cell_type'].astype(str))):
   row={'config':c,'cell_type':t}
   for h in ('h1','h2'):
    mask=o[h]['cell_type'].astype(str)==t;row[f'n_{h}']=int(mask.sum());row[f'median_bi_{h}']=float(np.median(o[h]['basin_index'][mask])) if mask.any() else np.nan;row[f'class_{h}']=majority(o[h]['cls'][mask]) if mask.any() else ''
   row['agree']=bool(row['class_h1']==row['class_h2']);bas.append(row)
 draws=pd.concat(allrows,ignore_index=True);pd.DataFrame(bas,columns=['config','cell_type','n_h1','n_h2','median_bi_h1','median_bi_h2','class_h1','class_h2','agree']).to_csv(OUT/'basin_table.csv',index=False);draws.to_csv(OUT/'draws.csv',index=False)
 # figures: rows=configs, columns=tested types
 fig,axes=plt.subplots(len(CONFIGS),len(types),figsize=(3.0*len(types),2.8*len(CONFIGS)),squeeze=False)
 for ri,c in enumerate(CONFIGS):
  o={h:readops(c,h) for h in ('h1','h2')}
  for ci,t in enumerate(types):
   ax=axes[ri,ci]
   for h,col in (('h1','C0'),('h2','C1')):ax.hist(o[h]['lam'][o[h]['cell_type'].astype(str)==t].ravel(),density=True,histtype='step',color=col,label=h)
   ax.set_title(f'{c}: {t}',fontsize=8);ax.set_xlabel('tangent eigenvalue')
  axes[ri,0].set_ylabel('density');axes[ri,-1].legend(fontsize=7)
  d=draws[draws.config==c];ins=axes[ri,-1].inset_axes([.48,.47,.48,.46]);ins.hist(d.loc[d.kind.str.startswith('null_'),'spectral_w2'],density=True,color='.7');ins.axvline(np.median(d.loc[d.kind=='pooled','spectral_w2']),color='C3');ins.set_title('null / pooled',fontsize=6)
 fig.tight_layout();fig.savefig(OUT/'fig_spectra.png',dpi=180);plt.close(fig)
 fig,axes=plt.subplots(1,len(CONFIGS),figsize=(4.2*len(CONFIGS),3.8),squeeze=False);table=pd.DataFrame(bas)
 for ax,c in zip(axes[0],CONFIGS):
  x=table[(table.config==c)&table.cell_type.isin(types)].set_index('cell_type').loc[types];xx,yy=x.median_bi_h1.to_numpy(),x.median_bi_h2.to_numpy();ax.scatter(xx,yy)
  for t,u,v in zip(types,xx,yy):ax.annotate(t,(u,v),fontsize=7)
  lo,hi=min(xx.min(),yy.min()),max(xx.max(),yy.max());ax.plot([lo,hi],[lo,hi],color='.3',ls='--');ax.set_title(c);ax.set_xlabel('median basin_index H1');ax.set_ylabel('median basin_index H2')
 fig.tight_layout();fig.savefig(OUT/'fig_basin.png',dpi=180);plt.close(fig)
 # selection after all config metrics
 outcomes=[]
 for c in CONFIGS:
  if gates[f'symmetry_h1_{c}']>1e-4 or gates[f'symmetry_h2_{c}']>1e-4:out='revise'
  elif c in ('pca30','fa30') and gates[f'frac_rank_ge10_{c}']<.80:out='excluded'
  elif min(gates[f'spectral_sep_h1_{c}'],gates[f'spectral_sep_h2_{c}'],gates[f'basin_sep_h1_{c}'],gates[f'basin_sep_h2_{c}'])<2:out='inconclusive'
  elif gates[f'spectral_ratio_{c}']>3 or gates[f'basin_ratio_{c}']>3:out='fails'
  elif gates[f'spectral_ratio_{c}']<=1.5 and gates[f'basin_ratio_{c}']<=1.5 and gates[f'class_tv_ratio_{c}']<=1.5:out='survives'
  else:out='inconclusive'
  outcomes.append(out)
 gates['selected_layer']=next((c for c,o in zip(CONFIGS,outcomes) if o=='survives'),'none');print(json.dumps(gates,sort_keys=True));return gates
def report():
 gated={**split_metrics(),**latent_metrics(),**compare(True)};files=['requirements.lock','split_table.csv','genes.txt','latent_pca30_h1.npy','latent_pca30_h2.npy','latent_fa30_h1.npy','latent_fa30_h2.npy','latent_fa10_h1.npy','latent_fa10_h2.npy','ops_pca30_h1.npz','ops_pca30_h2.npz','ops_fa30_h1.npz','ops_fa30_h2.npz','ops_fa10_h1.npz','ops_fa10_h2.npz','draw_indices.npz','basin_table.csv','draws.csv','fig_spectra.png','fig_basin.png','scripts/k0_run.py'];base='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/K0/scripts/k0_run.py';obj={'task_id':'K0','executor':'manus','branch':'task/K0','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':SEED,'started':'2026-10-15T08:26:08-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':gated,'deliverables':[{'path':x,'sha256':sha(OUT/x)} for x in files],'regeneration':{**{x:base+' --stage all' for x in files},**{x:base+' --stage report' for x in gated}},'deviations':[],'deviations_requested':[],'observations':[f'Input pbmc_10k_v3.h5ad SHA-256: {sha(SOURCE)}.'],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')};(OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('split','latent','drift','operators','compare','report','all'));s=p.parse_args().stage
 if s in ('split','all'):split()
 if s in ('latent','all'):latent()
 if s in ('drift','all'):drift()
 if s in ('operators','all'):operators()
 if s in ('compare','all'):compare()
 if s in ('report','all'):report()
if __name__=='__main__':main()
