#!/usr/bin/env python
import argparse, hashlib, json, subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import scanpy as sc
import torch
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'reports'/'K0b'; K0=ROOT/'reports'/'K0'; DATA=ROOT/'data'/'pbmc_10k_v3'
SOURCE=DATA/'pbmc_10k_v3.h5ad'; HALVES=DATA/'k0_halves.csv'; SEEDS=(20261014,20261016,20261017,20261018,20261019); FIT_DRAW_SEED=20261020; R=10
INPUTS=[HALVES,K0/'genes.txt',K0/'latent_pca30_h1.npy',K0/'latent_pca30_h2.npy',K0/'draw_indices.npz',K0/'report.json']
KNOWN=["tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]","tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility"]
CLASSES=('', 'holder','fork','drift')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def latent(h):return K0/f'latent_pca30_{h}.npy'
def ops(h):return OUT/f'ops_avg_{h}.npz'
def drift_path(h,s):return OUT/f'drift_{h}_{s}.pt'
def source_half(h):
 a=sc.read_h5ad(SOURCE);a=a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy();m=pd.read_csv(HALVES).set_index('barcode')['half'];keep=a.obs_names.astype(str).isin(m.index)&(a.obs_names.astype(str).map(m)==h);return a[keep].copy()
def labels(h):
 a=source_half(h);return np.asarray(a.obs['cell_type'].astype(str).to_numpy(),dtype=str),np.asarray(a.obs_names.astype(str).to_numpy(),dtype=str)
def cfg(seed):return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=seed,batch_size=256,loss='esm',lr=5e-4,beta=1.,esm_laplacian='hutchinson',device='cpu')
def tangent_basis(z):
 nbr=NearestNeighbors(n_neighbors=41).fit(z).kneighbors(z,return_distance=False)[:,1:];v=np.empty((len(z),30,R),dtype=np.float32)
 for i in range(len(z)):v[i]=np.linalg.eigh(np.cov(z[nbr[i]],rowvar=False))[1][:,-R:]
 return v
def symmetrize(j):return (j+np.swapaxes(j,1,2))/2
def classify(lam):
 tol=.1*float(np.median(np.abs(lam)));top=lam.max(axis=1);return (-top).astype(np.float64),np.where(top>tol,'fork',np.where(top<-tol,'holder','drift')).astype(str),tol
def save_drift(h,s,dr):torch.save({'state_dict':dr._model.state_dict(),'config':asdict(dr.config)},drift_path(h,s))
def fit_metrics():
 out=[]
 for h in ('H1','H2'):
  z=np.load(latent(h.lower())).astype(np.float32);v=tangent_basis(z).astype(np.float64);jmean=np.zeros((len(z),30,30),dtype=np.float64);lams=[];bis=[];clss=[];tols=[]
  for s in SEEDS:
   dr=EnergyDriftAdapter(30,config=cfg(s)).fit(z,verbose=True);sy=dr.symmetry_check(z,n_cells=100);save_drift(h,s,dr);j=symmetrize(dr.jacobian_autodiff(z).astype(np.float64));jmean+=j;lr=np.einsum('ndi,ndm,nmj->nij',v,j,v);lr=symmetrize(lr);lam=np.linalg.eigvalsh(lr).astype(np.float64);bi,cls,tol=classify(lam);lams.append(lam);bis.append(bi);clss.append(cls);tols.append(tol);out.append({'half':h,'seed':s,'median_asym':float(sy['median_asym']),'max_asym':float(sy['max_asym'])})
  jmean/=len(SEEDS);lamavg=np.linalg.eigvalsh(symmetrize(np.einsum('ndi,ndm,nmj->nij',v,jmean,v))).astype(np.float64);biavg,clsavg,tolavg=classify(lamavg);ct,bc=labels(h);np.savez_compressed(ops(h.lower()),lam_avg=lamavg,lam_seeds=np.stack(lams).astype(np.float64),basin_index_avg=biavg,cls_avg=np.asarray(clsavg,dtype=str),cell_type=ct,barcode=bc,basin_index_seeds=np.stack(bis).astype(np.float64),cls_seeds=np.asarray(clss,dtype=str),tol_seeds=np.asarray(tols,dtype=np.float64),tol_avg=np.asarray(tolavg,dtype=np.float64),symmetry_median=np.asarray([x['median_asym'] for x in out if x['half']==h],dtype=np.float64),seeds=np.asarray(SEEDS,dtype=np.int64))
 pd.DataFrame(out).to_csv(OUT/'fit_symmetry.csv',index=False);print(json.dumps({'symmetry_max':float(max(x['median_asym'] for x in out))},sort_keys=True))
def readops(h):
 with np.load(ops(h),allow_pickle=False) as z:return {k:z[k] for k in z.files}
def w2(a,b):return float(np.sqrt(ot.wasserstein_1d(a,b,p=2)))
def w1(a,b):return float(ot.wasserstein_1d(a,b,p=1))
def tv(a,b):return float(.5*sum(abs(np.mean(a==q)-np.mean(b==q)) for q in CLASSES))
def q(a,b):return 1. if b==0 and a==0 else (float('inf') if b==0 else float(a/b))
def load_plan():
 z=np.load(K0/'draw_indices.npz',allow_pickle=False);types=[str(t) for t in z['types']];m={t:int(x) for t,x in zip(types,z['m'])};plan=[]
 for d in range(200):plan.append({'kind':'pooled','draw':d,'h1':{t:z[f'p_h1_{d}_{i}'] for i,t in enumerate(types)},'h2':{t:z[f'p_h2_{d}_{i}'] for i,t in enumerate(types)}})
 for h in ('h1','h2'):
  for d in range(100):plan.append({'kind':f'null_{h}','draw':d,h:{t:z[f'n_{h}_{d}_{i}'] for i,t in enumerate(types)}})
 for h in ('h1','h2'):
  for d in range(100):
   pairs=[]
   for i,u in enumerate(types):
    for v in types[i+1:]:
     j=len(pairs);pairs.append((u,v))
   plan.append({'kind':f'sep_{h}','draw':d,h:{(u,v):(z[f's_{h}_{d}_{i}_a'],z[f's_{h}_{d}_{i}_b']) for i,(u,v) in enumerate(pairs)}})
 return types,m,plan
def arrays(o,mode,h,seedidx=None):
 if mode=='avg':return o[h]['lam_avg'],o[h]['basin_index_avg'],o[h]['cls_avg']
 return o[h]['lam_seeds'][seedidx],o[h]['basin_index_seeds'][seedidx],o[h]['cls_seeds'][seedidx]
def k0_stats(mode,seedidx=0):
 types,_,plan=load_plan();o={h:readops(h) for h in ('h1','h2')};rows=[]
 for rec in plan:
  kind=rec['kind'];sp=[];bi=[];cl=[]
  if kind=='pooled':
   l1,b1,c1=arrays(o,mode,'h1',seedidx);l2,b2,c2=arrays(o,mode,'h2',seedidx)
   for t in types:
    a,b=rec['h1'][t],rec['h2'][t];sp.append(w2(l1[a].ravel(),l2[b].ravel()));bi.append(w1(b1[a],b2[b]));cl.append(tv(c1[a],c2[b]))
  elif kind.startswith('null_'):
   h=kind[-2:];l,b,c=arrays(o,mode,h,seedidx)
   for t in types:
    x=rec[h][t];a,d=x[:len(x)//2],x[len(x)//2:];sp.append(w2(l[a].ravel(),l[d].ravel()));bi.append(w1(b[a],b[d]));cl.append(tv(c[a],c[d]))
  else:
   h=kind[-2:];l,b,c=arrays(o,mode,h,seedidx)
   for u,v in rec[h]:
    a,d=rec[h][(u,v)];sp.append(w2(l[a].ravel(),l[d].ravel()));bi.append(w1(b[a],b[d]))
  rows.append({'spectra':'avg' if mode=='avg' else 'seed20261014','kind':kind,'draw':rec['draw'],'spectral_w2':float(np.mean(sp)),'basin_w1':float(np.mean(bi)),'class_tv':float(np.mean(cl)) if cl else np.nan})
 return pd.DataFrame(rows)
def summaries(df):
 med=lambda mask,col:float(np.median(df.loc[mask,col]));sp,spn=med(df.kind=='pooled','spectral_w2'),med(df.kind.str.startswith('null_'),'spectral_w2');bi,bin_=med(df.kind=='pooled','basin_w1'),med(df.kind.str.startswith('null_'),'basin_w1');cl,cln=med(df.kind=='pooled','class_tv'),med(df.kind.str.startswith('null_'),'class_tv')
 return {'spectral_ratio':q(sp,spn),'basin_ratio':q(bi,bin_),'class_tv_ratio':q(cl,cln),'spectral_sep_h1':q(med(df.kind=='sep_h1','spectral_w2'),med(df.kind=='null_h1','spectral_w2')),'spectral_sep_h2':q(med(df.kind=='sep_h2','spectral_w2'),med(df.kind=='null_h2','spectral_w2')),'basin_sep_h1':q(med(df.kind=='sep_h1','basin_w1'),med(df.kind=='null_h1','basin_w1')),'basin_sep_h2':q(med(df.kind=='sep_h2','basin_w1'),med(df.kind=='null_h2','basin_w1')),'pooled_spectral':sp,'pooled_basin':bi,'null_spectral_h1':med(df.kind=='null_h1','spectral_w2'),'null_spectral_h2':med(df.kind=='null_h2','spectral_w2'),'null_basin_h1':med(df.kind=='null_h1','basin_w1'),'null_basin_h2':med(df.kind=='null_h2','basin_w1')}
def refit_draws():
 types,m,_=load_plan();o={h:readops(h) for h in ('h1','h2')};rng=np.random.default_rng(FIT_DRAW_SEED);rows=[]
 for h in ('h1','h2'):
  labels_h=o[h]['cell_type'].astype(str);indices={t:np.flatnonzero(labels_h==t) for t in types}
  for d in range(100):
   si=1+d%4;sp=[];bi=[]
   for t in types:
    x=rng.choice(indices[t],2*m[t],replace=False);a,b=x[:m[t]],x[m[t]:];sp.append(w2(o[h]['lam_seeds'][0,a].ravel(),o[h]['lam_seeds'][si,b].ravel()));bi.append(w1(o[h]['basin_index_seeds'][0,a],o[h]['basin_index_seeds'][si,b]))
   rows.append({'half':h.upper(),'draw':d,'seed_a':SEEDS[0],'seed_b':SEEDS[si],'spectral_w2':float(np.mean(sp)),'basin_w1':float(np.mean(bi))})
 out=pd.DataFrame(rows);out.to_csv(OUT/'refit_draws.csv',index=False);return out
def compare():
 seed=k0_stats('seed',0);avg=k0_stats('avg');stats=pd.concat([seed,avg],ignore_index=True);stats.to_csv(OUT/'k0_stats.csv',index=False);re=refit_draws();ss=summaries(seed);sa=summaries(avg);old=json.loads((K0/'report.json').read_text())['gated'];rels=[abs(ss[k]-old[f'{k}_pca30'])/abs(old[f'{k}_pca30']) for k in ('spectral_ratio','basin_ratio','class_tv_ratio')];g={'k0_repro_max_rel':float(max(rels)),'refit_w2_h1':float(np.median(re.loc[re.half=='H1','spectral_w2'])),'refit_w2_h2':float(np.median(re.loc[re.half=='H2','spectral_w2'])),'fitnoise_spectral_h1':q(float(np.median(re.loc[re.half=='H1','spectral_w2'])),ss['null_spectral_h1']),'fitnoise_spectral_h2':q(float(np.median(re.loc[re.half=='H2','spectral_w2'])),ss['null_spectral_h2']),'fitnoise_basin_h1':q(float(np.median(re.loc[re.half=='H1','basin_w1'])),ss['null_basin_h1']),'fitnoise_basin_h2':q(float(np.median(re.loc[re.half=='H2','basin_w1'])),ss['null_basin_h2']),'k0_vs_refit_spectral':q(ss['pooled_spectral'],float(np.median(re.spectral_w2))),'k0_vs_refit_basin':q(ss['pooled_basin'],float(np.median(re.basin_w1))),'spectral_ratio_avg':sa['spectral_ratio'],'basin_ratio_avg':sa['basin_ratio'],'class_tv_ratio_avg':sa['class_tv_ratio'],'spectral_sep_h1_avg':sa['spectral_sep_h1'],'spectral_sep_h2_avg':sa['spectral_sep_h2'],'basin_sep_h1_avg':sa['basin_sep_h1'],'basin_sep_h2_avg':sa['basin_sep_h2']}
 if min(g['spectral_sep_h1_avg'],g['spectral_sep_h2_avg'],g['basin_sep_h1_avg'],g['basin_sep_h2_avg'])<2:g['outcome_avg']='inconclusive'
 elif g['spectral_ratio_avg']>3 or g['basin_ratio_avg']>3:g['outcome_avg']='fails'
 elif g['spectral_ratio_avg']<=1.5 and g['basin_ratio_avg']<=1.5 and g['class_tv_ratio_avg']<=1.5:g['outcome_avg']='survives'
 else:g['outcome_avg']='inconclusive'
 # Fit-noise figure
 fig,axes=plt.subplots(1,2,figsize=(10,3.8),sharey=True)
 for ax,h in zip(axes,('h1','h2')):
  ax.hist(seed.loc[seed.kind==f'null_{h}','spectral_w2'],bins=22,density=True,histtype='step',label='K0 null');ax.hist(re.loc[re.half==h.upper(),'spectral_w2'],bins=22,density=True,histtype='step',label='refit');ax.hist(seed.loc[seed.kind=='pooled','spectral_w2'],bins=22,density=True,histtype='step',label='K0 pooled')
  for val,col in [(np.median(seed.loc[seed.kind==f'null_{h}','spectral_w2']),'C0'),(np.median(re.loc[re.half==h.upper(),'spectral_w2']),'C1'),(np.median(seed.loc[seed.kind=='pooled','spectral_w2']),'C2')]:ax.axvline(val,color=col,ls='--')
  ax.set_title(h.upper());ax.set_xlabel('spectral W2');ax.legend(fontsize=8)
 axes[0].set_ylabel('density');fig.tight_layout();fig.savefig(OUT/'fig_fitnoise.png',dpi=180);plt.close(fig)
 types,_,_=load_plan();o={h:readops(h) for h in ('h1','h2')};fig,axs=plt.subplots(1,len(types),figsize=(3*len(types),3.1),squeeze=False)
 for ax,t in zip(axs[0],types):
  for h,col in (('h1','C0'),('h2','C1')):ax.hist(o[h]['lam_avg'][o[h]['cell_type'].astype(str)==t].ravel(),bins=36,density=True,histtype='step',color=col,label=h)
  ax.set_title(t,fontsize=8);ax.set_xlabel('averaged restricted eigenvalue')
 axs[0][0].set_ylabel('density');axs[0][-1].legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'fig_spectra_avg.png',dpi=180);plt.close(fig)
 print(json.dumps(g,sort_keys=True));return g
def fits_gates():
 a=pd.read_csv(OUT/'fit_symmetry.csv');return {'symmetry_max':float(a.median_asym.max())}
def report():
 gated={**fits_gates(),**compare()};files=['requirements.lock','ops_avg_h1.npz','ops_avg_h2.npz','refit_draws.csv','k0_stats.csv','fig_fitnoise.png','fig_spectra_avg.png','scripts/k0b_run.py'];base='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/K0b/scripts/k0b_run.py';status='complete';stop=None
 if gated['symmetry_max']>1e-4:status='stopped';stop='0. Validity: symmetry_max > 1e-4 → REVISE'
 elif gated['k0_repro_max_rel']>1e-4:status='stopped';stop='0. Validity: k0_repro_max_rel > 1e-4 → REVISE (seed-20261014 fits do not reproduce K0)'
 obs=[f'Input SHA-256: {p.relative_to(ROOT)} = {sha(p)}.' for p in INPUTS]
 obj={'task_id':'K0b','executor':'manus','branch':'task/K0b','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':'2026-10-15T11:06:20-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':status,'gated':gated,'deliverables':[{'path':x,'sha256':sha(OUT/x)} for x in files],'regeneration':{**{x:base+' --stage all' for x in files},**{x:base+' --stage report' for x in gated}},'deviations':[],'deviations_requested':[],'observations':obs,'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')}
 if stop:obj['stop_rule']=stop
 (OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('fits','compare','report','all'));s=p.parse_args().stage
 if s in ('fits','all'):fit_metrics()
 if s in ('compare','all'):compare()
 if s in ('report','all'):report()
if __name__=='__main__':main()
