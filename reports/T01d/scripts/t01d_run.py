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
from scipy.stats import spearmanr
from sklearn.decomposition import FactorAnalysis
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP
from scjdo2.opspace import classify_basin_latent

ROOT=Path(__file__).resolve().parents[3]; DATA=ROOT/'data'/'pbmc_v2_v3'; OUT=ROOT/'reports'/'T01d'
SEED,D,K,FLOOR=20261014,30,40,0.05; CLASSES=('', 'holder','fork','drift')
KNOWN=["tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]","tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility"]

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def key(d): return 'v3thin' if d=='v3thin' else 'v2'
def path(d,kind): return OUT/f'{kind}_{key(d)}.'+''
def lat(d): return OUT/f'latent_fa_{key(d)}.npy'
def opsp(d): return OUT/f'ops_{key(d)}.npz'
def driftp(d): return OUT/f'drift_{key(d)}.pt'
def raw(d): return DATA/('pbmc_v3_thin.h5ad' if d=='v3thin' else 'pbmc_v2.h5ad')
def qc(d):
 a=sc.read_h5ad(raw(d)); return a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()
def depth(a): return np.asarray(a.X.sum(axis=1)).ravel().astype(np.int64)
def nz(a): return np.asarray((a.X>0).sum(axis=0)).ravel()

def thin_metrics():
 a2=qc('v2'); a3=sc.read_h5ad(DATA/'pbmc_v3.h5ad'); a3=a3[a3.obs['qc_pass'].to_numpy(dtype=bool)].copy(); at=qc('v3thin')
 if (len(a2),len(a3),len(at))!=(992,1113,1113): raise RuntimeError('QC cell count mismatch')
 table=pd.read_csv(OUT/'depth_table.csv'); tested=table[(table.n_v2>=50)&(table.n_v3>=50)]
 if len(tested)!=4: raise RuntimeError(f'Expected 4 tested types, got {len(tested)}')
 return {'depth_v2_median':float(np.median(depth(a2))),'depth_v3_median':float(np.median(depth(a3))),'depth_v3thin_median':float(np.median(depth(at))),'depth_match_max_rel':float((np.abs(tested.median_depth_v3thin-tested.median_depth_v2)/tested.median_depth_v2).max())}
def thin():
 a2=qc('v2'); a3full=sc.read_h5ad(DATA/'pbmc_v3.h5ad'); a3=a3full[a3full.obs['qc_pass'].to_numpy(dtype=bool)].copy()
 if (len(a2),len(a3))!=(992,1113): raise RuntimeError('Expected 992/1113 QC cells')
 d2,d3=depth(a2),depth(a3); l2=a2.obs['cell_type'].astype(str).to_numpy(); l3=a3.obs['cell_type'].astype(str).to_numpy(); target=np.empty(len(a3))
 for label in sorted(set(l3)):
  ids=np.flatnonzero(l3==label); reference=d2[l2==label]; reference=reference if len(reference) else d2
  ordered=ids[np.lexsort((a3.obs_names[ids].astype(str).to_numpy(),d3[ids]))]; q=(np.arange(len(ordered))+0.5)/len(ordered); target[ordered]=np.quantile(reference,q,method='linear')
 p=np.minimum(1.,target/d3); x=a3.X.tocsr(copy=True) if sparse.issparse(a3.X) else sparse.csr_matrix(a3.X)
 if not np.all(np.equal(x.data,np.floor(x.data))): raise RuntimeError('Input counts are noninteger')
 rng=np.random.default_rng(SEED)
 for c in range(x.shape[0]):
  if p[c]<1:
   s,e=x.indptr[c],x.indptr[c+1]; x.data[s:e]=rng.binomial(x.data[s:e].astype(np.int64),p[c]).astype(x.data.dtype)
 x.eliminate_zeros(); a3.X=x
 if not np.all(np.equal(x.data,np.floor(x.data))): raise RuntimeError('Thinned counts are noninteger')
 a3.write_h5ad(DATA/'pbmc_v3_thin.h5ad'); dt=depth(a3); rows=[]
 for label in sorted(set(l2)|set(l3)):
  m2,m3=l2==label,l3==label; rows.append({'cell_type':label,'n_v2':int(m2.sum()),'n_v3':int(m3.sum()),'median_depth_v2':float(np.median(d2[m2])) if m2.any() else np.nan,'median_depth_v3':float(np.median(d3[m3])) if m3.any() else np.nan,'median_depth_v3thin':float(np.median(dt[m3])) if m3.any() else np.nan})
 pd.DataFrame(rows).to_csv(OUT/'depth_table.csv',index=False); print(json.dumps(thin_metrics(),sort_keys=True))

def select(write=False):
 data={d:qc(d) for d in ('v2','v3thin')}; detected=sorted(set(data['v2'].var_names[nz(data['v2'])>=3].astype(str))&set(data['v3thin'].var_names[nz(data['v3thin'])>=3].astype(str)))
 if len(detected)<2000: raise RuntimeError('Too few detected genes')
 data={d:a[:,detected].copy() for d,a in data.items()}; both=ad.concat([data['v2'],data['v3thin']],axis=0,join='inner',merge='same',index_unique=None)
 sc.pp.highly_variable_genes(both,flavor='seurat_v3',n_top_genes=2000,batch_key=None); genes=sorted(both.var_names[both.var['highly_variable'].to_numpy()].astype(str))
 if len(genes)!=2000: raise RuntimeError(f'Expected 2000 HVGs, got {len(genes)}')
 if write: (OUT/'genes.txt').write_text('\n'.join(genes)+'\n')
 return data,detected,genes
def norm(a):
 a=a.copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
def ranks(z):
 z=np.asarray(z,dtype=np.float32); out=np.zeros(len(z),dtype=np.int32)
 for c in range(len(z)):
  dist=np.linalg.norm(z-z[c],axis=1); dist[c]=np.inf; ids=np.argpartition(dist,K)[:K]; ev=np.abs(np.linalg.eigvalsh(np.cov(z[ids],rowvar=False,bias=False))); out[c]=int((ev>FLOOR*float(ev.max())).sum())
 return out
def latent_metrics():
 genes=(OUT/'genes.txt').read_text().splitlines()
 if len(genes)!=2000 or genes!=sorted(genes): raise RuntimeError('Invalid genes.txt')
 data={d:qc(d) for d in ('v2','v3thin')}; detected=sorted(set(data['v2'].var_names[nz(data['v2'])>=3].astype(str))&set(data['v3thin'].var_names[nz(data['v3thin'])>=3].astype(str)))
 zero={d:int((nz(a[:,genes])==0).sum()) for d,a in data.items()}; r={d:ranks(np.load(lat(d))) for d in data}; allr=np.concatenate([r['v2'],r['v3thin']])
 return {'n_genes_detected_both':len(detected),'n_hvg_zero_v2':zero['v2'],'n_hvg_zero_v3':zero['v3thin'],'rank_fa_v2':float(np.median(r['v2'])),'rank_fa_v3':float(np.median(r['v3thin'])),'frac_fa_at_ceiling':float((allr>=28).mean())}
def latent():
 data,_,genes=select(write=True)
 for d in ('v2','v3thin'): np.save(lat(d),np.asarray(FactorAnalysis(D,random_state=SEED,max_iter=2000).fit_transform(norm(data[d][:,genes].copy())),dtype=np.float32))
 print(json.dumps(latent_metrics(),sort_keys=True))

def cfg(): return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=SEED,batch_size=256,loss='esm',lr=5e-4,beta=1.,esm_laplacian='hutchinson',device='cpu')
def inputs(d):
 z=np.asarray(np.load(lat(d)),dtype=np.float32); labels=qc(d).obs['cell_type'].astype(str).to_numpy()
 if z.shape!=(len(labels),D): raise RuntimeError('Latent input mismatch')
 return z,labels
def save_drift(d,x): torch.save({'state_dict':x._model.state_dict(),'dim':D,'config':asdict(x.config)},driftp(d))
def load_drift(d):
 saved=torch.load(driftp(d),map_location='cpu',weights_only=False); dr=EnergyDriftAdapter(D,config=EnergyDriftConfig(**saved['config'])); model=_EnergyMLP(D,dr.config.hidden,dr.config.depth).to('cpu'); model.load_state_dict(saved['state_dict']); model.eval(); dr._model,dr._device,dr._fitted=model,'cpu',True; return dr
def drift():
 for d in ('v2','v3thin'):
  z,_=inputs(d); save_drift(d,EnergyDriftAdapter(D,config=cfg()).fit(z,verbose=True))
def readops(d):
 with np.load(opsp(d),allow_pickle=False) as x: return {k:x[k] for k in x.files}
def writeops(d,x): np.savez_compressed(opsp(d),**x)
def operators():
 ans={}
 for d,g in (('v2','v2'),('v3thin','v3')):
  z,labels=inputs(d); dr=load_drift(d); J=dr.jacobian_autodiff(z).astype(np.float32); sy=dr.symmetry_check(z,n_cells=100); writeops(d,{'J':J,'eig':np.linalg.eigvalsh(J).astype(np.float32),'cell_type':np.asarray(labels,dtype=str),'symmetry_median':np.asarray(sy['median_asym'],dtype=np.float64),'symmetry_max':np.asarray(sy['max_asym'],dtype=np.float64),'symmetry_mean':np.asarray(sy['mean_asym'],dtype=np.float64)}); ans[f'symmetry_{g}']=float(sy['median_asym'])
 print(json.dumps(ans,sort_keys=True))
def majority(x):
 c=pd.Series(x).value_counts(); winners=c[c==c.max()].index.astype(str).tolist();
 if len(winners)!=1: return ''
 return winners[0]
def basin():
 data={}
 for d,keyd in (('v2','v2'),('v3thin','v3')):
  x=readops(d); z,labels=inputs(d); a=ad.AnnData(X=np.zeros((len(labels),1),dtype=np.float32)); a.obs['cell_type']=pd.Categorical(labels); a.obsm['X_latent']=z; a.uns['scjdo2']={'J_latent':x['J']}; tol=.1*float(np.median(np.abs(x['eig']))); classify_basin_latent(a,tol=tol,use_rep='X_latent',rule='lead_eig'); x['basin_index']=a.obs['op_basin_index'].to_numpy(dtype=np.float32); x['basin_class']=a.obs['op_behaviour_basin'].astype(str).to_numpy(dtype=str); x['tol']=np.asarray(tol,dtype=np.float64); writeops(d,x); data[keyd]=x
 rows=[]
 for t in sorted(set(data['v2']['cell_type'].astype(str))|set(data['v3']['cell_type'].astype(str))):
  row={'cell_type':t}
  for d in ('v2','v3'):
   m=data[d]['cell_type'].astype(str)==t; row[f'n_{d}']=int(m.sum()); row[f'median_bi_{d}']=float(np.median(data[d]['basin_index'][m])) if m.any() else np.nan; row[f'class_{d}']=majority(data[d]['basin_class'][m]) if m.any() else ''
  row['agree']=bool(row['class_v2']==row['class_v3']); rows.append(row)
 pd.DataFrame(rows,columns=['cell_type','n_v2','n_v3','median_bi_v2','median_bi_v3','class_v2','class_v3','agree']).to_csv(OUT/'basin_table.csv',index=False)

def w2(a,b): return float(np.sqrt(ot.wasserstein_1d(a,b,p=2)))
def w1(a,b): return float(ot.wasserstein_1d(a,b,p=1))
def tv(a,b): return float(.5*sum(abs(np.mean(a==c)-np.mean(b==c)) for c in CLASSES))
def quotient(a,b): return 1. if b==0 and a==0 else (float('inf') if b==0 else a/b)
def compare(write=True):
 data={'v2':readops('v2'),'v3':readops('v3thin')}; types=sorted([t for t in set(data['v2']['cell_type'].astype(str))&set(data['v3']['cell_type'].astype(str)) if (data['v2']['cell_type'].astype(str)==t).sum()>=50 and (data['v3']['cell_type'].astype(str)==t).sum()>=50]); ids={d:{t:np.flatnonzero(data[d]['cell_type'].astype(str)==t) for t in types} for d in ('v2','v3')}; m={t:min(len(ids['v2'][t]),len(ids['v3'][t]))//2 for t in types}; rng=np.random.default_rng(SEED); rows=[]
 for draw in range(200):
  a={t:rng.choice(ids['v2'][t],m[t],replace=False) for t in types}; b={t:rng.choice(ids['v3'][t],m[t],replace=False) for t in types}; sp=[]; bi=[]; cl=[]
  for t in types: sp.append(w2(data['v2']['eig'][a[t]].ravel(),data['v3']['eig'][b[t]].ravel())); bi.append(w1(data['v2']['basin_index'][a[t]],data['v3']['basin_index'][b[t]])); cl.append(tv(data['v2']['basin_class'][a[t]],data['v3']['basin_class'][b[t]]))
  rows.append(dict(kind='pooled',draw=draw,spectral_w2=np.mean(sp),basin_w1=np.mean(bi),class_tv=np.mean(cl)))
 for d in ('v2','v3'):
  for draw in range(100):
   sp=[];bi=[];cl=[]
   for t in types:
    c=rng.choice(ids[d][t],2*m[t],replace=False); a,b=c[:m[t]],c[m[t]:]; sp.append(w2(data[d]['eig'][a].ravel(),data[d]['eig'][b].ravel()));bi.append(w1(data[d]['basin_index'][a],data[d]['basin_index'][b]));cl.append(tv(data[d]['basin_class'][a],data[d]['basin_class'][b]))
   rows.append(dict(kind=f'null_{d}',draw=draw,spectral_w2=np.mean(sp),basin_w1=np.mean(bi),class_tv=np.mean(cl)))
 for d in ('v2','v3'):
  for draw in range(100):
   sp=[];bi=[]
   for first,second in itertools.combinations(types,2):
    size=min(m[first],m[second]);a=rng.choice(ids[d][first],size,replace=False);b=rng.choice(ids[d][second],size,replace=False);sp.append(w2(data[d]['eig'][a].ravel(),data[d]['eig'][b].ravel()));bi.append(w1(data[d]['basin_index'][a],data[d]['basin_index'][b]))
   rows.append(dict(kind=f'sep_{d}',draw=draw,spectral_w2=np.mean(sp),basin_w1=np.mean(bi),class_tv=np.nan))
 draws=pd.DataFrame(rows,columns=['kind','draw','spectral_w2','basin_w1','class_tv']);
 if write: draws.to_csv(OUT/'draws.csv',index=False)
 med=lambda mask,col:float(np.median(draws.loc[mask,col])); sp,spn=med(draws.kind=='pooled','spectral_w2'),med(draws.kind.str.startswith('null_'),'spectral_w2'); bi,bin_=med(draws.kind=='pooled','basin_w1'),med(draws.kind.str.startswith('null_'),'basin_w1'); cl,cln=med(draws.kind=='pooled','class_tv'),med(draws.kind.str.startswith('null_'),'class_tv'); table=pd.read_csv(OUT/'basin_table.csv'); lookup=table.set_index('cell_type').loc[types]; rho=float(spearmanr(lookup['median_bi_v2'],lookup['median_bi_v3']).statistic)
 gate={'n_cell_types_tested':len(types),'symmetry_v2':float(data['v2']['symmetry_median']),'symmetry_v3':float(data['v3']['symmetry_median']),'spectral_w2_pooled':sp,'spectral_w2_null':spn,'spectral_ratio':quotient(sp,spn),'basin_w1_pooled':bi,'basin_w1_null':bin_,'basin_ratio':quotient(bi,bin_),'class_tv_pooled':cl,'class_tv_null':cln,'class_tv_ratio':quotient(cl,cln),'spectral_sep_v2':quotient(med(draws.kind=='sep_v2','spectral_w2'),med(draws.kind=='null_v2','spectral_w2')),'spectral_sep_v3':quotient(med(draws.kind=='sep_v3','spectral_w2'),med(draws.kind=='null_v3','spectral_w2')),'basin_sep_v2':quotient(med(draws.kind=='sep_v2','basin_w1'),med(draws.kind=='null_v2','basin_w1')),'basin_sep_v3':quotient(med(draws.kind=='sep_v3','basin_w1'),med(draws.kind=='null_v3','basin_w1')),'basin_index_rho':rho,'basin_class_agreement':float(lookup['agree'].mean())}
 if write:
  fig,axes=plt.subplots(1,len(types),figsize=(3.1*len(types),3.3),squeeze=False);axes=axes[0]
  for ax,t in zip(axes,types):
   for d,c in (('v2','C0'),('v3','C1')): ax.hist(data[d]['eig'][data[d]['cell_type'].astype(str)==t].ravel(),density=True,histtype='step',color=c,label=d)
   ax.set_title(t);ax.set_xlabel('eigenvalue')
  axes[0].set_ylabel('density');axes[-1].legend(title='dataset');ins=axes[-1].inset_axes([.46,.48,.5,.45]);ins.hist(draws.loc[draws.kind.str.startswith('null_'),'spectral_w2'],density=True,color='.7');ins.axvline(sp,color='C3');ins.set_title('null / pooled',fontsize=7);fig.tight_layout();fig.savefig(OUT/'fig_spectra.png',dpi=180);plt.close(fig)
  fig,ax=plt.subplots(figsize=(4.4,4));x,y=lookup['median_bi_v2'].to_numpy(),lookup['median_bi_v3'].to_numpy();ax.scatter(x,y)
  for t,xx,yy in zip(types,x,y):ax.annotate(t,(xx,yy),fontsize=8)
  lo,hi=min(x.min(),y.min()),max(x.max(),y.max());ax.plot([lo,hi],[lo,hi],color='.3',linestyle='--');ax.set_xlabel('median basin_index v2');ax.set_ylabel('median basin_index v3-thin');ax.set_title(f'Spearman rho = {rho:.3f}');fig.tight_layout();fig.savefig(OUT/'fig_basin.png',dpi=180);plt.close(fig)
 print(json.dumps(gate,sort_keys=True)); return gate

def report():
 gate={**thin_metrics(),**latent_metrics(),**compare(True)}; files=['requirements.lock','depth_table.csv','genes.txt','latent_fa_v2.npy','latent_fa_v3thin.npy','ops_v2.npz','ops_v3thin.npz','basin_table.csv','draws.csv','fig_spectra.png','fig_basin.png','scripts/t01d_run.py'];base='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T01d/scripts/t01d_run.py'; obj={'task_id':'T01d','executor':'manus','branch':'task/T01d','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':SEED,'started':'2026-10-15T05:03:03-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':gate,'deliverables':[{'path':x,'sha256':sha(OUT/x)} for x in files],'regeneration':{**{x:base+' --stage all' for x in files},**{x:base+' --stage report' for x in gate}},'deviations':['Tied Plasma cells v3-thin majority (fork=4, drift=4) represented by empty class in basin_table.csv; this type has 8 cells and is outside tested types.'],'deviations_requested':[],'observations':[f"Thinned input pbmc_v3_thin.h5ad SHA-256: {sha(DATA/'pbmc_v3_thin.h5ad')}."],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')};(OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('thin','latent','drift','operators','basin','compare','report','all'));s=p.parse_args().stage
 if s in ('thin','all'):thin()
 if s in ('latent','all'):latent()
 if s in ('drift','all'):drift()
 if s in ('operators','all'):operators()
 if s in ('basin','all'):basin()
 if s in ('compare','all'):compare()
 if s in ('report','all'):report()
if __name__=='__main__':main()
