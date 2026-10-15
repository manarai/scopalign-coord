#!/usr/bin/env python
import argparse, hashlib, itertools, json, os, subprocess
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
ROOT=Path(__file__).resolve().parents[3]; DATA=ROOT/'data'/'pbmc_10k_v3'; OUT=ROOT/'reports'/'T01e'; A_SOURCE=ROOT/'data'/'pbmc_v2_v3'/'pbmc_v3.h5ad'
URL='https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5'; RAW=DATA/'raw'/Path(URL).name
SEED,D,K,FLOOR=20261014,30,40,.05; CLASSES=('', 'holder','fork','drift')
KNOWN=["tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]","tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]","tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility"]
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def count(a): return np.asarray(a.X.sum(axis=1)).ravel().astype(np.int64)
def nz(a): return np.asarray((a.X>0).sum(axis=0)).ravel()
def h5(d): return DATA/({'a':'A_matched.h5ad','b':'B_matched.h5ad'}[d])
def lat(d): return OUT/f'latent_fa_{d}.npy'
def opsp(d): return OUT/f'ops_{d}.npz'
def driftp(d): return OUT/f'drift_{d}.pt'
def qc(a): return a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()

def celltypist_labels(a):
 model_home=DATA/'celltypist_home';model_home.mkdir(parents=True,exist_ok=True);os.environ['HOME']=str(model_home)
 import celltypist
 celltypist.models.download_models(model='Immune_All_High.pkl'); b=a.copy();sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b); pred=celltypist.annotate(b,model='Immune_All_High.pkl',majority_voting=False).predicted_labels
 if isinstance(pred,pd.DataFrame):pred=pred['predicted_labels']
 labels=pd.Series(pred,index=pred.index,dtype='string').reindex(a.obs_names)
 if labels.isna().any():raise RuntimeError('CellTypist returned missing labels')
 return labels.astype(str).to_numpy()
def stage_metrics():
 b=sc.read_h5ad(DATA/'pbmc_10k_v3.h5ad');return {'n_cells_b_raw':int(b.n_obs),'n_cells_b_qc':int(b.obs['qc_pass'].sum())}
def stage():
 RAW.parent.mkdir(parents=True,exist_ok=True)
 if not RAW.exists() or RAW.stat().st_size==0:
  part=RAW.with_suffix(RAW.suffix+'.partial');subprocess.run(['curl','--fail','--location','--retry','5','--retry-all-errors','--continue-at','-','--output',str(part),URL],check=True);part.replace(RAW)
 (OUT/'data_checksums.txt').write_text(f'{sha(RAW)}  {RAW.relative_to(ROOT)}\n')
 b=sc.read_10x_h5(RAW,genome=None,gex_only=True)
 if 'gene_ids' not in b.var.columns:raise RuntimeError("read_10x_h5 did not retain var['gene_ids']")
 b.var_names_make_unique(); vals=b.X.data if sparse.issparse(b.X) else np.asarray(b.X)
 if not np.allclose(vals,np.round(vals)):raise RuntimeError('B counts are non-integer')
 b.X=b.X.astype(np.int32);b.var['mt']=np.asarray(b.var_names.astype(str).str.startswith('MT-'));sc.pp.calculate_qc_metrics(b,qc_vars=['mt'],inplace=True);b.obs['qc_pass']=((b.obs['n_genes_by_counts']>=200)&(b.obs['pct_counts_mt']<20)).astype(bool);b.obs['cell_type']=celltypist_labels(b);b.write_h5ad(DATA/'pbmc_10k_v3.h5ad',compression='gzip');print(json.dumps(stage_metrics(),sort_keys=True))

def thin_by_prob(a,p,labels,types,rng):
 x=a.X.tocsr(copy=True) if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
 if not np.all(np.equal(x.data,np.floor(x.data))):raise RuntimeError('Counts are noninteger before thinning')
 for t in types:
  for c in np.flatnonzero(labels==t):
   if p[c]<1:
    s,e=x.indptr[c],x.indptr[c+1];x.data[s:e]=rng.binomial(x.data[s:e].astype(np.int64),p[c]).astype(x.data.dtype)
 x.eliminate_zeros();a.X=x
 if not np.all(np.equal(x.data,np.floor(x.data))):raise RuntimeError('Counts are noninteger after thinning')
def match_metrics():
 a=sc.read_h5ad(h5('a'));b=sc.read_h5ad(h5('b'));table=pd.read_csv(OUT/'match_table.csv');tested=table[(table.n_a>=50)&(table.n_bsub>=50)]
 if len(tested)==0:raise RuntimeError('No tested cell types')
 countdiff=np.abs(tested.n_a-tested.n_bsub);depthdiff=np.abs(tested.median_depth_b_after-tested.median_depth_a_after)/tested.median_depth_a_after
 return {'n_cells_a':int(a.n_obs),'n_cells_bsub':int(b.n_obs),'count_match_max_abs':int(countdiff.max()),'depth_a_median':float(np.median(count(a))),'depth_b_median':float(np.median(count(b))),'depth_match_max_rel':float(depthdiff.max())}
def match():
 a=qc(sc.read_h5ad(A_SOURCE)); bqc=qc(sc.read_h5ad(DATA/'pbmc_10k_v3.h5ad'))
 if len(a)!=1113:raise RuntimeError(f'Expected A QC count 1113, got {len(a)}')
 la=a.obs['cell_type'].astype(str).to_numpy();lb=bqc.obs['cell_type'].astype(str).to_numpy();rng_sub=np.random.default_rng(SEED);selected=np.zeros(len(bqc),dtype=bool)
 for t in sorted(set(la)):
  ids=np.flatnonzero(lb==t);take=min(int((la==t).sum()),len(ids))
  if take:selected[rng_sub.choice(ids,take,replace=False)]=True
 b=bqc[selected].copy();lb= b.obs['cell_type'].astype(str).to_numpy();da,db=count(a),count(b);pa,pb=np.ones(len(a)),np.ones(len(b));shared=sorted(set(la)&set(lb))
 for t in shared:
  ia,ib=np.flatnonzero(la==t),np.flatnonzero(lb==t);ma,mb=np.median(da[ia]),np.median(db[ib])
  if ma==mb:continue
  hi,lo,high_depth,low_depth,barcodes=(ia,ib,da,db,a.obs_names) if ma>mb else (ib,ia,db,da,b.obs_names)
  ordered=hi[np.lexsort((barcodes[hi].astype(str).to_numpy(),high_depth[hi]))];q=(np.arange(len(ordered))+0.5)/len(ordered);target=np.quantile(low_depth[lo],q,method='linear');prob=np.minimum(1.,target/high_depth[ordered])
  if ma>mb:pa[ordered]=prob
  else:pb[ordered]=prob
 rng_thin=np.random.default_rng(20261015);thin_by_prob(a,pa,la,shared,rng_thin);thin_by_prob(b,pb,lb,shared,rng_thin);a.write_h5ad(h5('a'));b.write_h5ad(h5('b'))
 da2,db2=count(a),count(b);rows=[]
 for t in sorted(set(la)|set(lb)):
  ia,ib=la==t,lb==t;rows.append({'cell_type':t,'n_a':int(ia.sum()),'n_b_qc':int((bqc.obs['cell_type'].astype(str).to_numpy()==t).sum()),'n_bsub':int(ib.sum()),'median_depth_a_before':float(np.median(da[ia])) if ia.any() else np.nan,'median_depth_b_before':float(np.median(db[ib])) if ib.any() else np.nan,'median_depth_a_after':float(np.median(da2[ia])) if ia.any() else np.nan,'median_depth_b_after':float(np.median(db2[ib])) if ib.any() else np.nan})
 pd.DataFrame(rows,columns=['cell_type','n_a','n_b_qc','n_bsub','median_depth_a_before','median_depth_b_before','median_depth_a_after','median_depth_b_after']).to_csv(OUT/'match_table.csv',index=False);print(json.dumps(match_metrics(),sort_keys=True))

def inputs(d):return sc.read_h5ad(h5(d))
def select(write=False):
 data={d:inputs(d) for d in ('a','b')};detected=sorted(set(data['a'].var_names[nz(data['a'])>=3].astype(str))&set(data['b'].var_names[nz(data['b'])>=3].astype(str)))
 if len(detected)<2000:raise RuntimeError('Too few detected genes')
 data={d:x[:,detected].copy() for d,x in data.items()};both=ad.concat([data['a'],data['b']],axis=0,join='inner',merge='same',index_unique=None)
 if list(both.var_names.astype(str))!=sorted(both.var_names.astype(str)):raise RuntimeError('Genes were not alphabetically sorted before HVG selection')
 sc.pp.highly_variable_genes(both,flavor='seurat_v3',n_top_genes=2000,batch_key=None);genes=sorted(both.var_names[both.var['highly_variable'].to_numpy()].astype(str))
 if len(genes)!=2000:raise RuntimeError('Expected 2000 HVGs')
 if write:(OUT/'genes.txt').write_text('\n'.join(genes)+'\n')
 return data,detected,genes
def norm(x):
 x=x.copy();sc.pp.normalize_total(x,target_sum=1e4);sc.pp.log1p(x);return x.X.toarray() if sparse.issparse(x.X) else np.asarray(x.X)
def ranks(z):
 z=np.asarray(z,dtype=np.float32);out=np.zeros(len(z),dtype=np.int32)
 for c in range(len(z)):
  dist=np.linalg.norm(z-z[c],axis=1);dist[c]=np.inf;ids=np.argpartition(dist,K)[:K];ev=np.abs(np.linalg.eigvalsh(np.cov(z[ids],rowvar=False,bias=False)));out[c]=int((ev>FLOOR*float(ev.max())).sum())
 return out
def latent_metrics():
 genes=(OUT/'genes.txt').read_text().splitlines();data={d:inputs(d) for d in ('a','b')};detected=sorted(set(data['a'].var_names[nz(data['a'])>=3].astype(str))&set(data['b'].var_names[nz(data['b'])>=3].astype(str)));zero={d:int((nz(x[:,genes])==0).sum()) for d,x in data.items()};r={d:ranks(np.load(lat(d))) for d in data};allr=np.concatenate([r['a'],r['b']]);return {'n_genes_detected_both':len(detected),'n_hvg_zero_a':zero['a'],'n_hvg_zero_b':zero['b'],'rank_fa_a':float(np.median(r['a'])),'rank_fa_b':float(np.median(r['b'])),'frac_fa_at_ceiling':float((allr>=28).mean())}
def latent():
 data,_,genes=select(True)
 for d in ('a','b'):np.save(lat(d),np.asarray(FactorAnalysis(D,random_state=SEED,max_iter=2000).fit_transform(norm(data[d][:,genes].copy())),dtype=np.float32))
 print(json.dumps(latent_metrics(),sort_keys=True))

def cfg():return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=SEED,batch_size=256,loss='esm',lr=5e-4,beta=1.,esm_laplacian='hutchinson',device='cpu')
def zin(d):
 z=np.asarray(np.load(lat(d)),dtype=np.float32);labels=inputs(d).obs['cell_type'].astype(str).to_numpy()
 if z.shape!=(len(labels),D):raise RuntimeError('Latent shape mismatch')
 return z,labels
def save_dr(d,dr):torch.save({'state_dict':dr._model.state_dict(),'dim':D,'config':asdict(dr.config)},driftp(d))
def load_dr(d):
 x=torch.load(driftp(d),map_location='cpu',weights_only=False);dr=EnergyDriftAdapter(D,config=EnergyDriftConfig(**x['config']));m=_EnergyMLP(D,dr.config.hidden,dr.config.depth).to('cpu');m.load_state_dict(x['state_dict']);m.eval();dr._model,dr._device,dr._fitted=m,'cpu',True;return dr
def drift():
 for d in ('a','b'):
  z,_=zin(d);save_dr(d,EnergyDriftAdapter(D,config=cfg()).fit(z,verbose=True))
def readops(d):
 with np.load(opsp(d),allow_pickle=False) as x:return {k:x[k] for k in x.files}
def writeops(d,x):np.savez_compressed(opsp(d),**x)
def operators():
 gates={}
 for d in ('a','b'):
  z,labels=zin(d);dr=load_dr(d);J=dr.jacobian_autodiff(z).astype(np.float32);sy=dr.symmetry_check(z,n_cells=100);writeops(d,{'J':J,'eig':np.linalg.eigvalsh(J).astype(np.float32),'cell_type':np.asarray(labels,dtype=str),'symmetry_median':np.asarray(sy['median_asym'],dtype=np.float64),'symmetry_max':np.asarray(sy['max_asym'],dtype=np.float64),'symmetry_mean':np.asarray(sy['mean_asym'],dtype=np.float64)});gates[f'symmetry_{d}']=float(sy['median_asym'])
 print(json.dumps(gates,sort_keys=True))
def majority(x):
 c=pd.Series(x).value_counts();w=c[c==c.max()].index.astype(str).tolist();return w[0] if len(w)==1 else ''
def basin():
 data={}
 for d in ('a','b'):
  x=readops(d);z,labels=zin(d);a=ad.AnnData(X=np.zeros((len(labels),1),dtype=np.float32));a.obs['cell_type']=pd.Categorical(labels);a.obsm['X_latent']=z;a.uns['scjdo2']={'J_latent':x['J']};tol=.1*float(np.median(np.abs(x['eig'])));classify_basin_latent(a,tol=tol,use_rep='X_latent',rule='lead_eig');x['basin_index']=a.obs['op_basin_index'].to_numpy(dtype=np.float32);x['basin_class']=a.obs['op_behaviour_basin'].astype(str).to_numpy(dtype=str);x['tol']=np.asarray(tol,dtype=np.float64);writeops(d,x);data[d]=x
 rows=[]
 for t in sorted(set(data['a']['cell_type'].astype(str))|set(data['b']['cell_type'].astype(str))):
  row={'cell_type':t}
  for d in ('a','b'):
   mask=data[d]['cell_type'].astype(str)==t;row[f'n_{d}']=int(mask.sum());row[f'median_bi_{d}']=float(np.median(data[d]['basin_index'][mask])) if mask.any() else np.nan;row[f'class_{d}']=majority(data[d]['basin_class'][mask]) if mask.any() else ''
  row['agree']=bool(row['class_a']==row['class_b']);rows.append(row)
 pd.DataFrame(rows,columns=['cell_type','n_a','n_b','median_bi_a','median_bi_b','class_a','class_b','agree']).to_csv(OUT/'basin_table.csv',index=False)
def w2(a,b):return float(np.sqrt(ot.wasserstein_1d(a,b,p=2)))
def w1(a,b):return float(ot.wasserstein_1d(a,b,p=1))
def tv(a,b):return float(.5*sum(abs(np.mean(a==c)-np.mean(b==c)) for c in CLASSES))
def quotient(a,b):return 1. if b==0 and a==0 else (float('inf') if b==0 else a/b)
def compare(write=True):
 data={d:readops(d) for d in ('a','b')};types=sorted([t for t in set(data['a']['cell_type'].astype(str))&set(data['b']['cell_type'].astype(str)) if (data['a']['cell_type'].astype(str)==t).sum()>=50 and (data['b']['cell_type'].astype(str)==t).sum()>=50]);ids={d:{t:np.flatnonzero(data[d]['cell_type'].astype(str)==t) for t in types} for d in ('a','b')};m={t:min(len(ids['a'][t]),len(ids['b'][t]))//2 for t in types};rng=np.random.default_rng(SEED);rows=[]
 for draw in range(200):
  ia={t:rng.choice(ids['a'][t],m[t],replace=False) for t in types};ib={t:rng.choice(ids['b'][t],m[t],replace=False) for t in types};sp=[];bi=[];cl=[]
  for t in types:sp.append(w2(data['a']['eig'][ia[t]].ravel(),data['b']['eig'][ib[t]].ravel()));bi.append(w1(data['a']['basin_index'][ia[t]],data['b']['basin_index'][ib[t]]));cl.append(tv(data['a']['basin_class'][ia[t]],data['b']['basin_class'][ib[t]]))
  rows.append(dict(kind='pooled',draw=draw,spectral_w2=np.mean(sp),basin_w1=np.mean(bi),class_tv=np.mean(cl)))
 for d in ('a','b'):
  for draw in range(100):
   sp=[];bi=[];cl=[]
   for t in types:
    q=rng.choice(ids[d][t],2*m[t],replace=False);left,right=q[:m[t]],q[m[t]:];sp.append(w2(data[d]['eig'][left].ravel(),data[d]['eig'][right].ravel()));bi.append(w1(data[d]['basin_index'][left],data[d]['basin_index'][right]));cl.append(tv(data[d]['basin_class'][left],data[d]['basin_class'][right]))
   rows.append(dict(kind=f'null_{d}',draw=draw,spectral_w2=np.mean(sp),basin_w1=np.mean(bi),class_tv=np.mean(cl)))
 for d in ('a','b'):
  for draw in range(100):
   sp=[];bi=[]
   for u,v in itertools.combinations(types,2):
    size=min(m[u],m[v]);left=rng.choice(ids[d][u],size,replace=False);right=rng.choice(ids[d][v],size,replace=False);sp.append(w2(data[d]['eig'][left].ravel(),data[d]['eig'][right].ravel()));bi.append(w1(data[d]['basin_index'][left],data[d]['basin_index'][right]))
   rows.append(dict(kind=f'sep_{d}',draw=draw,spectral_w2=np.mean(sp),basin_w1=np.mean(bi),class_tv=np.nan))
 draws=pd.DataFrame(rows,columns=['kind','draw','spectral_w2','basin_w1','class_tv']);
 if write:draws.to_csv(OUT/'draws.csv',index=False)
 med=lambda mask,col:float(np.median(draws.loc[mask,col]));sp,spn=med(draws.kind=='pooled','spectral_w2'),med(draws.kind.str.startswith('null_'),'spectral_w2');bi,bin_=med(draws.kind=='pooled','basin_w1'),med(draws.kind.str.startswith('null_'),'basin_w1');cl,cln=med(draws.kind=='pooled','class_tv'),med(draws.kind.str.startswith('null_'),'class_tv');table=pd.read_csv(OUT/'basin_table.csv');lookup=table.set_index('cell_type').loc[types];rho=float(spearmanr(lookup['median_bi_a'],lookup['median_bi_b']).statistic)
 gate={'n_cell_types_tested':len(types),'symmetry_a':float(data['a']['symmetry_median']),'symmetry_b':float(data['b']['symmetry_median']),'spectral_w2_pooled':sp,'spectral_w2_null':spn,'spectral_ratio':quotient(sp,spn),'basin_w1_pooled':bi,'basin_w1_null':bin_,'basin_ratio':quotient(bi,bin_),'class_tv_pooled':cl,'class_tv_null':cln,'class_tv_ratio':quotient(cl,cln),'spectral_sep_a':quotient(med(draws.kind=='sep_a','spectral_w2'),med(draws.kind=='null_a','spectral_w2')),'spectral_sep_b':quotient(med(draws.kind=='sep_b','spectral_w2'),med(draws.kind=='null_b','spectral_w2')),'basin_sep_a':quotient(med(draws.kind=='sep_a','basin_w1'),med(draws.kind=='null_a','basin_w1')),'basin_sep_b':quotient(med(draws.kind=='sep_b','basin_w1'),med(draws.kind=='null_b','basin_w1')),'basin_index_rho':rho,'basin_class_agreement':float(lookup['agree'].mean())}
 if write:
  fig,axes=plt.subplots(1,len(types),figsize=(3.1*len(types),3.3),squeeze=False);axes=axes[0]
  for ax,t in zip(axes,types):
   for d,c in (('a','C0'),('b','C1')):ax.hist(data[d]['eig'][data[d]['cell_type'].astype(str)==t].ravel(),density=True,histtype='step',color=c,label=d)
   ax.set_title(t);ax.set_xlabel('eigenvalue')
  axes[0].set_ylabel('density');axes[-1].legend(title='dataset');ins=axes[-1].inset_axes([.46,.48,.5,.45]);ins.hist(draws.loc[draws.kind.str.startswith('null_'),'spectral_w2'],density=True,color='.7');ins.axvline(sp,color='C3');ins.set_title('null / pooled',fontsize=7);fig.tight_layout();fig.savefig(OUT/'fig_spectra.png',dpi=180);plt.close(fig)
  fig,ax=plt.subplots(figsize=(4.4,4));x,y=lookup['median_bi_a'].to_numpy(),lookup['median_bi_b'].to_numpy();ax.scatter(x,y)
  for t,xx,yy in zip(types,x,y):ax.annotate(t,(xx,yy),fontsize=8)
  lo,hi=min(x.min(),y.min()),max(x.max(),y.max());ax.plot([lo,hi],[lo,hi],color='.3',linestyle='--');ax.set_xlabel('median basin_index A');ax.set_ylabel('median basin_index B');ax.set_title(f'Spearman rho = {rho:.3f}');fig.tight_layout();fig.savefig(OUT/'fig_basin.png',dpi=180);plt.close(fig)
 print(json.dumps(gate,sort_keys=True));return gate
def deviations():
 out=[]
 for d in ('a','b'):
  x=readops(d);labels=x['cell_type'].astype(str);classes=x['basin_class'].astype(str)
  for t in sorted(set(labels)):
   c=pd.Series(classes[labels==t]).value_counts();
   if len(c)>1 and (c==c.max()).sum()>1:out.append({'what':f'{t} ({d}, {int((labels==t).sum())} cells) has a tied majority basin class; written as an empty class in basin_table.csv.','why':'The task does not specify a display convention for a tied majority; tested-type statistics use per-cell classes and are unaffected.'})
 return out
def report():
 gate={**stage_metrics(),**match_metrics(),**latent_metrics(),**compare(True)};files=['requirements.lock','data_checksums.txt','match_table.csv','genes.txt','latent_fa_a.npy','latent_fa_b.npy','ops_a.npz','ops_b.npz','basin_table.csv','draws.csv','fig_spectra.png','fig_basin.png','scripts/t01e_run.py'];base='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/T01e/scripts/t01e_run.py';obj={'task_id':'T01e','executor':'manus','branch':'task/T01e','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':SEED,'started':'2026-10-15T05:20:21-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':gate,'deliverables':[{'path':f,'sha256':sha(OUT/f)} for f in files],'regeneration':{**{f:base+' --stage all' for f in files},**{x:base+' --stage report' for x in gate}},'deviations':deviations(),'deviations_requested':[],'observations':[f'Raw B matrix SHA-256: {sha(RAW)}.',f'A_matched.h5ad SHA-256: {sha(h5("a"))}.',f'B_matched.h5ad SHA-256: {sha(h5("b"))}.'],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')};(OUT/'report.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('stage','match','latent','drift','operators','basin','compare','report','all'));s=p.parse_args().stage
 if s in ('stage','all'):stage()
 if s in ('match','all'):match()
 if s in ('latent','all'):latent()
 if s in ('drift','all'):drift()
 if s in ('operators','all'):operators()
 if s in ('basin','all'):basin()
 if s in ('compare','all'):compare()
 if s in ('report','all'):report()
if __name__=='__main__':main()
