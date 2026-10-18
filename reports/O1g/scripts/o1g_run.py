#!/usr/bin/env python
import argparse, hashlib, json, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
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

ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'reports'/'O1g'; CACHE=OUT/'cache'; DATA=ROOT/'data'/'pbmc_10k_v3'
SOURCE=DATA/'pbmc_10k_v3.h5ad'; HALVES=DATA/'k0_halves.csv'; K0=ROOT/'reports'/'K0'; GENES=K0/'genes.txt'; LAT1=K0/'latent_pca30_h1.npy'; LAT2=K0/'latent_pca30_h2.npy'
O1B=ROOT.parent/'scopalign-manus-O1b'/'reports'/'O1b'; O1B_DRIFTS=O1B/'cache'/'drifts'
D=30; FIRST5=(20261014,20261016,20261017,20261018,20261019); SEEDS20=(20261014,)+tuple(range(20261016,20261035)); NULL_SEEDS=(20261016,20261017,20261018,20261019); NS=(2000,2759); REPS=(0,1,2)
KNOWN=['tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]','tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]','tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]','tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility']
WORKER_PROBLEM=None

def chunks(f):
 while True:
  b=f.read(1<<20)
  if not b:return
  yield b
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in chunks(f):h.update(b)
 return h.hexdigest()
def fit_workers():return max(1,(os.cpu_count() or 1)-1)
def start_workers():return min(3,fit_workers())
def cfg(seed):return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=seed,batch_size=256,loss='esm',lr=5e-4,beta=1.0,esm_laplacian='hutchinson',device='cpu')
def local_fit(case,side,seed):return CACHE/'drifts'/f'{case}_{side}_{seed}.pt'
def old_fit(case,side,seed):return O1B_DRIFTS/f'{case}_{side}_{seed}.pt'
def load_saved(path):return torch.load(path,map_location='cpu',weights_only=False)
def row(case,side,seed,saved):return {'fit_id':f'{case}_{side}_{seed}','case':case,'side':side,'seed':int(seed),'median_asym':float(saved['symmetry']['median_asym']),'max_asym':float(saved['symmetry']['max_asym'])}
def make_meta():
 OUT.mkdir(parents=True,exist_ok=True); CACHE.mkdir(parents=True,exist_ok=True); (OUT/'.gitignore').write_text('cache/\n__pycache__/\nrun.log\nrun.pid\n')
 p=OUT/'run_meta.json'
 if not p.exists():p.write_text(json.dumps({'started':datetime.now(timezone.utc).isoformat()},sort_keys=True)+'\n')
 return json.loads(p.read_text())
def save_sym(rows):
 p=OUT/'fit_symmetry.csv'; old=pd.read_csv(p).to_dict('records') if p.exists() else []; d={x['fit_id']:x for x in old}; d.update({x['fit_id']:x for x in rows}); pd.DataFrame([d[k] for k in sorted(d)]).to_csv(p,index=False)
def fit_worker(item):
 torch.set_num_threads(1); case,side,seed,z,path=item; p=Path(path)
 if p.exists():return row(case,side,seed,load_saved(p))
 p.parent.mkdir(parents=True,exist_ok=True); drift=EnergyDriftAdapter(D,config=cfg(seed)).fit(np.asarray(z,dtype=np.float32),verbose=True); sym=drift.symmetry_check(np.asarray(z,dtype=np.float32),n_cells=100)
 saved={'state_dict':drift._model.state_dict(),'config':asdict(drift.config),'symmetry':{k:float(v) for k,v in sym.items()}}; torch.save(saved,p); return row(case,side,seed,saved)
def halves():
 a=sc.read_h5ad(SOURCE); a=a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy(); t=pd.read_csv(HALVES).set_index('barcode')['half']; names=a.obs_names.astype(str); m1=names.isin(t.index)&(names.map(t)=='H1');m2=names.isin(t.index)&(names.map(t)=='H2');return a[m1].copy(),a[m2].copy()
def dense(a):
 a=a.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
def pca_halves():
 h1,h2=halves();genes=[x.strip() for x in GENES.read_text().splitlines() if x.strip()]
 if len(genes)!=2000:raise RuntimeError('K0 genes list must contain exactly 2000 genes')
 p1=PCA(n_components=D,svd_solver='full',random_state=20261014);p2=PCA(n_components=D,svd_solver='full',random_state=20261014);z1=p1.fit_transform(dense(h1[:,genes]));z2=p2.fit_transform(dense(h2[:,genes]));u,_,vt=np.linalg.svd(p2.components_@p1.components_.T,full_matrices=False);q=(u@vt).astype(np.float64);chk=max(float(np.max(np.abs(z1-np.load(LAT1).astype(np.float64)))),float(np.max(np.abs(z2-np.load(LAT2).astype(np.float64)))));return z1.astype(np.float32),z2.astype(np.float32),q,chk,h1,h2
def case_inputs():
 z=np.load(LAT1).astype(np.float32)
 if len(z)!=5518:raise RuntimeError(f'Expected 5518 H1 cells, found {len(z)}')
 ans=[]
 for n in NS:
  for r in REPS:
   rng=np.random.default_rng(20261030+r);perm=rng.permutation(5518);q=np.asarray(ortho_group.rvs(D,random_state=20261031+r),dtype=np.float64);x=z[np.sort(perm[:n])];y=(z[np.sort(perm[n:2*n])].astype(np.float64)@q.T).astype(np.float32);ans.append((f'R_n{n}_r{r}',x,y))
 x,y,_,_,_,_=pca_halves();ans.append(('F',x,y));return ans
def fit_stage():
 make_meta(); rows=[];reused=[];todo=[]
 for case,x,y in case_inputs():
  for side,z in (('X',x),('Y',y)):
   for seed in SEEDS20:
    lp=local_fit(case,side,seed)
    if lp.exists():rows.append(row(case,side,seed,load_saved(lp)))
    elif seed in FIRST5:
     op=old_fit(case,side,seed)
     if not op.exists():raise RuntimeError(f'Required O1b fit is missing: {op}')
     saved=load_saved(op);rows.append(row(case,side,seed,saved));reused.append({'case':case,'side':side,'seed':int(seed),'source':str(op),'sha256':sha(op)})
    else:todo.append((case,side,seed,z,str(lp)))
 if todo:
  with ProcessPoolExecutor(max_workers=fit_workers()) as ex:rows.extend(list(ex.map(fit_worker,todo)))
 save_sym(rows);pd.DataFrame(reused,columns=['case','side','seed','source','sha256']).sort_values(['case','side','seed']).to_csv(OUT/'reused_o1b_fits.csv',index=False);print(json.dumps({'fit_rows':len(rows),'new_fits':len(todo),'reused_o1b_fits':len(reused),'fit_workers':fit_workers()},sort_keys=True),flush=True)
def build_model(saved):
 m=_EnergyMLP(D,int(saved['config']['hidden']),int(saved['config']['depth']));m.load_state_dict(saved['state_dict']);m=m.to(dtype=torch.float64,device='cpu').eval()
 for p in m.parameters():p.requires_grad_(False)
 return m
def get_models(case):
 mx=[];my=[];rows=[]
 for seed in SEEDS20:
  for side,bucket in (('X',mx),('Y',my)):
   lp=local_fit(case,side,seed);op=old_fit(case,side,seed)
   if lp.exists():s=load_saved(lp)
   elif seed in FIRST5 and op.exists():s=load_saved(op)
   else:raise RuntimeError(f'Missing fit after fit stage: {case}_{side}_{seed}')
   bucket.append(build_model(s));rows.append(row(case,side,seed,s))
 save_sym(rows);return mx,my
def fields(model,x,need_j,track):
 z=x.clone().requires_grad_(True);e=model(z);g=torch.autograd.grad(e.sum(),z,create_graph=(need_j or track),retain_graph=(need_j or track))[0];u=-g
 if not need_j:return (u,None) if track else (u.detach(),None)
 rs=[]
 for j in range(D):rs.append(-torch.autograd.grad(g[:,j].sum(),z,create_graph=track,retain_graph=(track or j<D-1))[0])
 j=torch.stack(rs,dim=1);return (u,j) if track else (u.detach(),j.detach())
def ensemble(models,x,need_j,track):
 us=[];js=[]
 for m in models:
  u,j=fields(m,x,need_j,track);us.append(u)
  if need_j:js.append(j)
 return torch.stack(us).mean(0),torch.stack(js).mean(0) if need_j else None
def individual_ops(model,z):
 out=[]
 for start in range(0,len(z),128):
  _,j=fields(model,torch.as_tensor(np.asarray(z[start:start+128],dtype=np.float64),dtype=torch.float64),True,False);out.append(j.cpu().numpy())
 return np.concatenate(out)
def starts():
 yield 0,np.eye(D,dtype=np.float64)
 for i in range(1,20):yield i,np.asarray(ortho_group.rvs(D,random_state=20263000+i),dtype=np.float64)
def start_overlap():
 best=0;ss=[q for _,q in starts()]
 for r in REPS:
  q=np.asarray(ortho_group.rvs(D,random_state=20261031+r),dtype=np.float64)
  for q0 in ss:best=max(best,int(np.linalg.norm(q0-q,ord='fro')<1e-8))
 return int(best)
class Problem:
 def __init__(self,mx,my,x,y,r):
  self.mx,self.my=mx,my;self.x=torch.as_tensor(np.asarray(x,dtype=np.float64),dtype=torch.float64);self.y=torch.as_tensor(np.asarray(y,dtype=np.float64),dtype=torch.float64);rng=np.random.default_rng(20261040+r);self.xa=self.x[rng.choice(len(x),size=min(256,len(x)),replace=False)];self.yb=self.y[rng.choice(len(y),size=min(256,len(y)),replace=False)];self.ux,_=ensemble(mx,self.x,False,False);self.uy,_=ensemble(my,self.y,False,False);_,self.jx=ensemble(mx,self.xa,True,False);_,self.jy=ensemble(my,self.yb,True,False);self.du=torch.sum(self.ux*self.ux);self.dv=torch.sum(self.uy*self.uy);self.dx=torch.sum(self.jx*self.jx);self.dy=torch.sum(self.jy*self.jy)
 def loss(self,q,track):
  uy,_=ensemble(self.my,self.x@q.T,False,track);ux,_=ensemble(self.mx,self.y@q,False,track);dl=(torch.sum((uy-self.ux@q.T)**2)/self.du+torch.sum((ux-self.uy@q)**2)/self.dv)/2;_,jy=ensemble(self.my,self.xa@q.T,True,track);_,jx=ensemble(self.mx,self.yb@q,True,track);qj=torch.einsum('ij,njk,lk->nil',q,self.jx,q);qtj=torch.einsum('ji,njk,kl->nil',q,self.jy,q);hl=(torch.sum((jy-qj)**2)/self.dx+torch.sum((jx-qtj)**2)/self.dy)/2;return dl+hl,dl,hl
class Closure:
 def __init__(self,p,q0,s):self.p,self.q0,self.s=p,q0,s
 def __call__(self):
  self.s.grad=None;q=self.q0@torch.matrix_exp(self.s-self.s.T);v,_,_=self.p.loss(q,True);v.backward();return v
def init_start_worker(problem):
 global WORKER_PROBLEM;torch.set_num_threads(1);WORKER_PROBLEM=problem
def optimise_start(item):
 torch.set_num_threads(1);start,q0n=item;q0=torch.as_tensor(q0n,dtype=torch.float64);s=torch.zeros((D,D),dtype=torch.float64,requires_grad=True);opt=torch.optim.LBFGS([s],lr=1.0,max_iter=500,tolerance_grad=1e-10,tolerance_change=1e-12,line_search_fn='strong_wolfe');opt.step(Closure(WORKER_PROBLEM,q0,s));q=(q0@torch.matrix_exp(s-s.T)).detach();v,dl,hl=WORKER_PROBLEM.loss(q,False);return {'start':int(start),'L_hat':float(v.detach()),'D':float(dl.detach()),'H':float(hl.detach()),'Q':q.cpu().numpy()}
def sortkey(x):return(x['L_hat'],x['start'])
def estimate(problem):
 with ProcessPoolExecutor(max_workers=start_workers(),initializer=init_start_worker,initargs=(problem,)) as ex:ans=list(ex.map(optimise_start,starts()))
 return min(ans,key=sortkey),ans
def test_r_stage():
 make_meta();z=np.load(LAT1).astype(np.float32);p=OUT/'testR.csv';rows=pd.read_csv(p).to_dict('records') if p.exists() else [];done={(int(x['n']),int(x['replicate'])) for x in rows}
 for n in NS:
  for r in REPS:
   if (n,r) in done:continue
   rng=np.random.default_rng(20261030+r);perm=rng.permutation(5518);qstar=np.asarray(ortho_group.rvs(D,random_state=20261031+r),dtype=np.float64);x=z[np.sort(perm[:n])];y=(z[np.sort(perm[n:2*n])].astype(np.float64)@qstar.T).astype(np.float32);mx,my=get_models(f'R_n{n}_r{r}');problem=Problem(mx,my,x,y,r);best,_=estimate(problem);lstar,_,_=problem.loss(torch.as_tensor(qstar,dtype=torch.float64),False);record={'n':int(n),'replicate':int(r),'err':float(np.linalg.norm(best['Q']-qstar,ord='fro')/np.sqrt(60.0)),'L_hat':float(best['L_hat']),'L_star':float(lstar.detach()),'best_start':int(best['start']),'D':float(best['D']),'H':float(best['H'])};rows.append(record);pd.DataFrame(rows).sort_values(['n','replicate']).to_csv(p,index=False);print(json.dumps(record,sort_keys=True),flush=True)
 df=pd.DataFrame(rows).sort_values(['n','replicate']);df.to_csv(p,index=False);meta={'start_overlap':start_overlap(),'frac_hat_below_star':float(np.mean(df.L_hat.to_numpy()<=df.L_star.to_numpy()))};(OUT/'testR_meta.json').write_text(json.dumps(meta,sort_keys=True)+'\n');print(json.dumps(meta,sort_keys=True),flush=True)
def residual(a,b,w):
 w=np.asarray(w,dtype=np.float64);w=w/w.sum();d=a-b;return float(np.sqrt(np.sum(w*np.sum(d*d,axis=(1,2)))/np.sum(w*np.sum(b*b,axis=(1,2)))))
def frame_ratio(z1,z2,l1,l2,m1,m2,q):
 groups=sorted(t for t in set(l1)&set(l2) if np.sum(l1==t)>=100 and np.sum(l2==t)>=100)
 if not groups:raise RuntimeError('No T02c centroid clusters for frame residual')
 w=np.asarray([min(np.sum(l1==t),np.sum(l2==t)) for t in groups],dtype=np.float64);c1=np.stack([centroid_operator(z1[l1==t],m1[20261014][l1==t],n_cells_centroid=30).astype(np.float64) for t in groups]);c2=np.stack([centroid_operator(z2[l2==t],m2[20261014][l2==t],n_cells_centroid=30).astype(np.float64) for t in groups]);num=residual(np.einsum('ij,kjl,ml->kim',q,c1,q),c2,w);rng=np.random.default_rng(20261022);vals=[]
 for z,l,m in ((z1,l1,m1),(z2,l2,m2)):
  for seed in NULL_SEEDS:
   a=[];b=[]
   for t in groups:
    idx=np.flatnonzero(l==t);perm=rng.permutation(idx);h=len(perm)//2;a.append(centroid_operator(z[perm[:h]],m[20261014][perm[:h]],n_cells_centroid=30).astype(np.float64));b.append(centroid_operator(z[perm[h:2*h]],m[seed][perm[h:2*h]],n_cells_centroid=30).astype(np.float64))
   vals.append(residual(np.stack(a),np.stack(b),w))
 den=float(np.median(vals));return float(num/den) if den else float('inf')
def umap(rep):
 a=ad.AnnData(X=np.zeros((len(rep),1),dtype=np.float32));a.obsm['view']=np.asarray(rep,dtype=np.float64);sc.pp.neighbors(a,n_neighbors=15,use_rep='view',random_state=20261014);sc.tl.umap(a,min_dist=.5,random_state=20261014);return np.asarray(a.obsm['X_umap'],dtype=np.float64)
def mixing(rep,ds):
 nb=NearestNeighbors(n_neighbors=31).fit(rep).kneighbors(rep,return_distance=False)[:,1:];return float(np.mean(np.mean(ds[nb]!=ds[:,None],axis=1))/.5)
def views(z1,z2,qhat,qload,h1,h2,ops):
 l1=h1.obs['cell_type'].astype(str).to_numpy();l2=h2.obs['cell_type'].astype(str).to_numpy();ds=np.asarray(['H1']*len(z1)+['H2']*len(z2),dtype=str);types=np.concatenate([l1,l2]).astype(str);tested=sorted(t for t in set(l1)&set(l2) if np.sum(l1==t)>=50 and np.sum(l2==t)>=50);mask=np.isin(types,tested);m1={s:ops[('H1',s)] for s in FIRST5};m2={s:ops[('H2',s)] for s in FIRST5};reps={'separate_frames':(np.vstack([z1,z2]),np.eye(D)),'E5b_estimate':(np.vstack([z1@qhat.T,z2]),qhat),'known_answer':(np.vstack([z1@qload.T,z2]),qload)};rows=[];coords={}
 for name,(rep,q) in reps.items():
  mix=mixing(rep,ds);sil=float(silhouette_score(rep[mask],types[mask],metric='euclidean'));fr=frame_ratio(z1.astype(np.float64),z2.astype(np.float64),l1,l2,m1,m2,q);rows.append({'view':name,'mix':mix,'sil':sil,'frame_ratio':fr});coords[name]=umap(rep)
 pd.DataFrame(rows).to_csv(OUT/'views.csv',index=False);fig,axes=plt.subplots(3,2,figsize=(10,13));colors=plt.get_cmap('tab10')
 for i,rowx in enumerate(rows):
  xy=coords[rowx['view']];ax=axes[i,0]
  for v,c in (('H1','C0'),('H2','C1')):
   mm=ds==v;ax.scatter(xy[mm,0],xy[mm,1],s=4,alpha=.7,color=c,label=v,linewidths=0)
  ax.set_title(f"{rowx['view']}: half; mix={rowx['mix']:.3f}, sil={rowx['sil']:.3f}, frame={rowx['frame_ratio']:.3f}");ax.legend(markerscale=2,fontsize=8);ax.set_xticks([]);ax.set_yticks([]);ax=axes[i,1]
  for j,t in enumerate(sorted(set(types))):
   mm=types==t;ax.scatter(xy[mm,0],xy[mm,1],s=4,alpha=.7,color=colors(j%10),label=t,linewidths=0)
  ax.set_title(f"{rowx['view']}: cell type; mix={rowx['mix']:.3f}, sil={rowx['sil']:.3f}, frame={rowx['frame_ratio']:.3f}");ax.legend(markerscale=2,fontsize=7,ncol=2);ax.set_xticks([]);ax.set_yticks([])
 fig.tight_layout();fig.savefig(OUT/'fig_umap.png',dpi=180);plt.close(fig);return {x['view']:x for x in rows}
def test_f_stage():
 make_meta();x,y,qload,chk,h1,h2=pca_halves();mx,my=get_models('F');best,_=estimate(Problem(mx,my,x,y,0));np.save(OUT/'Q_hat_F.npy',best['Q']);np.save(OUT/'Q_load_F.npy',qload);ops={}
 for side,models,z in (('H1',mx,x),('H2',my,y)):
  for seed,model in zip(FIRST5,models[:5]):
   p=CACHE/'operators'/f'F_{side}_{seed}.npy';p.parent.mkdir(parents=True,exist_ok=True)
   if not p.exists():np.save(p,individual_ops(model,z))
   ops[(side,seed)]=np.load(p)
 v=views(x,y,best['Q'],qload,h1,h2,ops);e=v['E5b_estimate'];n=v['separate_frames'];rec={'mix_F':float(e['mix']),'mix_F_noA':float(n['mix']),'sil_F':float(e['sil']),'sil_F_noA':float(n['sil']),'frame_ratio_F':float(e['frame_ratio']),'frame_ratio_F_noA':float(n['frame_ratio']),'err_F':float(np.linalg.norm(best['Q']-qload,ord='fro')/np.sqrt(60.0)),'best_start':int(best['start']),'L_hat':float(best['L_hat']),'latent_check_max':float(chk)};pd.DataFrame([rec]).to_csv(OUT/'testF.csv',index=False);print(json.dumps(rec,sort_keys=True),flush=True)
def report_stage():
 make_meta();subprocess.run([sys.executable,'-m','pip','freeze'],check=True,stdout=(OUT/'requirements.lock').open('w'));tr=pd.read_csv(OUT/'testR.csv');tf=pd.read_csv(OUT/'testF.csv').iloc[0];sym=pd.read_csv(OUT/'fit_symmetry.csv');meta=json.loads((OUT/'testR_meta.json').read_text());med={n:float(tr.loc[tr.n==n,'err'].median()) for n in NS};small=[n for n in NS if med[n]<=.10];tfpass=bool(float(tf.mix_F)>=.80 and float(tf.frame_ratio_F)<=1.5 and float(tf.sil_F)>=.9*float(tf.sil_F_noA));g={'err_R_E5b_2000':med[2000],'err_R_E5b_2759':med[2759],'frac_hat_below_star':float(meta['frac_hat_below_star']),'start_overlap':int(meta['start_overlap']),'mix_F':float(tf.mix_F),'mix_F_noA':float(tf.mix_F_noA),'sil_F':float(tf.sil_F),'sil_F_noA':float(tf.sil_F_noA),'frame_ratio_F':float(tf.frame_ratio_F),'err_F':float(tf.err_F),'symmetry_max':float(sym.median_asym.max()),'latent_check_max':float(tf.latent_check_max),'e5b_outcome':'works' if small and tfpass else 'fails','min_n_e5b':int(min(small)) if small else 0};fig,ax=plt.subplots(figsize=(7,4.5));ax.plot(list(NS),[med[n] for n in NS],marker='o',label='E5b');old=ROOT/'reports'/'O1b'/'testR.csv'
 if old.exists():
  o=pd.read_csv(old);avail=[n for n in NS if n in set(o.n)];ax.plot(avail,[float(o.loc[o.n==n,'err'].median()) for n in avail],marker='o',label='O1b E5')
 ax.axhline(.10,color='black',linestyle='--',linewidth=1,label='0.10');ax.set_xscale('log');ax.set_xticks(list(NS));ax.set_xticklabels([str(n) for n in NS]);ax.set_xlabel('cells per copy');ax.set_ylabel('median rotation error');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'fig_recovery.png',dpi=180);plt.close(fig)
 files=['requirements.lock','testR.csv','testF.csv','Q_hat_F.npy','views.csv','fig_umap.png','fig_recovery.png','scripts/o1g_run.py'];cmd='cd "$HOME/Documents/scopalign-manus-O1g" && .venv/bin/python reports/O1g/scripts/o1g_run.py';regen={p:cmd+' --stage all' for p in files};regen.update({'err_R_E5b_2000':cmd+' --stage testR','err_R_E5b_2759':cmd+' --stage testR','frac_hat_below_star':cmd+' --stage testR','start_overlap':cmd+' --stage testR','mix_F':cmd+' --stage testF','mix_F_noA':cmd+' --stage testF','sil_F':cmd+' --stage testF','sil_F_noA':cmd+' --stage testF','frame_ratio_F':cmd+' --stage testF','err_F':cmd+' --stage testF','symmetry_max':cmd+' --stage fits','latent_check_max':cmd+' --stage testF','e5b_outcome':cmd+' --stage report','min_n_e5b':cmd+' --stage report'});rm=json.loads((OUT/'run_meta.json').read_text());arts=[{'path':'fit_symmetry.csv','sha256':sha(OUT/'fit_symmetry.csv')},{'path':'reused_o1b_fits.csv','sha256':sha(OUT/'reused_o1b_fits.csv')},{'path':'testR_meta.json','sha256':sha(OUT/'testR_meta.json')},{'path':'Q_load_F.npy','sha256':sha(OUT/'Q_load_F.npy')}];rep={'task_id':'O1g','executor':'manus','branch':'task/O1g','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':rm['started'],'finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':g,'deliverables':[{'path':p,'sha256':sha(OUT/p)} for p in files],'artifacts':arts,'regeneration':regen,'deviations':[{'what':'limited concurrent E5 start workers to three','why':'The 20-model ensemble exhausted the local Mac memory after the completed fit stage; all 20 fixed starts, seeds, and float64 optimizer settings remain unchanged.'}],'deviations_requested':[],'observations':[f'Input SHA-256: {SOURCE.relative_to(ROOT)} = {sha(SOURCE)}.',f'Input SHA-256: {HALVES.relative_to(ROOT)} = {sha(HALVES)}.',f'Input SHA-256: {GENES.relative_to(ROOT)} = {sha(GENES)}.',f'Input SHA-256: {LAT1.relative_to(ROOT)} = {sha(LAT1)}.',f'Input SHA-256: {LAT2.relative_to(ROOT)} = {sha(LAT2)}.','The reused_o1b_fits.csv artifact records the source path and SHA-256 for every reused O1b fit.'],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')};(OUT/'report.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n');print(json.dumps(g,sort_keys=True),flush=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('fits','testR','testF','report','all'));stage=p.parse_args().stage
 if stage in ('fits','all'):fit_stage()
 if stage in ('testR','all'):test_r_stage()
 if stage in ('testF','all'):test_f_stage()
 if stage in ('report','all'):report_stage()
if __name__=='__main__':main()
