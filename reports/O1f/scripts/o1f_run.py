#!/usr/bin/env python
import argparse, hashlib, json, os, subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import ortho_group
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
import torch

ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'reports'/'O1f'; CASES=OUT/'cases'; O1=ROOT/'reports'/'O1'; K0=ROOT/'reports'/'K0'; DATA=ROOT/'data'/'pbmc_10k_v3'
D=30; NS=(500,1000,2000,2759); REPS=(0,1,2); CASE_INPUTS=tuple(sorted((O1/'cases').glob('case_*.npz'))); TESTF=O1/'testF.npz'; SOURCE=DATA/'pbmc_10k_v3.h5ad'; HALVES=DATA/'k0_halves.csv'; LAT1=K0/'latent_pca30_h1.npy'; LAT2=K0/'latent_pca30_h2.npy'; O1R=O1/'testR.csv'
KNOWN=['tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]','tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]','tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]','tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def gatesp():return OUT/'gates.json'
def gates():return json.loads(gatesp().read_text()) if gatesp().exists() else {}
def saveg(g):gatesp().write_text(json.dumps(g,sort_keys=True)+'\n')
def c_in(n,r):return O1/'cases'/f'case_{n}_{r}.npz'
def c_out(n,r):return CASES/f'case_{n}_{r}.npz'
def c_rand(n,r):return CASES/f'case_{n}_{r}_random.npz'
def f_result():return OUT/'testF_result.npz'
def sym(a):return .5*(a+np.swapaxes(a,-1,-2))
def groups(mx,my,tx,ty):
 g=[]
 for t in sorted(set(tx)&set(ty)):
  ix=np.flatnonzero(tx==t);iy=np.flatnonzero(ty==t)
  if len(ix)>=20 and len(iy)>=20:g.append((t,ix,iy))
 if not g:raise RuntimeError('no matched groups')
 gx=np.stack([mx[ix].mean(axis=0) for _,ix,_ in g]);gy=np.stack([my[iy].mean(axis=0) for _,_,iy in g]);w=np.asarray([min(len(ix),len(iy)) for _,ix,iy in g],dtype=np.float64)
 return [x[0] for x in g],gx,gy,w
def torch_obj(gx,gy,w):
 x=torch.as_tensor(gx,dtype=torch.float64);y=torch.as_tensor(gy,dtype=torch.float64);wt=torch.as_tensor(w,dtype=torch.float64);den=torch.sum(wt*torch.sum(y*y,dim=(1,2)))
 def obj(q):return torch.sum(wt*torch.sum((q@x@q.T-y)**2,dim=(1,2)))/den
 return obj
def optimise(q0,obj,max_iter):
 q0=torch.as_tensor(q0,dtype=torch.float64);s=torch.zeros((D,D),dtype=torch.float64,requires_grad=True);opt=torch.optim.LBFGS([s],lr=1,max_iter=max_iter,tolerance_grad=1e-10,tolerance_change=1e-12,line_search_fn='strong_wolfe')
 def qnow():return q0@torch.matrix_exp(s-s.T)
 def closure():opt.zero_grad();v=obj(qnow());v.backward();return v
 opt.step(closure)
 with torch.no_grad():q=qnow().numpy();loss=float(obj(torch.as_tensor(q,dtype=torch.float64)))
 return q,loss
def eigen_candidate(gx,gy,w,k):
 _,vx=np.linalg.eigh(gx[k]);_,vy=np.linalg.eigh(gy[k]);obj=torch_obj(gx,gy,w)
 for _ in range(10):
  changed=False
  for j in range(D):
   q=vy@vx.T;old=float(obj(torch.as_tensor(q,dtype=torch.float64)));vx[:,j]*=-1;new=float(obj(torch.as_tensor(vy@vx.T,dtype=torch.float64)))
   if new<old:changed=True
   else:vx[:,j]*=-1
  if not changed:break
 return vy@vx.T
def search(gx,gy,w,tag):
 rp=CASES/f'{tag}_random.npz';obj=torch_obj(gx,gy,w)
 if rp.exists():
  with np.load(rp,allow_pickle=False) as z:starts=z['Q'];losses=z['loss']
 else:
  rec=[]
  for i in range(2000):
   q,l=optimise(ortho_group.rvs(D,random_state=20262000+i),obj,50);rec.append((l,i,q))
  rec.sort(key=lambda x:(x[0],x[1]));top=rec[:50];starts=np.stack([x[2] for x in top]);losses=np.asarray([x[0] for x in top]);np.savez_compressed(rp,Q=starts,loss=losses)
 candidates=[(f'random_{i}',starts[i]) for i in range(50)]+[(f'eigen_{k}',eigen_candidate(gx,gy,w,k)) for k in range(len(w))]
 best=(None,np.inf,'')
 for source,q0 in candidates:
  q,l=optimise(q0,obj,500)
  if l<best[1]:best=(q,l,source)
 return best[0],best[1],best[2]
def one_case(n,r):
 op=c_out(n,r)
 if op.exists():return str(op)
 with np.load(c_in(n,r),allow_pickle=False) as z:
  mx=sym(z['M_X'].astype(np.float64).mean(axis=0));my=sym(z['M_Y'].astype(np.float64).mean(axis=0));tx=z['type_X'].astype(str);ty=z['type_Y'].astype(str);target=z['Q_star'].astype(np.float64)
 labels,gx,gy,w=groups(mx,my,tx,ty);q,lhat,source=search(gx,gy,w,f'case_{n}_{r}');obj=torch_obj(gx,gy,w);lstar=float(obj(torch.as_tensor(target,dtype=torch.float64)));er=float(np.linalg.norm(q-target,'fro')/np.sqrt(60.));np.savez_compressed(op,n=np.asarray(n),replicate=np.asarray(r),groups=np.asarray(labels,dtype=str),Q_hat=q,L_hat=np.asarray(lhat),L_star=np.asarray(lstar),err=np.asarray(er),best_source=np.asarray(source),reached_truth=np.asarray(lhat<=1.01*lstar),alt_minimum=np.asarray((lhat<.99*lstar) and (er>.10)))
 return str(op)
def run_case_pair(pair):return one_case(*pair)
def testR():
 CASES.mkdir(parents=True,exist_ok=True);todo=[(n,r) for n in NS for r in REPS if not c_out(n,r).exists()]
 if todo:
  with ProcessPoolExecutor(max_workers=4) as ex:list(ex.map(run_case_pair,todo))
 rows=[]
 for n in NS:
  for r in REPS:
   with np.load(c_out(n,r),allow_pickle=False) as z:rows.append({'n':n,'replicate':r,'err':float(z['err']),'L_hat':float(z['L_hat']),'L_star':float(z['L_star']),'best_source':str(z['best_source']),'reached_truth':bool(z['reached_truth']),'alt_minimum':bool(z['alt_minimum'])})
 df=pd.DataFrame(rows);df.to_csv(OUT/'testR.csv',index=False);g=gates()
 for n in NS:g[f'err_R_E9_{n}']=float(np.median(df[df.n==n].err))
 g['frac_reached_truth']=float(df.reached_truth.mean());g['n_alt_minimum']=int(df.alt_minimum.sum());saveg(g);print(json.dumps({k:g[k] for k in g if k.startswith('err_R') or k in ('frac_reached_truth','n_alt_minimum')},sort_keys=True))
def halves_labels():
 a=sc.read_h5ad(SOURCE);a=a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy();tab=pd.read_csv(HALVES).set_index('barcode')['half'];names=a.obs_names.astype(str);h1=names.isin(tab.index)&(names.map(tab)=='H1');h2=names.isin(tab.index)&(names.map(tab)=='H2');return a[h1].obs['cell_type'].astype(str).to_numpy(),a[h2].obs['cell_type'].astype(str).to_numpy()
def mix_sil(z,ds,labels):
 nn=NearestNeighbors(n_neighbors=31).fit(z).kneighbors(z,return_distance=False)[:,1:];mix=float(np.mean(np.mean(ds[nn]!=ds[:,None],axis=1))/.5);good=sorted(t for t in set(labels) if np.sum(labels==t)>=50);sil=float(silhouette_score(z[np.isin(labels,good)],labels[np.isin(labels,good)]));return mix,sil
def umap(z):
 a=sc.AnnData(X=np.zeros((len(z),1),dtype=np.float32));a.obsm['view']=z;sc.pp.neighbors(a,n_neighbors=15,use_rep='view',random_state=20261014);sc.tl.umap(a,min_dist=.5,random_state=20261014);return a.obsm['X_umap']
def views(qhat,qload):
 z1,z2=np.load(LAT1).astype(np.float64),np.load(LAT2).astype(np.float64);t1,t2=halves_labels();ds=np.asarray(['H1']*len(z1)+['H2']*len(z2));labels=np.concatenate([t1,t2]);items=[('no_rotation',np.vstack([z1,z2])),('E9',np.vstack([z1@qhat.T,z2])),('Q_load',np.vstack([z1@qload.T,z2]))];rows=[];fig,axes=plt.subplots(3,2,figsize=(10,12));cmap=plt.get_cmap('tab10')
 for i,(name,z) in enumerate(items):
  m,s=mix_sil(z,ds,labels);rows.append({'view':name,'mix':m,'sil':s});xy=umap(z)
  for d,c in (('H1','C0'),('H2','C1')):axes[i,0].scatter(xy[ds==d,0],xy[ds==d,1],s=4,alpha=.7,color=c,label=d,linewidths=0)
  axes[i,0].set_title(f'{name}: dataset; mix={m:.3f}, sil={s:.3f}');axes[i,0].legend(fontsize=8);axes[i,0].set_xticks([]);axes[i,0].set_yticks([])
  for j,t in enumerate(sorted(set(labels))):axes[i,1].scatter(xy[labels==t,0],xy[labels==t,1],s=4,alpha=.7,color=cmap(j%10),label=t,linewidths=0)
  axes[i,1].set_title(f'{name}: cell type; mix={m:.3f}, sil={s:.3f}');axes[i,1].legend(fontsize=7,ncol=2);axes[i,1].set_xticks([]);axes[i,1].set_yticks([])
 fig.tight_layout();fig.savefig(OUT/'fig_umap.png',dpi=180);plt.close(fig);pd.DataFrame(rows).to_csv(OUT/'views.csv',index=False)
def testF():
 if f_result().exists():
  with np.load(f_result(),allow_pickle=False) as z:q=z['Q_hat'];lhat=float(z['L_hat']);lstar=float(z['L_star']);er=float(z['err']);source=str(z['best_source']);qload=z['Q_load']
 else:
  with np.load(TESTF,allow_pickle=False) as z:mx=sym(z['Mbar_X'].astype(np.float64));my=sym(z['Mbar_Y'].astype(np.float64));qload=z['Q_load'].astype(np.float64)
  tx,ty=halves_labels();labels,gx,gy,w=groups(mx,my,tx,ty);q,lhat,source=search(gx,gy,w,'testF');obj=torch_obj(gx,gy,w);lstar=float(obj(torch.as_tensor(qload,dtype=torch.float64)));er=float(np.linalg.norm(q-qload,'fro')/np.sqrt(60.));np.savez_compressed(f_result(),groups=np.asarray(labels,dtype=str),Q_hat=q,Q_load=qload,L_hat=np.asarray(lhat),L_star=np.asarray(lstar),err=np.asarray(er),best_source=np.asarray(source),reached_truth=np.asarray(lhat<=1.01*lstar))
 pd.DataFrame([{'err':er,'L_hat':lhat,'L_star':lstar,'best_source':source,'reached_truth':lhat<=1.01*lstar}]).to_csv(OUT/'testF.csv',index=False);views(q,qload);g=gates();g['err_F_E9']=er;saveg(g);print(json.dumps({'err_F_E9':er},sort_keys=True))
def report():
 g=gates();works=[n for n in NS if g[f'err_R_E9_{n}']<=.10];g['e9_outcome']='works' if g['err_F_E9']<=.25 and works else 'fails';g['min_n_e9']=int(min(works)) if works else 0;saveg(g);d=pd.read_csv(OUT/'testR.csv');fig,ax=plt.subplots(figsize=(7,4.5));e=d.groupby('n').err.median();ax.plot(e.index,e.values,marker='o',label='E9');o=pd.read_csv(O1R)
 for est in sorted(o.estimator.unique()):
  x=o[o.estimator==est].groupby('n').err.median();ax.plot(x.index,x.values,alpha=.45,ls='--',label=est)
 ax.axhline(.10,color='C3',ls=':',label='0.10');ax.set_xscale('symlog',linthresh=500);ax.set_xlabel('cells per copy');ax.set_ylabel('median Test R error');ax.legend(ncol=3,fontsize=7);fig.tight_layout();fig.savefig(OUT/'fig_recovery.png',dpi=180);plt.close(fig)
 files=['requirements.lock','testR.csv','testF.csv','Q_hat_F.npy','views.csv','fig_umap.png','fig_recovery.png','scripts/o1f_run.py'];np.save(OUT/'Q_hat_F.npy',np.load(f_result(),allow_pickle=False)['Q_hat']);cmd='cd "$HOME/Documents/scopalign-manus-O1f" && .venv/bin/python reports/O1f/scripts/o1f_run.py';reg={p:cmd+' --stage all' for p in files};reg.update({f'err_R_E9_{n}':cmd+' --stage testR' for n in NS});reg.update({'err_F_E9':cmd+' --stage testF','frac_reached_truth':cmd+' --stage testR','n_alt_minimum':cmd+' --stage testR','e9_outcome':cmd+' --stage report','min_n_e9':cmd+' --stage report'})
 inputs=list(CASE_INPUTS)+[TESTF,SOURCE,HALVES,LAT1,LAT2,O1R];rep={'task_id':'O1f','executor':'manus','branch':'task/O1f','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':'2026-10-15T23:56:00-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':g,'deliverables':[{'path':p,'sha256':sha(OUT/p)} for p in files],'regeneration':reg,'deviations':[],'deviations_requested':[],'observations':[f'Input SHA-256: {p.relative_to(ROOT)} = {sha(p)}.' for p in inputs],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')};(OUT/'report.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n');print(json.dumps(g,sort_keys=True))
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('testR','testF','report','all'));a=p.parse_args().stage
 if a in ('testR','all'):testR()
 if a in ('testF','all'):testF()
 if a in ('report','all'):report()
if __name__=='__main__':main()
