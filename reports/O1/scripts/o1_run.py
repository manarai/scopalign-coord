#!/usr/bin/env python
import argparse, hashlib, json, subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.linalg import orthogonal_procrustes
from scipy.stats import ortho_group
from sklearn.decomposition import PCA
import torch
from scjdo2.energy_drift import EnergyDriftAdapter, EnergyDriftConfig, _EnergyMLP

ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'reports'/'O1'; CASES=OUT/'cases'; MODELS=OUT/'models'; DATA=ROOT/'data'/'pbmc_10k_v3'; K0=ROOT/'reports'/'K0'
SOURCE=DATA/'pbmc_10k_v3.h5ad'; HALVES=DATA/'k0_halves.csv'; GENES=K0/'genes.txt'; LAT1=K0/'latent_pca30_h1.npy'; LAT2=K0/'latent_pca30_h2.npy'
INPUTS=(SOURCE,HALVES,GENES,LAT1,LAT2); SEEDS=(20261014,20261016,20261017,20261018,20261019); NS=(500,1000,2000,2759); REPS=(0,1,2); ESTS=('E0','E1','E2','E3','E4'); D=30
KNOWN=['tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]','tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]','tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]','tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def gatesp(): return OUT/'gates.json'
def gates(): return json.loads(gatesp().read_text()) if gatesp().exists() else {}
def saveg(x): gatesp().write_text(json.dumps(x,sort_keys=True)+'\n')
def cfg(seed): return EnergyDriftConfig(hidden=96,depth=3,n_epochs=1500,seed=seed,batch_size=256,loss='esm',lr=5e-4,beta=1.,esm_laplacian='hutchinson',device='cpu')
def rawqc():
 a=sc.read_h5ad(SOURCE);return a[a.obs['qc_pass'].to_numpy(dtype=bool)].copy()
def halves():
 a=rawqc(); tab=pd.read_csv(HALVES).set_index('barcode')['half']; names=a.obs_names.astype(str); h1=names.isin(tab.index)&(names.map(tab)=='H1');h2=names.isin(tab.index)&(names.map(tab)=='H2');return a[h1].copy(),a[h2].copy()
def norm(a):
 a=a.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);return a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
def cpath(n,r):return CASES/f'case_{n}_{r}.npz'
def modelpath(tag,side,seed):return MODELS/f'{tag}_{side}_{seed}.pt'
def symm(j):return (j+np.swapaxes(j,1,2))/2.0
def store_sym(rows):
 p=OUT/'fit_symmetry.csv'; old=pd.read_csv(p) if p.exists() else pd.DataFrame();pd.concat([old,pd.DataFrame(rows)],ignore_index=True).to_csv(p,index=False)
def load_drift(path):
 x=torch.load(path,map_location='cpu',weights_only=False);d=EnergyDriftAdapter(D,config=EnergyDriftConfig(**x['config']));m=_EnergyMLP(D,d.config.hidden,d.config.depth).to('cpu');m.load_state_dict(x['state_dict']);m.eval();d._model,d._device,d._fitted=m,'cpu',True;return d
def fit_fields(z,tag,side):
 ms=[];us=[];rows=[]
 for seed in SEEDS:
  path=modelpath(tag,side,seed)
  d=load_drift(path) if path.exists() else EnergyDriftAdapter(D,config=cfg(seed)).fit(z.astype(np.float32),verbose=True)
  if not path.exists(): torch.save({'state_dict':d._model.state_dict(),'config':asdict(d.config)},path)
  ck=d.symmetry_check(z.astype(np.float32),n_cells=100)
  ms.append(symm(d.jacobian_autodiff(z.astype(np.float32)).astype(np.float64)));us.append(d.velocity(z.astype(np.float32)).astype(np.float64));rows.append({'case':tag,'side':side,'seed':seed,'median_asym':float(ck['median_asym']),'max_asym':float(ck['max_asym'])})
 return np.stack(ms),np.stack(us),rows
def group_indices(tx,ty):
 return [(t,np.flatnonzero(tx==t),np.flatnonzero(ty==t)) for t in sorted(set(tx)&set(ty)) if (tx==t).sum()>=20 and (ty==t).sum()>=20]
def as_t(x):return torch.as_tensor(x,dtype=torch.float64,device='cpu')
def qof(q0,s):return q0@torch.matrix_exp(s-s.T)
def loss_pair(q0,mx,my,ux=None,uy=None):
 mx=as_t(mx);my=as_t(my);ux=None if ux is None else as_t(ux);uy=None if uy is None else as_t(uy);denm=torch.sum(my*my);denu=None if uy is None else torch.sum(uy*uy)
 def f(q):
  qq=q@mx@q.T;h=torch.sum((qq-my)**2)/denm
  return h if ux is None else h+torch.sum((q@ux.T-uy.T)**2)/denu
 return f
def loss_type(q0,mx,my,w):
 mx=as_t(mx);my=as_t(my);w=as_t(w);den=torch.sum(w*torch.sum(my*my,dim=(1,2)))
 def f(q):return torch.sum(w*torch.sum((q@mx@q.T-my)**2,dim=(1,2)))/den
 return f
def minimise(q0,fun):
 s=torch.zeros((D,D),dtype=torch.float64,requires_grad=True);opt=torch.optim.LBFGS([s],lr=1,max_iter=500,tolerance_grad=1e-10,tolerance_change=1e-12,line_search_fn='strong_wolfe')
 def closure():opt.zero_grad();v=fun(qof(as_t(q0),s));v.backward();return v
 opt.step(closure)
 with torch.no_grad():q=qof(as_t(q0),s).cpu().numpy();v=float(fun(as_t(q)).cpu())
 return q,v
def e0(mx,my,tx,ty):
 gs=group_indices(tx,ty);ax=np.stack([mx[0,ix].mean(axis=0) for _,ix,_ in gs]);ay=np.stack([my[0,iy].mean(axis=0) for _,_,iy in gs]);w=np.array([min(len(ix),len(iy)) for _,ix,iy in gs],float);fun=loss_type(None,ax,ay,w);best=(None,np.inf)
 for i in range(20):
  q0=np.eye(D) if i==0 else ortho_group.rvs(D,random_state=20261021+i);q,v=minimise(q0,fun)
  if v<best[1]:best=(q,v)
 return best[0],1,best[1]
def pair_list(mx,my,tx,ty,initial,q=None,ux=None,uy=None):
 ii=[];jj=[]
 for _,ix,iy in group_indices(tx,ty):
  if initial:
   ex=np.linalg.eigvalsh(mx[ix]);ey=np.linalg.eigvalsh(my[iy]);c=np.sum(ex*ex,axis=1)[:,None]+np.sum(ey*ey,axis=1)[None,:]-2*ex@ey.T
  else:
   a=q@mx[ix]@q.T;af=a.reshape(len(ix),-1);bf=my[iy].reshape(len(iy),-1);c=(np.sum(af*af,axis=1)[:,None]+np.sum(bf*bf,axis=1)[None,:]-2*af@bf.T)
   if ux is not None:
    av=(q@ux[ix].T).T;bv=uy[iy];mm=np.median(np.sum(my[iy]*my[iy],axis=(1,2)));mu=np.median(np.sum(bv*bv,axis=1));c=c/mm+(np.sum(av*av,axis=1)[:,None]+np.sum(bv*bv,axis=1)[None,:]-2*av@bv.T)/mu
  p=ot.emd(ot.unif(len(ix)),ot.unif(len(iy)),c);ii.append(ix);jj.append(iy[np.argmax(p,axis=1)])
 return np.concatenate(ii),np.concatenate(jj)
def procrustes(ux,uy,ii,jj):
 u,_,vt=np.linalg.svd(uy[jj].T@ux[ii],full_matrices=False);return u@vt
def icp(mx,my,ux,uy,tx,ty,kind):
 base_mx,base_my=(mx[0],my[0]) if kind in ('E1','E2') else (mx.mean(axis=0),my.mean(axis=0));base_ux,base_uy=(ux[0],uy[0]) if kind=='E2' else (ux.mean(axis=0),uy.mean(axis=0))
 if kind=='E4':hmx,hmy=gapfield(base_mx),gapfield(base_my)
 else:hmx,hmy=base_mx,base_my
 ii,jj=pair_list(base_mx,base_my,tx,ty,True)
 q=None;last=np.column_stack([ii,jj]);final=np.nan
 for it in range(1,11):
  fun=loss_pair(None,hmx[ii],hmy[jj],None,None) if kind=='E1' else loss_pair(None,hmx[ii],hmy[jj],base_ux[ii],base_uy[jj])
  if q is None:
   if kind=='E1':
    best=(None,np.inf)
    for s in range(20):
     q0=np.eye(D) if s==0 else ortho_group.rvs(D,random_state=20261021+s);qq,v=minimise(q0,fun)
     if v<best[1]:best=(qq,v)
    q,final=best
   else:q,final=minimise(procrustes(base_ux,base_uy,ii,jj),fun)
  else:q,final=minimise(q,fun)
  ni,nj=pair_list(hmx,hmy,tx,ty,False,q,base_ux if kind!='E1' else None,base_uy if kind!='E1' else None)
  cur=np.column_stack([ni,nj])
  if np.array_equal(cur,last):return q,it,final
  ii,jj,last=ni,nj,cur
 return q,10,final
def gapfield(m):
 vals,vec=np.linalg.eigh(m);out=np.zeros_like(m)
 for i in range(len(m)):
  v=vals[i];gap=np.array([np.min(np.abs(v[k]-np.delete(v,k))) for k in range(D)]);g=np.minimum(1.,gap/(.1*np.abs(v)));idx=np.argsort(np.abs(v))[-10:];vv=vec[i][:,idx];out[i]=(vv*(g[idx]*v[idx]))@vv.T
 return out
def estimate(mx,my,ux,uy,tx,ty,e):
 return e0(mx,my,tx,ty) if e=='E0' else icp(mx,my,ux,uy,tx,ty,e)
def err(q,target):return float(np.linalg.norm(q-target,'fro')/np.sqrt(60.))
def write_case(n,r,z,types):
 rng=np.random.default_rng(20261030+r);perm=rng.permutation(5518);sx=np.sort(perm[:n]);sy=np.sort(perm[n:2*n]);qs=ortho_group.rvs(D,random_state=20261031+r);x=z[sx].astype(np.float32);y=(z[sy]@qs.T).astype(np.float32);mx,ux,a=fit_fields(x,f'{n}_{r}','X');my,uy,b=fit_fields(y,f'{n}_{r}','Y');np.savez_compressed(cpath(n,r),S_X=sx,S_Y=sy,Q_star=qs,M_X=mx.astype(np.float32),M_Y=my.astype(np.float32),u_X=ux.astype(np.float32),u_Y=uy.astype(np.float32),type_X=types[sx].astype(str),type_Y=types[sy].astype(str));store_sym(a+b)
def resultpath(n,r): return CASES/f'case_{n}_{r}_results.csv'
def testR():
 CASES.mkdir(parents=True,exist_ok=True);MODELS.mkdir(exist_ok=True);z=np.load(LAT1).astype(np.float64);h1,_=halves();types=h1.obs['cell_type'].astype(str).to_numpy();rows=[]
 for n in NS:
  for r in REPS:
   if not cpath(n,r).exists():write_case(n,r,z,types)
   rp=resultpath(n,r);done=pd.read_csv(rp) if rp.exists() else pd.DataFrame(columns=['n','replicate','estimator','err','n_iter','final_loss'])
   with np.load(cpath(n,r),allow_pickle=False) as x:mx=x['M_X'].astype(np.float64);my=x['M_Y'].astype(np.float64);ux=x['u_X'].astype(np.float64);uy=x['u_Y'].astype(np.float64);tx=x['type_X'].astype(str);ty=x['type_Y'].astype(str);target=x['Q_star'].astype(np.float64)
   for e in ESTS:
    if e in set(done.estimator.astype(str)): continue
    q,it,loss=estimate(mx,my,ux,uy,tx,ty,e);done=pd.concat([done,pd.DataFrame([{'n':n,'replicate':r,'estimator':e,'err':err(q,target),'n_iter':it,'final_loss':loss}])],ignore_index=True);done.to_csv(rp,index=False)
   rows.append(done)
 df=pd.concat(rows,ignore_index=True);df.to_csv(OUT/'testR.csv',index=False);g=gates()
 for e in ESTS:
  for n in NS:g[f'err_R_{e}_{n}']=float(np.median(df[(df.estimator==e)&(df.n==n)].err))
 q0=ortho_group.rvs(D,random_state=20261031);g['err_random']=float(np.median([err(ortho_group.rvs(D,random_state=20261100+i),q0) for i in range(100)]));sp=pd.read_csv(OUT/'fit_symmetry.csv');g['symmetry_max']=float(sp.median_asym.max());saveg(g);print(json.dumps({k:g[k] for k in g if k.startswith('err_R_') or k in ('err_random','symmetry_max')},sort_keys=True))
def testF():
 MODELS.mkdir(exist_ok=True);h1,h2=halves();genes=GENES.read_text().splitlines();p1=PCA(30,svd_solver='full',random_state=20261014);p2=PCA(30,svd_solver='full',random_state=20261014);z1=p1.fit_transform(norm(h1[:,genes]));z2=p2.fit_transform(norm(h2[:,genes]));old1=np.load(LAT1);old2=np.load(LAT2);check=float(max(np.abs(z1-old1).max(),np.abs(z2-old2).max()));u,_,vt=np.linalg.svd(p2.components_@p1.components_.T,full_matrices=False);ql=u@vt;mx,ux,a=fit_fields(z1.astype(np.float32),'F','X');my,uy,b=fit_fields(z2.astype(np.float32),'F','Y');store_sym(a+b);tx=h1.obs['cell_type'].astype(str).to_numpy();ty=h2.obs['cell_type'].astype(str).to_numpy();rows=[];pack={'Q_load':ql,'Mbar_X':mx.mean(axis=0).astype(np.float32),'Mbar_Y':my.mean(axis=0).astype(np.float32),'ubar_X':ux.mean(axis=0).astype(np.float32),'ubar_Y':uy.mean(axis=0).astype(np.float32)}
 for e in ESTS:
  q,it,loss=estimate(mx,my,ux,uy,tx,ty,e);pack[f'Q_hat_{e}']=q;rows.append({'estimator':e,'err':err(q,ql),'n_iter':it,'final_loss':loss})
 np.savez_compressed(OUT/'testF.npz',**pack);pd.DataFrame(rows).to_csv(OUT/'testF.csv',index=False);g=gates();g['latent_check_max']=check
 for row in rows:g[f"err_F_{row['estimator']}"]=float(row['err'])
 sp=pd.read_csv(OUT/'fit_symmetry.csv');g['symmetry_max']=float(sp.median_asym.max());saveg(g);print(json.dumps({k:g[k] for k in g if k.startswith('err_F_') or k in ('latent_check_max','symmetry_max')},sort_keys=True))
def report():
 g=gates();r=pd.read_csv(OUT/'testR.csv');f=pd.read_csv(OUT/'testF.csv').set_index('estimator');selected='none';mn=0
 for nreq in (1000,2759):
  for e in ESTS:
   if float(f.loc[e,'err'])<=.25 and float(g[f'err_R_{e}_{nreq}'])<=.10: selected=e;mn=next((n for n in NS if g[f'err_R_{e}_{n}']<=.10),0);break
  if selected!='none':break
 g['selected_estimator']=selected;g['min_n_selected']=int(mn);saveg(g)
 fig,ax=plt.subplots(figsize=(7,4.5))
 for e in ESTS:
  d=r[r.estimator==e].groupby('n',sort=True).err.median();ax.plot(d.index,d.values,marker='o',label=e)
 ax.axhline(.10,color='C3',ls='--',label='0.10');ax.axhline(g['err_random'],color='.4',ls=':',label='random median');ax.set_xscale('symlog',linthresh=500);ax.set_xlabel('cells per copy');ax.set_ylabel('median Test R error');ax.legend(ncol=2,fontsize=8);fig.tight_layout();fig.savefig(OUT/'fig_recovery.png',dpi=180);plt.close(fig)
 files=['requirements.lock','testR.csv','testF.csv','testF.npz','fig_recovery.png','scripts/o1_run.py'];casefiles=[str(p.relative_to(OUT)) for p in sorted(CASES.glob('case_*.npz'))];cmd='cd "$HOME/Documents/scopalign-manus" && .venv/bin/python reports/O1/scripts/o1_run.py';regen={p:cmd+' --stage all' for p in files}
 for e in ESTS:
  for n in NS:regen[f'err_R_{e}_{n}']=cmd+' --stage testR'
  regen[f'err_F_{e}']=cmd+' --stage testF'
 regen.update({'err_random':cmd+' --stage testR','symmetry_max':cmd+' --stage testF','latent_check_max':cmd+' --stage testF','selected_estimator':cmd+' --stage report','min_n_selected':cmd+' --stage report'})
 report={'task_id':'O1','executor':'manus','branch':'task/O1','commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'seed':20261014,'started':'2026-10-15T17:21:00-06:00','finished':datetime.now(timezone.utc).isoformat(),'status':'complete','gated':g,'deliverables':[{'path':p,'sha256':sha(OUT/p)} for p in files],'artifacts':[{'path':p,'sha256':sha(OUT/p)} for p in casefiles],'regeneration':regen,'deviations':[{'what':'Resumed a previously interrupted O1 run from the first missing Test R case (case_2759_1), preserving the ten existing case arrays.','why':'The prior run stopped after persisting 10 of 12 deterministic case arrays; the task owner directed continuation without regenerating saved cases.'},{'what':'Restarted the O1 pipeline as a detached caffeinate-protected process using --stage all.','why':'The foreground continuation ended before the second missing case and downstream stages completed; the task owner directed a resumable detached restart.'}],'deviations_requested':[],'observations':[f'Input SHA-256: {p.relative_to(ROOT)} = {sha(p)}.' for p in INPUTS],'known_failures':KNOWN,'lockfile_sha256':sha(OUT/'requirements.lock')}
 if g['symmetry_max']>1e-4 or g['latent_check_max']>1e-4:report['validity']='REVISE'
 (OUT/'report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(g,sort_keys=True))
def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=('testR','testF','report','all'));a=p.parse_args().stage
 if a in ('testR','all'):testR()
 if a in ('testF','all'):testF()
 if a in ('report','all'):report()
if __name__=='__main__':main()
