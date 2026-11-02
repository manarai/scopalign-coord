#!/usr/bin/env python
"""O1m: fixed axis/sign search at d=30, synthetic bridge and PBMC transfer."""
import argparse, hashlib, importlib.util, json, os, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.linalg import expm
from scipy.stats import ortho_group
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from scjdo2.energy_drift import EnergyDriftAdapter, _EnergyMLP

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "O1m"
CACHE, ART = OUT / "cache", OUT / "artifacts"
O1B = ROOT.parent / "scopalign-manus-O1b" / "reports" / "O1b"
O1G_DIR = ROOT.parent / "scopalign-manus-O1g" / "reports" / "O1g"
O1J_DIR = ROOT.parent / "scopalign-manus-O1j" / "reports" / "O1j"
O1K_DIR = ROOT.parent / "scopalign-manus-O1k" / "reports" / "O1k"
O1L_DIR = ROOT.parent / "scopalign-manus-O1l" / "reports" / "O1l"
O1G_SCRIPT = O1G_DIR / "scripts" / "o1g_run.py"
O1L_SCRIPT = O1L_DIR / "scripts" / "o1l_run.py"
D, K, N = 30, 4, 3000
GAMMAS, REPS = (1.0, 1.05, 1.15), (0, 1, 2)
FIRST5 = (20261014, 20261016, 20261017, 20261018, 20261019)
KNOWN = [
 "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
 "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
 "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
 "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
GATES = ("symmetry_max","fits_match","l_star_repro_max","negctrl_err_shape_d30","errA_d30_g105","errA_d30_g115","errB_d30_g105","errB_d30_g115","err_R_2759","err_R_top5_2759","err_R_top10_2759","frac_hat_below_star_R","basin_radius_R","mix_F","mix_F_noA","sil_F","sil_F_noA","frame_ratio_F","o1m_outcome")
SYNTH_COLS = ("rung","gamma","n","r","err","err_shape","L_hat","L_star","best_start","screen_rank_best","sign_sweeps")
WORKER_PROBLEM = None
MAIN_THREADS = torch.get_num_threads()

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def atomic_path(p, suffix=".tmp"):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);f=tempfile.NamedTemporaryFile(dir=p.parent,prefix="."+p.name,suffix=suffix,delete=False);f.close();return Path(f.name)
def atomic_json(p,v):
    q=atomic_path(p);q.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n");os.replace(q,p)
def atomic_csv(p,d):
    q=atomic_path(p);d.to_csv(q,index=False);os.replace(q,p)
def atomic_npz(p,**x):
    q=atomic_path(p,".npz");np.savez(q,**x);os.replace(q,p)
def atomic_torch(p,x):
    q=atomic_path(p,".pt");torch.save(x,q);os.replace(q,p)
def make_dirs():
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(exist_ok=True);ART.mkdir(exist_ok=True)
    (OUT/".gitignore").write_text("cache/\n__pycache__/\nrun.log\nrun.pid\n.mplconfig/\n.numba/\n")
    meta=OUT/"run_meta.json"
    if not meta.exists(): atomic_json(meta,{"started":datetime.now(timezone.utc).isoformat()})
    return json.loads(meta.read_text())
def imp(name,path):
    spec=importlib.util.spec_from_file_location(name,Path(path).resolve());
    if spec is None or spec.loader is None: raise RuntimeError(f"Cannot import {path}")
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def o1g(): return imp("o1m_o1g",O1G_SCRIPT)
def o1l(): return imp("o1m_o1l",O1L_SCRIPT)
def stag(msg,**kw): print(json.dumps({"event":msg,**kw},sort_keys=True),flush=True)
def tag(g): return f"g{g:.2f}".replace(".","p")
def sid(g,r): return f"S_{tag(g)}_n{N}_r{r}"

# Exact 30-d O1k construction, with O1m's issued d=30 spectrum.
def landscape():
    rng=np.random.default_rng(20264200);means=np.zeros((K,D),dtype=np.float64)
    for k in range(K): means[k,:3]=3.0*rng.standard_normal(3)
    rots=np.stack([np.asarray(ortho_group.rvs(D,random_state=20264210+k),dtype=np.float64) for k in range(K)])
    qs=np.stack([np.asarray(ortho_group.rvs(D,random_state=20264100+r),dtype=np.float64) for r in REPS])
    return {"means":means,"rotations":rots,"q_stars":qs}
def sigma2(g): return np.asarray(g**(15.5-np.arange(1,D+1)),dtype=np.float64)
def sample(g,r):
    l=landscape();rng=np.random.default_rng(20264300+int(r));s=np.sqrt(sigma2(g));q=l["q_stars"][int(r)]
    xl=rng.integers(0,K,N);xe=rng.standard_normal((N,D));yl=rng.integers(0,K,N);ye=rng.standard_normal((N,D))
    def draw(labels,eps): return l["means"][labels]+np.einsum("nij,nj->ni",l["rotations"][labels],eps*s)
    x=draw(xl,xe);yp=draw(yl,ye);return x,yp@q.T,xl,yl,q
class ExactEnergy(torch.nn.Module):
    def __init__(self,g,qbase=None):
        super().__init__()
        l=landscape();qbase=np.eye(D) if qbase is None else np.asarray(qbase,dtype=np.float64)
        prec=np.einsum("kij,j,klj->kil",l["rotations"],1.0/sigma2(g),l["rotations"]);prec=.5*(prec+np.swapaxes(prec,1,2))
        self.register_buffer("means",torch.as_tensor(l["means"],dtype=torch.float64));self.register_buffer("precision",torch.as_tensor(prec,dtype=torch.float64));self.register_buffer("qbase",torch.as_tensor(qbase,dtype=torch.float64))
        self.register_buffer("lognorm",torch.as_tensor(-.5*(D*np.log(2*np.pi)+np.sum(np.log(sigma2(g))))-np.log(K),dtype=torch.float64))
    def forward(self,x):
        z=x@self.qbase;delta=z[:,None,:]-self.means[None,:,:];q=torch.einsum("nki,kij,nkj->nk",delta,self.precision,delta);return -torch.logsumexp(self.lognorm-.5*q,dim=1)
def exact_models(g,q): return [ExactEnergy(g).eval()],[ExactEnergy(g,q).eval()]

# Fitting preserves O1g config/settings, only dimension and synthetic inputs differ.
def fit_path(case,side,seed): return CACHE/"drifts"/f"{case}_{side}_{seed}.pt"
def save_row(case,side,seed,s): return {"fit_id":f"{case}_{side}_{seed}","case":case,"side":side,"seed":int(seed),"median_asym":float(s["symmetry"]["median_asym"]),"max_asym":float(s["symmetry"]["max_asym"])}
def fit_worker(item):
    torch.set_num_threads(1);case,side,seed,z,path=item;path=Path(path)
    if path.exists(): return save_row(case,side,seed,torch.load(path,map_location="cpu",weights_only=False))
    O=o1g();drift=EnergyDriftAdapter(D,config=O.cfg(int(seed))).fit(np.asarray(z,dtype=np.float32),verbose=True);sym=drift.symmetry_check(np.asarray(z,dtype=np.float32),n_cells=100)
    s={"state_dict":drift._model.state_dict(),"config":asdict(drift.config),"symmetry":{k:float(v) for k,v in sym.items()}};atomic_torch(path,s);return save_row(case,side,seed,s)
def fits_stage():
    make_dirs();pending=[];rows=[]
    for g in (1.05,1.15):
      for r in REPS:
       case=sid(g,r);x,y,*_=sample(g,r)
       for side,z in (("X",x),("Y",y)):
        for seed in FIRST5:
         p=fit_path(case,side,seed)
         if p.exists(): rows.append(save_row(case,side,seed,torch.load(p,map_location="cpu",weights_only=False)))
         else: pending.append((case,side,seed,z,str(p)))
    if pending:
      with ProcessPoolExecutor(max_workers=3) as ex: rows+=list(ex.map(fit_worker,pending))
    pd.DataFrame(sorted(rows,key=lambda x:x["fit_id"])).to_csv(OUT/"fit_symmetry.csv",index=False);stag("fits_finished",rows=len(rows),new=len(pending))
def build(saved):
    m=_EnergyMLP(D,int(saved["config"]["hidden"]),int(saved["config"]["depth"]));m.load_state_dict(saved["state_dict"]);m=m.to(dtype=torch.float64,device="cpu").eval()
    for p in m.parameters():p.requires_grad_(False)
    return m
def synth_models(g,r):
    case=sid(g,r);mx=[];my=[]
    for seed in FIRST5:
      for side,b in (("X",mx),("Y",my)):
       p=fit_path(case,side,seed)
       if not p.exists(): raise RuntimeError(f"missing synthetic fit {p}")
       b.append(build(torch.load(p,map_location="cpu",weights_only=False)))
    return mx,my

# E5s-30 search: exact O1g Problem/Closure, O1j complement sweep, O1l sign screen/refine.
def problem(modelsx,modelsy,x,y,r):
    p=o1g().Problem(modelsx,modelsy,x,y,int(r));p.d=D;return p
def axes(z):
    val,u=np.linalg.eigh(np.cov(np.asarray(z,dtype=np.float64),rowvar=False));return val[::-1],u[:,::-1]
def qaxis(ux,uy,top,comp): return uy@np.diag(np.concatenate([np.asarray(top,dtype=np.float64),np.asarray(comp,dtype=np.float64)]))@ux.T
def flow(p,q): return float(p.loss(torch.as_tensor(q,dtype=torch.float64),False)[1].detach())
def complement(p,ux,uy):
    signs=np.ones(D-10,dtype=np.float64);top=np.ones(10,dtype=np.float64);used=0
    for sweep in range(1,6):
      changed=False
      for loc in range(D-10):
       a=signs.copy();b=signs.copy();a[loc]=-1.;b[loc]=1.;choice=-1. if flow(p,qaxis(ux,uy,top,a))<flow(p,qaxis(ux,uy,top,b)) else 1.
       if choice!=signs[loc]: signs[loc]=choice;changed=True
      used=sweep
      if not changed:break
    return signs,used
def screen_refine(p,ux,uy):
    comp,sw=complement(p,ux,uy);qs=[]
    for index in range(1024):
      top=np.asarray([1. if ((index>>bit)&1)==0 else -1. for bit in range(10)],dtype=np.float64);qs.append(qaxis(ux,uy,top,comp))
    losses=o1l().screen_losses(p,np.stack(qs),q_batch=2)
    order=np.asarray(sorted(range(1024),key=lambda i:(float(losses[i]),int(i))),dtype=np.int64)[:8]
    best,answers=o1l().refine(p,[(int(i),qs[int(i)]) for i in order])
    full_order=sorted(range(1024),key=lambda i:(float(losses[i]),int(i)))
    rank=int(full_order.index(int(best["start"])))+1
    return best,np.asarray(losses,dtype=np.float64),comp,sw,rank

def case_path(kind,key): return CACHE/kind/f"{key}.npz"
def complete(p,keys):
    try:
      with np.load(p,allow_pickle=False) as d:return all(k in d.files for k in keys)
    except Exception:return False
def savecase(p,row,q,**extra): atomic_npz(p,Q=np.asarray(q,dtype=np.float64),**{k:np.asarray(v) for k,v in row.items()},**{k:np.asarray(v) for k,v in extra.items()})
def loadcase(p,cols):
    with np.load(p,allow_pickle=False) as d:return {k:np.asarray(d[k]).item() for k in cols},np.asarray(d["Q"],dtype=np.float64)
def synth_one(rung,g,r):
    key=f"{rung}_{sid(g,r)}";path=case_path("synth",key)
    if complete(path,SYNTH_COLS):return False
    x,y,_,_,qstar=sample(g,r);mx,my=exact_models(g,qstar) if rung=="A" else synth_models(g,r);p=problem(mx,my,x,y,r);_,ux=axes(x);_,uy=axes(y);best,losses,comp,sw,rank=screen_refine(p,ux,uy);q=np.asarray(best["Q"],dtype=np.float64);ls=float(p.loss(torch.as_tensor(qstar,dtype=torch.float64),False)[0].detach());delta=q-qstar
    shape=np.diag(np.r_[np.zeros(3),np.ones(D-3)])
    row={"rung":rung,"gamma":float(g),"n":N,"r":int(r),"err":float(np.linalg.norm(delta,"fro")/np.sqrt(2*D)),"err_shape":float(np.linalg.norm(delta@shape,"fro")/np.sqrt(2*(D-3))),"L_hat":float(best["L_hat"]),"L_star":ls,"best_start":int(best["start"]),"screen_rank_best":rank,"sign_sweeps":int(sw)}
    savecase(path,row,q,screen_losses=losses,complement_signs=comp);stag("synth_case",**{k:row[k] for k in ("rung","gamma","r","err","err_shape","L_hat","L_star")});return True
def synth_stage():
    make_dirs()
    for g in GAMMAS:
      for r in REPS:synth_one("A",g,r)
    for g in (1.05,1.15):
      for r in REPS:synth_one("B",g,r)
    rows=[]
    for rung,gs in (("A",GAMMAS),("B",(1.05,1.15))):
      for g in gs:
       for r in REPS:
        p=case_path("synth",f"{rung}_{sid(g,r)}")
        row,_=loadcase(p,SYNTH_COLS);rows.append(row)
    atomic_csv(OUT/"synth_cases.csv",pd.DataFrame(rows,columns=SYNTH_COLS).sort_values(["rung","gamma","r"],kind="stable"))

def manifest():
    p=O1J_DIR/"fits_used.csv"
    if not p.exists(): raise FileNotFoundError(p)
    return pd.read_csv(p)
def models_manifest(case):
    f=manifest();mx=[];my=[];used=[]
    for seed in tuple(o1g().SEEDS20):
      for side,b in (("X",mx),("Y",my)):
       q=f.loc[(f.case.astype(str)==str(case))&(f.side.astype(str)==side)&(f.seed.astype(int)==int(seed))]
       if len(q)!=1:raise RuntimeError(f"manifest missing {case} {side} {seed}")
       path=Path(q.iloc[0].path);want=str(q.iloc[0].sha256)
       if not path.is_file() or sha(path)!=want:raise RuntimeError(f"manifest hash mismatch {path}")
       b.append(build(torch.load(path,map_location="cpu",weights_only=False)));used.append(path)
    return mx,my,used
def check_stage():
    make_dirs();f=manifest();needed=[]
    for case in [f"R_n2759_r{r}" for r in REPS]+["F"]:
      for seed in tuple(o1g().SEEDS20):
       for side in ("X","Y"):
        q=f.loc[(f.case.astype(str)==case)&(f.side.astype(str)==side)&(f.seed.astype(int)==int(seed))]
        if len(q)!=1:raise RuntimeError(f"O1j manifest cardinality {case} {side} {seed}")
        p=Path(q.iloc[0].path);needed.append({"case":case,"side":side,"seed":int(seed),"path":str(p),"sha256":sha(p) if p.exists() else None,"expected":str(q.iloc[0].sha256)})
    match=int(all(x["sha256"]==x["expected"] for x in needed))
    syn=[]
    for p in (O1G_SCRIPT,O1L_SCRIPT,O1J_DIR/"fits_used.csv",O1K_DIR/"scripts"/"o1k_run.py"):
      syn.append({"path":str(p),"sha256":sha(p)})
    atomic_json(ART/"inputs.json",{"fits_match":match,"manifest":str(O1J_DIR/"fits_used.csv"),"fits":needed,"reference_scripts":syn})
    stag("check_finished",fits_match=match,fit_count=len(needed))
def pbmc_r(r):
    O=o1g();z=np.load(O.LAT1).astype(np.float32);perm=np.random.default_rng(20261030+int(r)).permutation(5518);q=np.asarray(ortho_group.rvs(D,random_state=20261031+int(r)),dtype=np.float64);x=z[np.sort(perm[:2759])];y=(z[np.sort(perm[2759:5518])].astype(np.float64)@q.T).astype(np.float32);return x,y,q
def pbmc_one(r):
    key=f"R_n2759_r{r}";path=case_path("pbmc",key);cols=("n","replicate","err","err_5","err_10","L_hat","L_star","best_start","screen_rank_best","sign_sweeps")
    if complete(path,cols):return False
    x,y,qstar=pbmc_r(r);mx,my,_=models_manifest(key);p=problem(mx,my,x,y,r);_,ux=axes(x);_,uy=axes(y);best,losses,comp,sw,rank=screen_refine(p,ux,uy);q=np.asarray(best["Q"],dtype=np.float64);ordering=np.argsort(-np.var(x.astype(np.float64),axis=0),kind="stable");delta=q-qstar;ls=float(p.loss(torch.as_tensor(qstar,dtype=torch.float64),False)[0].detach())
    row={"n":2759,"replicate":int(r),"err":float(np.linalg.norm(delta,"fro")/np.sqrt(60)),"err_5":float(np.linalg.norm(delta[:,ordering[:5]],"fro")/np.sqrt(10)),"err_10":float(np.linalg.norm(delta[:,ordering[:10]],"fro")/np.sqrt(20)),"L_hat":float(best["L_hat"]),"L_star":ls,"best_start":int(best["start"]),"screen_rank_best":rank,"sign_sweeps":int(sw)}
    savecase(path,row,q,screen_losses=losses,complement_signs=comp);stag("pbmc_R_case",**row);return True
def f_ops():
    O=o1g();out={}
    for side in ("H1","H2"):
      for seed in FIRST5:
       p=O.CACHE/"operators"/f"F_{side}_{seed}.npy"
       if not p.exists():raise FileNotFoundError(f"Missing read-only O1g operator cache {p}")
       out[(side,seed)]=np.load(p)
    return out
def umap(rep): return o1g().umap(rep)
def mix(rep,ds): return o1g().mixing(rep,ds)
def fig_umap(x,y,qhat,qload,h1,h2,ops):
    O=o1g();q_e5b=np.load(O1G_DIR/"Q_hat_F.npy").astype(np.float64);l1=h1.obs["cell_type"].astype(str).to_numpy();l2=h2.obs["cell_type"].astype(str).to_numpy();ds=np.asarray(["H1"]*len(x)+["H2"]*len(y));types=np.r_[l1,l2];tested=sorted(t for t in set(l1)&set(l2) if (l1==t).sum()>=50 and (l2==t).sum()>=50);mask=np.isin(types,tested);m1={s:ops[("H1",s)] for s in FIRST5};m2={s:ops[("H2",s)] for s in FIRST5}
    reps=[("no_rotation",np.eye(D)),("E5b",q_e5b),("E5s_30",qhat),("known_answer",qload)];rows=[];coords=[]
    for name,q in reps:
      rep=np.vstack([x@q.T,y]);rows.append({"view":name,"mix":mix(rep,ds),"sil":float(silhouette_score(rep[mask],types[mask],metric="euclidean")),"frame_ratio":float(O.frame_ratio(x.astype(np.float64),y.astype(np.float64),l1,l2,m1,m2,q))});coords.append(umap(rep))
    atomic_csv(OUT/"testF.csv",pd.DataFrame(rows));fig,ax=plt.subplots(4,2,figsize=(10,17));colors=plt.get_cmap("tab10")
    for i,(row,xy) in enumerate(zip(rows,coords)):
      for label,c in (("H1","C0"),("H2","C1")):
       m=ds==label;ax[i,0].scatter(xy[m,0],xy[m,1],s=4,alpha=.7,color=c,label=label,linewidths=0)
      for j,t in enumerate(sorted(set(types))):
       m=types==t;ax[i,1].scatter(xy[m,0],xy[m,1],s=4,alpha=.7,color=colors(j%10),label=t,linewidths=0)
      for a in ax[i]:a.set_xticks([]);a.set_yticks([]);a.legend(markerscale=2,fontsize=7,ncol=2)
      ax[i,0].set_title(f"{row['view']}: mix={row['mix']:.3f}; sil={row['sil']:.3f}; frame={row['frame_ratio']:.3f}");ax[i,1].set_title(row['view'])
    fig.tight_layout();fig.savefig(OUT/"fig_umap.png",dpi=180);plt.close(fig);return {x["view"]:x for x in rows}
def pbmc_f():
    p=case_path("pbmc","F");cols=("mix_F","mix_F_noA","sil_F","sil_F_noA","frame_ratio_F","L_hat","best_start")
    if complete(p,cols):return False
    O=o1g();x,y,qload,chk,h1,h2=O.pca_halves();mx,my,_=models_manifest("F");pr=problem(mx,my,x,y,0);_,ux=axes(x);_,uy=axes(y);best,losses,comp,sw,rank=screen_refine(pr,ux,uy);q=np.asarray(best["Q"],dtype=np.float64);views=fig_umap(x,y,q,qload,h1,h2,f_ops());row={"mix_F":float(views["E5s_30"]["mix"]),"mix_F_noA":float(views["no_rotation"]["mix"]),"sil_F":float(views["E5s_30"]["sil"]),"sil_F_noA":float(views["no_rotation"]["sil"]),"frame_ratio_F":float(views["E5s_30"]["frame_ratio"]),"L_hat":float(best["L_hat"]),"best_start":int(best["start"]),"screen_rank_best":rank,"sign_sweeps":int(sw),"latent_check_max":float(chk)}
    savecase(p,row,q,screen_losses=losses,complement_signs=comp);np.save(ART/"Q_hat_F.npy",q);stag("pbmc_F",**row);return True
def write_pbmc():
    rows=[];qs=[];cols=("n","replicate","err","err_5","err_10","L_hat","L_star","best_start","screen_rank_best","sign_sweeps")
    for r in REPS:
      a,q=loadcase(case_path("pbmc",f"R_n2759_r{r}"),cols);rows.append(a);qs.append(q)
    atomic_csv(OUT/"testR.csv",pd.DataFrame(rows).sort_values("replicate"));atomic_npz(ART/"Q_hats.npz",Q_hats=np.stack(qs),reps=np.asarray(REPS))
def pbmc_stage():
    make_dirs()
    for r in REPS:pbmc_one(r)
    write_pbmc();pbmc_f()

def basin_one(r):
    key=f"R_n2759_r{r}";pout=case_path("basin",key)
    if complete(pout,("i","t","theta","final_err","final_L","returned")):return False
    x,y,qstar=pbmc_r(r);mx,my,_=models_manifest(key);pr=problem(mx,my,x,y,r);rows=[]
    for i in range(1,6):
      a=np.random.default_rng(20264400+i).standard_normal((D,D));a=a-a.T;a=a/np.linalg.norm(a,"fro")
      for target in (0.05,0.10,0.20,0.40):
       lo,hi=0.,1.
       dist=lambda z:float(np.linalg.norm(qstar@(expm(z*a)-np.eye(D)),"fro")/np.sqrt(60))
       while dist(hi)<target:hi*=2
       while hi-lo>1e-6:
        mid=(lo+hi)/2
        if dist(mid)<target:lo=mid
        else:hi=mid
       theta=(lo+hi)/2;q0=qstar@expm(theta*a);L=o1l();L.WORKER_PROBLEM=pr;ans=L.optimise_start_inprocess((0,q0));err=float(np.linalg.norm(np.asarray(ans["Q"])-qstar,"fro")/np.sqrt(60));final_l=float(ans["L_hat"]);rows.append((i,target,theta,err,final_l,int(err<=.10)))
    atomic_npz(pout,i=np.asarray([z[0] for z in rows]),t=np.asarray([z[1] for z in rows]),theta=np.asarray([z[2] for z in rows]),final_err=np.asarray([z[3] for z in rows]),final_L=np.asarray([z[4] for z in rows]),returned=np.asarray([z[5] for z in rows]));stag("basin_case",replicate=r);return True
def basin_stage():
    make_dirs();rows=[]
    for r in REPS:basin_one(r)
    for r in REPS:
      with np.load(case_path("basin",f"R_n2759_r{r}"),allow_pickle=False) as d:
       rows += [{"replicate":r,"i":int(i),"t":float(t),"theta":float(th),"final_err":float(e),"final_L":float(l),"returned":int(ret)} for i,t,th,e,l,ret in zip(d["i"],d["t"],d["theta"],d["final_err"],d["final_L"],d["returned"])]
    atomic_csv(OUT/"basin.csv",pd.DataFrame(rows).sort_values(["replicate","i","t"]))

def figures(s,b):
    fig,ax=plt.subplots(figsize=(7,4.5))
    for rung,c in (("A","C0"),("B","C1")):
      q=s[s.rung==rung]
      for g in sorted(q.gamma.unique()):
       d=q[np.isclose(q.gamma,g)];ax.scatter(np.full(len(d),g),d.err,color=c,alpha=.6);ax.plot([g],[d.err.median()],marker="o",color=c)
    ax.axhline(.1,color="black",ls=":");ax.set_xlabel("gamma");ax.set_ylabel("rotation error");fig.tight_layout();fig.savefig(OUT/"fig_ladder.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.5))
    for r in REPS:
      d=b[b.replicate==r];ax.plot([.05,.1,.2,.4],[d[np.isclose(d.t,t)].returned.mean() for t in (.05,.1,.2,.4)],marker="o",label=f"r={r}")
    ax.set_xlabel("start distance");ax.set_ylabel("return fraction");ax.set_ylim(-.03,1.03);ax.legend();fig.tight_layout();fig.savefig(OUT/"fig_basin.png",dpi=180);plt.close(fig)
def requirement():
    p=OUT/"requirements.lock";q=atomic_path(p)
    with q.open("w") as f:subprocess.run([sys.executable,"-m","pip","freeze"],check=True,stdout=f)
    os.replace(q,p)
def checks(paths):return [{"path":x,"sha256":sha(OUT/x)} for x in paths]
def checkpoint(stage):
    if os.environ.get("O1M_AUTO_CHECKPOINT","1")!="1":return
    subprocess.run(["git","-C",str(ROOT),"add","reports/O1m"],check=True)
    z=subprocess.run(["git","-C",str(ROOT),"diff","--cached","--quiet"]).returncode
    if z==0:return
    if z!=1:raise RuntimeError("git diff failure")
    subprocess.run(["git","-C",str(ROOT),"commit","-m",f"task: O1m {stage} stage"],check=True);subprocess.run(["git","-C",str(ROOT),"push","origin","task/O1m"],check=True)
def report_stage():
    make_dirs();s=pd.read_csv(OUT/"synth_cases.csv");tr=pd.read_csv(OUT/"testR.csv");tf=pd.read_csv(OUT/"testF.csv");te=tf.loc[tf["view"]=="E5s_30"].iloc[0];tn=tf.loc[tf["view"]=="no_rotation"].iloc[0];b=pd.read_csv(OUT/"basin.csv");fs=pd.read_csv(OUT/"fit_symmetry.csv");inp=json.loads((ART/"inputs.json").read_text());old=pd.read_csv(O1G_DIR/"testR.csv")
    lrep=max(abs(float(tr.loc[tr.replicate==r,"L_star"].iloc[0])-float(old.loc[(old.n==2759)&(old.replicate==r),"L_star"].iloc[0])) for r in REPS)
    neg=float(s[(s.rung=="A")&np.isclose(s.gamma,1.0)].err_shape.median());a105=float(s[(s.rung=="A")&np.isclose(s.gamma,1.05)].err.median());a115=float(s[(s.rung=="A")&np.isclose(s.gamma,1.15)].err.median());b105=float(s[(s.rung=="B")&np.isclose(s.gamma,1.05)].err.median());b115=float(s[(s.rung=="B")&np.isclose(s.gamma,1.15)].err.median())
    radii=[]
    for r in REPS:
      d=b[b.replicate==r];good=[t for t in (.05,.1,.2,.4) if int(d[np.isclose(d.t,t)].returned.sum())>=4];radii.append(max(good) if good else 0.)
    all_sym=max(float(fs.median_asym.max()),float(manifest().median_asym.max()))
    valid=all_sym<=1e-4 and inp["fits_match"]==1 and lrep<=1e-8 and neg>=.30
    rp=float(tr.err.median());fpass=float(te["mix"])>=.8 and float(te["frame_ratio"])<=1.5 and float(te["sil"])>=.9*float(tn["sil"])
    outcome="works" if rp<=.1 and fpass else "functional_only" if fpass else "search_d30" if a115>.1 else "dimension_limited" if b115>.1 else "data_limited"
    g={"symmetry_max":all_sym,"fits_match":int(inp["fits_match"]),"l_star_repro_max":lrep,"negctrl_err_shape_d30":neg,"errA_d30_g105":a105,"errA_d30_g115":a115,"errB_d30_g105":b105,"errB_d30_g115":b115,"err_R_2759":rp,"err_R_top5_2759":float(tr.err_5.median()),"err_R_top10_2759":float(tr.err_10.median()),"frac_hat_below_star_R":float(np.mean(tr.L_hat<=tr.L_star)),"basin_radius_R":float(np.median(radii)),"mix_F":float(te["mix"]),"mix_F_noA":float(tn["mix"]),"sil_F":float(te["sil"]),"sil_F_noA":float(tn["sil"]),"frame_ratio_F":float(te["frame_ratio"]),"o1m_outcome":outcome}
    if tuple(g)!=GATES:raise RuntimeError("gated keys")
    figures(s,b);requirement();m=json.loads((OUT/"run_meta.json").read_text());cmd='cd "$HOME/Documents/scopalign-manus-O1m" && .venv/bin/python reports/O1m/scripts/o1m_run.py';deliv=("requirements.lock","synth_cases.csv","testR.csv","testF.csv","basin.csv","fit_symmetry.csv","fig_ladder.png","fig_basin.png","fig_umap.png","scripts/o1m_run.py");arts=("artifacts/inputs.json","artifacts/screen.csv","artifacts/Q_hats.npz","artifacts/Q_hat_F.npy")
    # Consolidate all screen vectors only after every individual cache is complete.
    screen=[]
    for p in sorted((CACHE/"synth").glob("*.npz"))+sorted((CACHE/"pbmc").glob("*.npz")):
      with np.load(p,allow_pickle=False) as d:
       if "screen_losses" in d:screen.append({"case":p.stem,"n_screen":int(len(d["screen_losses"])),"min_loss":float(np.min(d["screen_losses"]))})
    atomic_csv(ART/"screen.csv",pd.DataFrame(screen))
    rep={"task_id":"O1m","executor":"manus","branch":"task/O1m","commit":subprocess.check_output(["git","-C",str(ROOT),"rev-parse","HEAD"],text=True).strip(),"seed":20261014,"started":m["started"],"finished":datetime.now(timezone.utc).isoformat(),"status":"complete","gated":g,"deliverables":checks(deliv),"artifacts":checks(arts),"regeneration":{x:cmd+" --stage all" for x in list(deliv)+list(arts)+list(GATES)},"deviations":[{"what":"Corrected task-local ExactEnergy base-module initialization before any synthetic case was written.","why":"The initial wrapper could not register its fixed buffers; completed 60 drift fits were retained and neither fitting configuration nor analysis inputs changed."},{"what":"Replaced a screen-rank output bookkeeping expression with deterministic Python lexicographic ordering by (loss, start index).","why":"The original NumPy indexing expression failed after the full screen and refinements completed; the prescribed ranking and selected starts are unchanged."},{"what":"Added final_L to each basin CSV/cache row before the basin stage.","why":"The pre-registered procedure requires recording each start’s final error and objective; no completed analysis cache is altered."},{"what":"Evaluated the 1024 full-loss screen candidates in bounded batches of two rather than one at a time.","why":"The one-at-a-time PBMC screen had no completed checkpoint after four hours. Each candidate remains evaluated under the same full float64 D+H loss, with the same starts, lexicographic ranking, and eight LBFGS refinements; batching only changes execution scheduling while respecting the local ensemble-memory limit."}],"deviations_requested":[],"observations":[],"known_failures":KNOWN,"lockfile_sha256":sha(OUT/"requirements.lock")}
    atomic_json(OUT/"report.json",rep);stag("report_finished",validity=bool(valid),gated=g)
def main():
    p=argparse.ArgumentParser();p.add_argument("--stage",required=True,choices=("check","fits","synth","pbmc","basin","report","all"));s=p.parse_args().stage
    for name,fn in (("check",check_stage),("fits",fits_stage),("synth",synth_stage),("pbmc",pbmc_stage),("basin",basin_stage),("report",report_stage)):
      if s in (name,"all"):
       fn()
       if s=="all":checkpoint(name)
if __name__=="__main__":main()
