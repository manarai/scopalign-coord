import importlib.util, numpy as np, torch, os
p=os.path.expanduser('~/Documents/scopalign-manus-O1g/reports/O1g/scripts/o1g_run.py')
s=importlib.util.spec_from_file_location('o1g',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
from scipy.stats import ortho_group
torch.set_num_threads(8)
r=0;n=2759;z=np.load(m.LAT1).astype(np.float32);rng=np.random.default_rng(20261030+r);perm=rng.permutation(5518)
qs=np.asarray(ortho_group.rvs(30,random_state=20261031+r));x=z[np.sort(perm[:n])].astype(np.float64);y=(z[np.sort(perm[n:2*n])].astype(np.float64)@qs.T)
mx,my=m.get_models(f'R_n{n}_r{r}')
def stats(models,c):
    t=torch.as_tensor(c);u,_=m.ensemble(models,t,False,False);u=u.numpy()
    sub=np.random.default_rng(0).choice(len(c),256,replace=False)
    _,J=m.ensemble(models,torch.as_tensor(c[sub]),True,False);J=J.numpy()
    U=u.T@u/len(u); JJ=np.einsum('nij,nkj->ik',J,J)/len(J); Jm=J.mean(0); C=np.cov(c.T)
    return {'U':U,'JJ':JJ,'Jmean':-Jm,'S':U/np.trace(U)+JJ/np.trace(JJ),'pos_cov':C}
SX,SY=stats(mx,x),stats(my,y)
def top(M,k): w,v=np.linalg.eigh((M+M.T)/2); i=np.argsort(-w); return v[:,i[:k]], w[i]
for name in SX:
    out=[]
    for k in (3,5,8):
        vx,wx=top(SX[name],k); vy,wy=top(SY[name],k)
        sv=np.linalg.svd(vy.T@(qs@vx),compute_uv=False); out.append((k, round(float(np.sqrt(max(0,1-sv.min()**2))),3), round(float(wx[k-1]/wx[k]),2)))
    print(name, '(k, sub_k, gap)', out, flush=True)
