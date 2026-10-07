#!/usr/bin/env python3
"""Paired Gaussian benchmark; every newly generated asset remains in this folder.

Use a Python environment with the dependencies documented in README.md.
Run --check for saved-output verification or --run to compute missing fits. Input seeds,
initial supports, masses and oracle draws are serialized once for every case.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
BASE=Path(__file__).resolve().parent
# Plotting is handled by manuscript-figures; this computation has no font cache.
import argparse, csv, hashlib, json, platform, shutil, subprocess, time, warnings
try:
    import resource
except ImportError:
    resource=None
import numpy as np
import scipy
import ot
from scipy.special import xlogy
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_info, threadpool_limits
threadpool_limits(1)

N=300; NTRUTH=1000; CAP=300; INNER_CAP=20000
OBJ_TOL=1e-8; RES_TOL=1e-7; INNER_TOL=1e-9
REGS=(.1,.3,1.)


def sqrtm(A):
    val,vec=np.linalg.eigh((A+A.T)/2)
    return (vec*np.sqrt(np.maximum(val,0)))@vec.T
def bary_cov(covs):
    S=np.mean(covs,axis=0)
    for _ in range(1000):
        R=sqrtm(S); inv=np.linalg.inv(R)
        A=np.mean([sqrtm(R@C@R) for C in covs],axis=0)
        T=inv@A@A@inv
        if np.linalg.norm(S-T)<=1e-12*(1+np.linalg.norm(S)): return (T+T.T)/2
        S=(T+T.T)/2
    raise RuntimeError('Gaussian covariance reference did not converge')
def distances(X,Y): return np.maximum(np.sum(X*X,1)[:,None]+np.sum(Y*Y,1)[None,:]-2*X@Y.T,0)
def bures(A,B): return np.sqrt(max(float(np.trace(A)+np.trace(B)-2*np.trace(sqrtm(sqrtm(A)@B@sqrtm(A)))),0))
def write_rows(path,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def loadrow(path):
    with open(path) as f:return next(csv.DictReader(f))
def create_case(rep,m,start=0):
    cid=f'rep{rep:02d}_m{m:03d}'+(f'_start{start}' if start else '')
    d=BASE/'inputs'/cid;d.mkdir(parents=True,exist_ok=True)
    if (d/'problem.npz').exists():
        data=dict(np.load(d/'problem.npz'));return cid,data
    rng=np.random.default_rng(700000+rep)
    signs=np.array([[1,1],[1,-1],[-1,1],[-1,-1]])
    means=[];covs=[];measures=[]
    for sign in signs:
        mu=10*sign+rng.normal(size=2);A=rng.normal(size=(4,2));cov=A.T@A
        means.append(mu);covs.append(cov);measures.append(rng.multivariate_normal(mu,cov,size=N))
    means=np.array(means);covs=np.array(covs);measures=np.array(measures)
    S=bary_cov(covs);mu=means.mean(0)
    oracle=np.random.default_rng(900000+rep).multivariate_normal(mu,S,size=2000)
    init_seed=800000+1000*rep+m+10000000*start
    Z0=KMeans(n_clusters=m,n_init=5,max_iter=100,random_state=init_seed).fit(np.vstack(measures)).cluster_centers_
    scale2=float(np.mean([np.mean(np.sum((Y-Y.mean(0))**2,1)) for Y in measures]))
    data=dict(measures=measures,means=means,covs=covs,bary_mean=mu,bary_cov=S,oracle=oracle,Z0=Z0,
              input_masses=np.full((4,N),1/N),bary_masses=np.full(m,1/m),outer_weights=np.full(4,.25),scale2=scale2,
              data_seed=700000+rep,init_seed=init_seed,oracle_seed=900000+rep)
    np.savez_compressed(d/'problem.npz',**data)

    return cid,data

def export_r_inputs(cid,data):
    """CSV duplicates are generated only for an R rerun; NPZ is canonical."""
    folder=BASE/'inputs'/cid
    np.savetxt(folder/'init.csv',data['Z0'],delimiter=',')
    for i,Y in enumerate(data['measures'],1):
        np.savetxt(folder/f'measure{i}.csv',Y,delimiter=',')

def transport_state(X,data,ratio=0):
    b=data['bary_masses'];M=np.zeros_like(X);obj=0.;marginal=0.;n_inner=0;inner_ok=True
    for Y,a,w in zip(data['measures'],data['input_masses'],data['outer_weights']):
        C=distances(X,Y)
        if ratio:
            eps=ratio*float(data['scale2'])
            # Fixed row/column shifts preserve the regularized optimal coupling.
            Csolve=C-C.min(axis=1,keepdims=True);Csolve=Csolve-Csolve.min(axis=0,keepdims=True)
            with warnings.catch_warnings(record=True) as captured:
                P,log=ot.sinkhorn(b,a,Csolve,eps,method='sinkhorn_stabilized',numItermax=INNER_CAP,
                                  stopThr=INNER_TOL,warn=True,log=True)
            n_inner=max(n_inner,int(log.get('niter',log.get('n_iter',0))))
            inner_ok=inner_ok and not any('did not converge' in str(s.message) for s in captured)
            obj+=w*(np.sum(P*C)+eps*np.sum(xlogy(P,P)-P))
        else:
            P=ot.emd(b,a,C,numItermax=100000,numThreads=1)
            obj+=w*np.sum(P*C)
        if not np.isfinite(P).all():raise RuntimeError('nonfinite transport plan')
        marginal=max(marginal,float(np.max(np.abs(P.sum(1)-b))),float(np.max(np.abs(P.sum(0)-a))))
        M+=w*(P@Y)/b[:,None]
    residual=float(np.sum(b[:,None]*(M-X)**2)/data['scale2'])
    return M,float(obj),residual,marginal,n_inner,inner_ok

def run_python(cid,data,rep,m,method,ratio):
    (BASE/'runs').mkdir(exist_ok=True)
    out=BASE/'runs'/f'{cid}_{method}'
    if Path(str(out)+'.csv').exists():return
    X=data['Z0'].copy();previous=None;history=[];status='iteration_limit';all_inner=True
    t=time.perf_counter()
    for it in range(CAP+1):
        M,obj,res,err,ni,ok=transport_state(X,data,ratio)
        all_inner &= ok
        change=np.inf if previous is None else abs(obj-previous)/max(1,abs(obj),abs(previous))
        history.append(dict(iter=it,objective=obj,residual=res,relative_change=change,marginal_error=err,max_inner_iterations=ni,inner_converged=ok))
        if it and it%20==0:
            write_rows(str(out)+'_partial.csv',history)
            print(cid,method,'checkpoint',it,f'res={res:.2g}',flush=True)
        if change<=OBJ_TOL and res<=RES_TOL and err<=1e-8 and ok:
            status='converged';break
        if it==CAP:break
        previous=obj;X=M
    elapsed=time.perf_counter()-t
    np.savetxt(str(out)+'_support.csv',X,delimiter=',')
    write_rows(str(out)+'_history.csv',history)
    # External exact scoring is excluded from optimization runtime for all solvers.
    _,exact_obj,exact_res,_,_,_=transport_state(X,data)
    row=dict(id=cid,rep=rep,support=m,method=method,reg_ratio=ratio,runtime_sec=elapsed,niter=it,
             objective=exact_obj,residual=res,relative_change=change,marginal_error=err,status=status,
             scale2=float(data['scale2']),regularized_objective=obj if ratio else '',
             exact_residual=exact_res,all_inner_converged=all_inner,
             max_marginal_error=max(h['marginal_error'] for h in history),
             max_objective_increase=max([history[i]['objective']-history[i-1]['objective'] for i in range(1,len(history))],default=0))
    write_rows(str(out)+'.csv',[row])
    print(cid,method,it,status,f'{elapsed:.2f}s',f'res={res:.2g}',flush=True)

def evaluate(cid,data,method):
    out=BASE/'runs'/f'{cid}_{method}'
    row=loadrow(str(out)+'.csv')
    X=np.loadtxt(str(out)+'_support.csv',delimiter=',');b=data['bary_masses']
    center=b@X;cov=(X-center).T@(b[:,None]*(X-center))
    truth=data['oracle'][:NTRUTH]
    row.update(semi_w2=float(np.sqrt(ot.emd2(b,ot.unif(len(truth)),distances(X,truth)))),
               mean_error=float(np.linalg.norm(center-data['bary_mean'])),cov_error=float(bures(cov,data['bary_cov'])))
    if method=='pot_exact':
        native=ot.lp.free_support_barycenter(list(data['measures']),list(data['input_masses']),data['Z0'],b,
                weights=data['outer_weights'],numItermax=int(row['niter']),stopThr=-1.,numThreads=1)
        row['pot_native_max_difference']=float(np.max(np.abs(native-X)))
    else:row['pot_native_max_difference']=''
    # Uniform output schema.
    for key in ('regularized_objective','exact_residual','all_inner_converged'):
        row.setdefault(key,'')
    return row

def run_benchmark(assemble_only=False):
    reps=int(os.environ.get('REPS','10'));supports=[int(s) for s in os.environ.get('SUPPORTS','50,100,200').split(',')]
    jobs=[];cases=[]
    for rep in range(1,reps+1):
        for m in supports:
            cid,data=create_case(rep,m);jobs.append((cid,data,rep,m));cases.append(dict(id=cid,rep=rep,support=m,scale2=float(data['scale2'])))
    write_rows(BASE/'cases.csv',cases)
    if not assemble_only:
        for cid,data,_,_ in jobs: export_r_inputs(cid,data)
        subprocess.run(['Rscript',str(BASE/'exact_r.R'),str(BASE)],check=True,env=os.environ.copy())
    allrows=[]
    for cid,data,rep,m in jobs:
        if not assemble_only:
            for method,ratio in [('pot_exact',0)]+[(f'sinkhorn_{r:g}',r) for r in REGS]:
                run_python(cid,data,rep,m,method,ratio)
        for method in ['r_exact','pot_exact']+[f'sinkhorn_{r:g}' for r in REGS]:allrows.append(evaluate(cid,data,method))
        keys=sorted(set().union(*(r.keys() for r in allrows)))
        write_rows(BASE/'paired_results.csv',[{k:r.get(k,'') for k in keys} for r in allrows])
    hardware={}
    if platform.system()=='Darwin' and shutil.which('system_profiler'):
        try:
            raw=subprocess.run(['system_profiler','SPHardwareDataType','-json'],capture_output=True,text=True,check=True)
            h=json.loads(raw.stdout)['SPHardwareDataType'][0]
            hardware={k:h[k] for k in ('chip_type','physical_memory','number_processors','machine_model') if k in h}
        except (OSError,subprocess.CalledProcessError,ValueError,KeyError):pass
    pools=[{k:v for k,v in item.items() if k!='filepath'} for item in threadpool_info()]
    metadata=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,POT=ot.__version__,hardware=hardware,
      machine=platform.machine(),platform=platform.platform(),processor=platform.processor(),
      threads={k:os.environ[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')},
      threadpools=pools,n_inputs=4,n_per_input=N,n_oracle=NTRUTH,outer_iteration_cap=CAP,inner_iteration_cap=INNER_CAP,
      objective_tol=OBJ_TOL,normalized_residual_tol=RES_TOL,inner_tol=INNER_TOL,regularization_ratios=REGS,repetitions=reps,supports=supports,
      time_scope='Optimization plus convergence diagnostics; input generation, initialization, external exact scoring, and oracle scoring excluded.')
    if resource is not None:
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        metadata['process_peak_rss_bytes']=rss if platform.system()=='Darwin' else rss*1024
    (BASE/'rerun_metadata.json').write_text(json.dumps(metadata,indent=2))
    (BASE/'input-checksums.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(BASE))+'\n' for p in sorted((BASE/'inputs').glob('*/problem.npz'))))
    print('Completed paired results',len(allrows),flush=True)
def check_saved():
    with (BASE/'paired_results.csv').open() as f: rows=list(csv.DictReader(f))
    with (BASE/'multistart_results.csv').open() as f: multi=list(csv.DictReader(f))
    assert len(rows)==150 and len(multi)==15
    lookup={(int(r['rep']),int(r['support']),r['method']):r for r in rows}
    assert len(lookup)==150
    assert all(float(r['max_marginal_error'])<1e-8 for r in rows)
    assert all(r['all_inner_converged'].lower()=='true' for r in rows if r['method'].startswith('sinkhorn'))
    exact_diff=max(abs(float(lookup[(r,m,'r_exact')]['objective'])-float(lookup[(r,m,'pot_exact')]['objective'])) for r in range(1,11) for m in (50,100,200))
    assert exact_diff<1e-10
    for line in (BASE/'input-checksums.sha256').read_text().splitlines():
        checksum,rel=line.split('  ',1)
        assert hashlib.sha256((BASE/rel).read_bytes()).hexdigest()==checksum, rel
    # Recompute one case per support size from saved supports, without fitting.
    rescored=[]
    for m in (50,100,200):
        cid=f'rep01_m{m:03d}';data=dict(np.load(BASE/'inputs'/cid/'problem.npz'))
        row=evaluate(cid,data,'r_exact');saved=lookup[(1,m,'r_exact')]
        for key in ('semi_w2','mean_error','cov_error'):
            assert abs(float(row[key])-float(saved[key]))<1e-9,(cid,key)
        rescored.append(cid)
    report=dict(paired_fits=len(rows),multistart_fits=len(multi),input_checksums='verified',
                exact_objective_max_difference=exact_diff,rescored_saved_cases=rescored,
                all_inner_solves_converged=True,maximum_marginal_error=max(float(r['max_marginal_error']) for r in rows))
    (BASE/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--run',action='store_true',help='Compute missing fits and rescore all fits (may take several minutes).')
    modes.add_argument('--assemble-only',action='store_true',help='Rescore existing supports without fitting.')
    modes.add_argument('--check',action='store_true',help='Verify saved outputs and three representative accuracy calculations (default).')
    args=parser.parse_args()
    if args.run or args.assemble_only:run_benchmark(assemble_only=args.assemble_only)
    else:check_saved()
