"""Three predeclared H2 tests. Run only in the authorized hosted Linux job."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import base64,hashlib,io,json,platform,resource,signal,subprocess,sys,time,traceback,zlib
from contextlib import redirect_stdout
from pathlib import Path
START=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3))
resource.setrlimit(resource.RLIMIT_CPU,(720,730))
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
def deadline(*args): raise TimeoutError('720-second scientific deadline; no rerun')
signal.signal(signal.SIGALRM,deadline);signal.alarm(720)
import numpy as np
import scipy,pyscf,ase,h5py
from scipy.linalg import eigh
from pyscf import gto,scf,ao2mo,fci,lib
from ase import Atoms
from ase.units import Bohr
from threadpoolctl import threadpool_limits,threadpool_info
OUT=Path('c2/results');OUT.mkdir(parents=True,exist_ok=True)
P=json.loads(Path('c2/protocol.json').read_text())
REF=json.loads(Path('reference.json').read_text())
lib.num_threads(1)
def dump(name,data):
    (OUT/name).write_text(json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+'\n')
def scalar(x): return float(x)
def spin(ci,n): return scalar(fci.spin_op.spin_square(ci,n,(1,1))[0])
def mol_h2(R,basis):
    a=Atoms('H2',positions=[[0,0,-R*Bohr/2],[0,0,R*Bohr/2]])
    return gto.M(atom=[(s,tuple(p/Bohr)) for s,p in zip(a.get_chemical_symbols(),a.positions)],unit='Bohr',basis=basis,spin=0,verbose=0,max_memory=1600)
def integrals(m,C):
    h=C.T@(m.intor_symmetric('int1e_kin')+m.intor_symmetric('int1e_nuc'))@C
    n=C.shape[1]
    return h,ao2mo.kernel(m,C,compact=False).reshape((n,)*4)
def explicit_matrix(h,eri):
    n=len(h);I=np.eye(n)
    return np.kron(h,I)+np.kron(I,h)+eri.transpose(0,2,1,3).reshape(n*n,n*n)
def sector(H,n,sign):
    i,j=np.triu_indices(n,0 if sign==1 else 1);a=i*n+j;b=j*n+i
    weights=np.where(i==j,.5,1/np.sqrt(2.))
    block=(H[np.ix_(a,a)]+sign*H[np.ix_(a,b)]+sign*H[np.ix_(b,a)]+H[np.ix_(b,b)])*weights[:,None]*weights[None,:]
    return block,i,j,weights
def lift(v,i,j,w,n,sign=1):
    c=np.zeros((n,n));np.add.at(c,(i,j),v*w);np.add.at(c,(j,i),sign*v*w)
    return c
def residual(solver,h,eri,c,E,n,enuc):
    heff=solver.absorb_h1e(h,eri,n,(1,1),.5)
    return solver.contract_2e(heff,c,n,(1,1))-(E-enuc)*c
def test(label,basis,R):
    t=time.monotonic();m=mol_h2(R,basis);S=m.intor_symmetric('int1e_ovlp');s,U=eigh(S)
    if s[0]<=1e-10:raise ValueError('Do not silently truncate the basis')
    mf=scf.RHF(m);mf.conv_tol=1e-13;mf.max_cycle=100;mf.kernel()
    if not mf.converged:raise RuntimeError('RHF orbital generation failed')
    C=mf.mo_coeff;n=C.shape[1];enuc=m.energy_nuc();h,eri=integrals(m,C)
    H=explicit_matrix(h,eri)
    sing,i,j,w=sector(H,n,1);es,vs=eigh(sing,subset_by_index=[0,1],driver='evr',check_finite=True)
    dense_ci=lift(vs[:,0],i,j,w,n)
    trip,_,_,_=sector(H,n,-1);et=eigh(trip,subset_by_index=[0,0],driver='evr',eigvals_only=True)
    E_dense=scalar(es[0]+enuc)
    direct=fci.direct_spin1.FCI(m);direct.verbose=4
    direct.conv_tol=1e-14;direct.conv_tol_residual=1e-10;direct.lindep=1e-16
    direct.max_cycle=200;direct.max_space=80;direct.pspace_size=0;direct.davidson_only=True;direct.max_memory=1600
    rng=np.random.default_rng(20261015)
    c0=rng.normal(size=(n,n));c0=(c0+c0.T)*1e-3;c0[0,0]+=1;c0/=np.linalg.norm(c0)
    log=io.StringIO();direct.stdout=log
    E,ci=direct.kernel(h,eri,n,(1,1),ci0=c0,ecore=enuc)
    (OUT/f'{label}_solver.log').write_text(log.getvalue())
    r_direct=scalar(np.linalg.norm(residual(direct,h,eri,ci,E,n,enuc)))
    r_matrix=scalar(np.linalg.norm(H@ci.ravel()-(E-enuc)*ci.ravel()))
    rd=scalar(np.linalg.norm(H@dense_ci.ravel()-es[0]*dense_ci.ravel()))
    rd_direct=scalar(np.linalg.norm(residual(direct,h,eri,dense_ci,E_dense,n,enuc)))
    random=rng.normal(size=(n,n));random/=np.linalg.norm(random)
    eff=direct.absorb_h1e(h,eri,n,(1,1),.5)
    action_delta=scalar(np.linalg.norm(H@random.ravel()-direct.contract_2e(eff,random,n,(1,1)).ravel()))
    row={'test':label,'basis':basis,'R_bohr':R,'norb':n,'product_dimension':n*n,'singlet_dimension':len(es) if False else n*(n+1)//2,
         'iterative_energy_total_Eh':scalar(E),'dense_lowest_singlet_total_Eh':E_dense,
         'dense_lowest_triplet_total_Eh':scalar(et[0]+enuc),'singlet_next_total_Eh':scalar(es[1]+enuc),
         'solver_converged':bool(direct.converged),'iterative_S2':spin(ci,n),'dense_S2':spin(dense_ci,n),
         'iterative_residual_direct_Eh':r_direct,'iterative_residual_matrix_Eh':r_matrix,
         'dense_residual_matrix_Eh':rd,'dense_residual_direct_Eh':rd_direct,
         'energy_path_delta_Eh':scalar(E-E_dense),'arbitrary_vector_operator_delta_Eh':action_delta,
         'orthogonality_max':scalar(np.max(np.abs(C.T@S@C-np.eye(n)))),
         'H_hermiticity_max':scalar(np.max(np.abs(H-H.T))),
         'overlap_min_eigenvalue':scalar(s[0]),'nuclear_repulsion_Eh':enuc,
         'BO_reference_total_Eh':REF['points_Eh'][str(R)],'BO_difference_Eh':E_dense-REF['points_Eh'][str(R)],
         'integral_hashes':{'h_sha256':hashlib.sha256(h.tobytes()).hexdigest(),'eri_sha256':hashlib.sha256(eri.tobytes()).hexdigest()},
         'initial_dense_overlap_abs':scalar(abs(np.vdot(c0,dense_ci))),
         'iterative_dense_overlap_abs':scalar(abs(np.vdot(ci,dense_ci)))}
    if label=='A':
        atom=gto.M(atom='H 0 0 0',basis=basis,spin=1,verbose=0)
        ha=atom.intor_symmetric('int1e_kin')+atom.intor_symmetric('int1e_nuc')
        endpoint=2*scalar(eigh(ha,atom.intor_symmetric('int1e_ovlp'),eigvals_only=True)[0])
        full=eigh(H,eigvals_only=True)
        Co=U/np.sqrt(s)[None,:];ho,go=integrals(m,Co)
        old=fci.direct_spin0.FCI(m);old.conv_tol=1e-12;old.max_cycle=200;old.max_space=32;old.max_memory=1200
        eo,co=old.kernel(ho,go,n,(1,1),ecore=enuc)
        hd=old.make_hdiag(ho,go,n,(1,1),compress=True)
        addr,h0=old.pspace(ho,go,n,(1,1),hd,min(400,len(hd)))
        pw,pv=eigh(h0);projection=[]
        for k in range(len(pw)):
            c=np.empty(n*n);c[addr]=pv[:,k];c=c.reshape(n,n);norm=scalar(np.linalg.norm((c+c.T)*.5))
            projection.append({'index':k,'energy_total_Eh':scalar(pw[k]+enuc),'singlet_projection_norm':norm,'old_check_accepts':abs(norm-1)<=1e-6})
        accepted=next(q for q in projection if q['old_check_accepts'])
        dump('A_old_pspace_spectrum.json',projection)
        H_old=explicit_matrix(ho,go)
        row['original_path_reproduction']={'energy_total_Eh':scalar(eo),'converged':bool(old.converged),'residual_Eh':scalar(np.linalg.norm(residual(old,ho,go,co,eo,n,enuc))),
             'S2':spin(co,n),'first_accepted_pspace_root':accepted,'independent_old_representation_lowest_Eh':scalar(eigh(H_old,eigvals_only=True)[0]+enuc),
             'pspace_matrix_spectrum_vs_explicit_max_Eh':scalar(np.max(np.abs(pw-eigh(H_old,eigvals_only=True))))}
        row['same_basis_atom_endpoint_Eh']=endpoint;row['ground_minus_endpoint_Eh']=E_dense-endpoint
        row['full_product_lowest_total_Eh']=scalar(full[0]+enuc)
        row['full_product_vs_singlet_Eh']=scalar(full[0]+enuc-E_dense)
        dump('A_full_product_eigenvalues.json',[scalar(e+enuc) for e in full])
    gates={'converged':bool(direct.converged),
           'residuals':max(r_direct,r_matrix,rd,rd_direct)<=P['residual_gate_Eh'],
           'spin':max(abs(row['iterative_S2']),abs(row['dense_S2']))<=P['spin_S2_abs_max'],
           'energy_paths':abs(E-E_dense)<=P['energy_path_agreement_Eh'],
           'operator_action':action_delta<=P['operator_action_agreement_Eh'],
           'representation':max(row['orthogonality_max'],row['H_hermiticity_max'])<=P['orthogonality_and_hermiticity_max'],
           'lowest_spin_sector':E_dense<=et[0]+enuc+P['ground_vs_triplet_tolerance_Eh'],
           'variational_BO':row['BO_difference_Eh']>=-1e-8}
    if label=='A':
        gates['asymptote']=abs(row['ground_minus_endpoint_Eh'])<=P['A_endpoint_tolerance_Eh']
        gates['full_space_lowest']=abs(row['full_product_vs_singlet_Eh'])<=1e-9
    row['gates']=gates;row['numerical_pass']=all(gates.values());row['wall_seconds']=time.monotonic()-t
    dump(f'{label}_eigenvectors.json',{'iterative_CI':ci.tolist(),'dense_CI':dense_ci.tolist(),'representation':'RHF MO; alpha row beta column'})
    dump(f'{label}_result.json',row)
    print('C2_RESULT '+json.dumps(row,sort_keys=True),flush=True)
    return row
def main():
    dump('environment.json',{'python':sys.version,'platform':platform.platform(),'os_release':Path('/etc/os-release').read_text(),
         'pyscf':pyscf.__version__,'ase':ase.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'h5py':h5py.__version__,
         'threads':threadpool_info(),'lib_threads':lib.num_threads(),'affinity':sorted(os.sched_getaffinity(0)),
         'RLIMIT_AS':resource.getrlimit(resource.RLIMIT_AS),'RLIMIT_CPU':resource.getrlimit(resource.RLIMIT_CPU),
         'runner':{k:os.environ.get(k) for k in ('GITHUB_RUN_ID','GITHUB_RUN_ATTEMPT','GITHUB_SHA','ImageOS','ImageVersion')},
         'source_hashes':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ('c2/run_certification.py','c2/protocol.json','c2/STATIC_AUDIT.txt','requirements.txt','reference.json','.github/workflows/h2-once.yml')}})
    (OUT/'requirements.freeze.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True))
    rows=[]
    for spec in P['tests']:rows.append(test(*spec))
    ok=all(r['numerical_pass'] for r in rows) and all(p['num_threads']==1 for p in threadpool_info())
    dump('qa.json',{'tests_completed':len(rows),'NUMERICAL_SOLVER_QUALIFIED':'PASS' if ok else 'FAIL',
         'H2_PHYSICAL_BASELINE':'QUALIFIED' if ok else 'NOT_QUALIFIED','baseline_scope':'ONLY three finite-basis fixed-geometry BO points, not microhartree BO, Re, De, a full curve or storage',
         'ORIGINAL_DISCOVERY_CAPABILITY':'NOT_TESTED','microhartree_BO_at_B_and_C':all(abs(r['BO_difference_Eh'])<=1e-6 for r in rows[1:]),
         'E01_C2':'GO' if ok else 'NO-GO','historical_E01_C':'NO-GO unchanged'})
    return 0 if ok else 2
if __name__=='__main__':
    code=1
    try:
        with threadpool_limits(limits=1):code=main()
    except Exception as exc:
        dump('execution_error.json',{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()})
        traceback.print_exc()
    finally:
        signal.alarm(0)
        u=resource.getrusage(resource.RUSAGE_SELF)
        dump('resources.json',{'wall_seconds':time.monotonic()-START,'user_CPU_seconds':u.ru_utime,'system_CPU_seconds':u.ru_stime,
             'max_RSS_KiB':u.ru_maxrss,'threads_final':threadpool_info(),'completed_results':len(list(OUT.glob('*_result.json')))})
        for path in sorted(OUT.iterdir()):
            if not path.is_file():continue
            raw=path.read_bytes();data=base64.b64encode(zlib.compress(raw,9)).decode();parts=[data[i:i+8000] for i in range(0,len(data),8000)]
            for i,part in enumerate(parts):print('C2_FILE '+json.dumps({'name':path.name,'sha256':hashlib.sha256(raw).hexdigest(),'part':i,'parts':len(parts),'data':part}),flush=True)
    sys.exit(code)
