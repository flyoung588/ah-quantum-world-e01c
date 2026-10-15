"""Real two-electron Coulomb FCI. Execute only on the authorized hosted runner."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import base64,csv,hashlib,json,platform,resource,signal,subprocess,sys,time,traceback,zlib
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import pyscf,ase,h5py
from pyscf import gto,ao2mo,fci,scf,lib
from ase import Atoms
from ase.units import Bohr,Hartree
from threadpoolctl import threadpool_info,threadpool_limits

OUT=Path('results');OUT.mkdir(exist_ok=True)
START=time.monotonic();POINTS=[];CACHE={}
REF=json.loads(Path('reference.json').read_text())
lib.num_threads(1)
if hasattr(os,'sched_setaffinity'):
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})

def dump(name,value):
    (OUT/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def deadline(*args):raise TimeoutError('660 second computation ceiling; no retry permitted')

def molecule(R,basis):
    atoms=Atoms('H2',positions=[[0,0,-R*Bohr/2],[0,0,R*Bohr/2]])
    return gto.M(atom=[(s,tuple(p/Bohr)) for s,p in zip(atoms.get_chemical_symbols(),atoms.positions)],
        basis=basis,unit='Bohr',charge=0,spin=0,verbose=0,max_memory=1200)

def energy(R,basis,*,tol=1e-12,representation='overlap',use_cache=True):
    key=(basis,round(float(R),12),tol,representation)
    if use_cache and key in CACHE:return CACHE[key]
    t=time.monotonic();mol=molecule(float(R),basis)
    S=mol.intor_symmetric('int1e_ovlp');s,U=eigh(S)
    if s[0]<=1e-10:raise RuntimeError('Linear dependence cutoff reached, do not silently change basis')
    C=U/np.sqrt(s)[None,:]
    hcore=mol.intor_symmetric('int1e_kin')+mol.intor_symmetric('int1e_nuc')
    rhf_E=None
    if representation=='RHF_MO':
        mf=scf.RHF(mol);mf.conv_tol=1e-13;mf.max_cycle=100;mf.kernel()
        if not mf.converged:raise RuntimeError('Cross-check RHF not converged')
        C=mf.mo_coeff;rhf_E=float(mf.e_tot)
    n=C.shape[1];h1=C.T@hcore@C
    eri=ao2mo.kernel(mol,C,compact=False).reshape((n,)*4)
    solver=fci.direct_spin0.FCI(mol) if representation=='overlap' else fci.direct_spin1.FCI(mol)
    solver.conv_tol=tol;solver.max_cycle=200;solver.max_space=32;solver.max_memory=1200;solver.verbose=0
    E,ci=solver.kernel(h1,eri,n,(1,1),ecore=mol.energy_nuc())
    heff=solver.absorb_h1e(h1,eri,n,(1,1),.5)
    residual=solver.contract_2e(heff,ci,n,(1,1))-(float(E)-mol.energy_nuc())*ci
    ss=float(fci.spin_op.spin_square(ci,n,(1,1))[0])
    row={'basis':basis,'R_bohr':float(R),'R_Angstrom':float(R*Bohr),'norb':n,
        'energy_total_Eh':float(E),'nuclear_repulsion_Eh':float(mol.energy_nuc()),
        'FCI_converged':bool(solver.converged),'FCI_tol':tol,'S2':ss,
        'CI_norm':float(np.linalg.norm(ci)),'residual_norm_Eh':float(np.linalg.norm(residual)),
        'overlap_min_eigenvalue':float(s[0]),'orthogonality_error':float(np.max(np.abs(C.T@S@C-np.eye(n)))),
        'representation':representation,'RHF_energy_Eh':rhf_E,'wall_seconds':time.monotonic()-t}
    POINTS.append(row)
    with (OUT/'all_evaluations.jsonl').open('a') as f:f.write(json.dumps(row,sort_keys=True)+'\n')
    print('E01_POINT '+json.dumps(row,sort_keys=True),flush=True)
    CACHE[key]=row
    return row

def H_atom(basis):
    m=gto.M(atom='H 0 0 0',basis=basis,unit='Bohr',spin=1,verbose=0)
    h=m.intor_symmetric('int1e_kin')+m.intor_symmetric('int1e_nuc')
    return float(eigh(h,m.intor_symmetric('int1e_ovlp'),eigvals_only=True)[0])

def main():
    environment={'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),
        'pyscf':pyscf.__version__,'ase':ase.__version__,'numpy':np.__version__,'scipy':scipy.__version__,
        'h5py':h5py.__version__,'ase_Bohr_Angstrom':Bohr,'ase_Hartree_eV':Hartree,
        'lib_num_threads':lib.num_threads(),'threadpools':threadpool_info(),
        'cpu_affinity':sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else None,
        'os_release':Path('/etc/os-release').read_text(),'kernel':platform.release(),
        'runner':{k:os.environ.get(k) for k in ('GITHUB_RUN_ID','GITHUB_RUN_ATTEMPT','GITHUB_SHA','RUNNER_OS','RUNNER_ARCH','ImageOS','ImageVersion')},
        'source_hashes':{n:hashlib.sha256(Path(n).read_bytes()).hexdigest() for n in ('run_h2.py','requirements.txt','protocol.json','reference.json','.github/workflows/h2-once.yml')}}
    dump('environment.json',environment)
    (OUT/'requirements.freeze.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True))
    bases=['cc-pVDZ','cc-pVTZ','cc-pVQZ'];curve=[];summary=[];checks=[]
    distances=sorted(map(float,REF['points_Eh']))
    for basis in bases:
        for R in distances:
            row=energy(R,basis);reference=REF['points_Eh'][str(R)]
            curve.append({**row,'reference_Eh':reference,'error_Eh':row['energy_total_Eh']-reference})
        opt=minimize_scalar(lambda r:energy(r,basis)['energy_total_Eh'],bounds=(1.3,1.5),method='bounded',options={'xatol':2e-7,'maxiter':24})
        r=float(opt.x);emin=energy(r,basis)['energy_total_Eh'];atomic=H_atom(basis)
        locations=[]
        for step in (.002,.001):
            em=energy(r-step,basis)['energy_total_Eh'];ep=energy(r+step,basis)['energy_total_Eh']
            curvature=(em+ep-2*emin)/(step*step)
            vertex=r+(em-ep)*step/(2*(em+ep-2*emin))
            locations.append({'step_bohr':step,'vertex_bohr':vertex,'curvature_Eh_per_bohr2':curvature})
        anchor=energy(1.4011,basis)
        loose=energy(1.4011,basis,tol=1e-10,use_cache=False)
        cross=energy(1.4011,basis,representation='RHF_MO',use_cache=False)
        check={'basis':basis,'same_basis_representation_delta_Eh':cross['energy_total_Eh']-anchor['energy_total_Eh'],
            'tol_1e10_vs_1e12_delta_Eh':loose['energy_total_Eh']-anchor['energy_total_Eh'],
            'FCI_below_RHF':cross['energy_total_Eh']<=cross['RHF_energy_Eh']+1e-9}
        checks.append(check)
        summary.append({'basis':basis,'Re_bohr':r,'Re_Angstrom':r*Bohr,'Emin_Eh':emin,
            'optimization_success':bool(opt.success),'optimization_nfev':int(opt.nfev),
            'H_atom_same_basis_Eh':atomic,'H_atom_basis_error_Eh':atomic+.5,
            'De_basis_consistent_Eh':2*atomic-emin,'De_basis_consistent_eV':(2*atomic-emin)*Hartree,
            'depth_exact_endpoint_diagnostic_Eh':-1-emin,
            'depth_vs_reference_point_Eh':2*atomic-emin-REF['reference_depth_at_1_4011_Eh'],
            'Re_minus_rounded_reference_bohr':r-1.4011,
            'anchor_error_Eh':anchor['energy_total_Eh']-REF['points_Eh']['1.4011'],
            'local_parabolic_checks':locations,
            'Re_step_sensitivity_bohr':abs(locations[0]['vertex_bohr']-locations[1]['vertex_bohr'])})
        dump('basis_summaries.json',summary);dump('numerical_checks.json',checks)
    with (OUT/'curve.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(curve[0]));writer.writeheader();writer.writerows(curve)
    dump('curve.json',curve)
    convergence=[]
    for R in distances:
        rr=[next(p for p in curve if p['basis']==b and p['R_bohr']==R) for b in bases]
        convergence.append({'R_bohr':R,'TZ_minus_DZ_Eh':rr[1]['energy_total_Eh']-rr[0]['energy_total_Eh'],
            'QZ_minus_TZ_Eh':rr[2]['energy_total_Eh']-rr[1]['energy_total_Eh'],
            'QZ_error_Eh':rr[2]['error_Eh'],'absolute_error_decreases_DZ_TZ_QZ':abs(rr[2]['error_Eh'])<abs(rr[1]['error_Eh'])<abs(rr[0]['error_Eh'])})
    numerical=all(p['FCI_converged'] and abs(p['S2'])<=1e-7 and p['residual_norm_Eh']<=1e-7 for p in POINTS)
    representation=all(abs(c['same_basis_representation_delta_Eh'])<=1e-8 and abs(c['tol_1e10_vs_1e12_delta_Eh'])<=1e-8 and c['FCI_below_RHF'] for c in checks)
    valid_variational=all(p['error_Eh']>=-1e-8 for p in curve)
    optimization=all(s['optimization_success'] and 1.3<s['Re_bohr']<1.5 for s in summary)
    pools_ok=all(t['num_threads']==1 for t in threadpool_info()) and lib.num_threads()==1
    qz=summary[-1];high_precision=abs(qz['anchor_error_Eh'])<=1e-6
    qa={'computational_completion':'PASS','all_FCI_converged_spin_and_residual':numerical,
        'alternative_representation_and_tolerance':representation,'variational_lower_bound_sanity':valid_variational,
        'all_optimized_in_bracket':optimization,'single_thread_verified':pools_ok,
        'reference_point_count':len(distances),'curve_row_count':len(curve),'total_evaluations':len(POINTS),
        'near_equilibrium_1mEh_diagnostic':abs(qz['anchor_error_Eh'])<=.001,
        'near_equilibrium_1microEh_diagnostic':high_precision,
        'literature_precision_reproduced':False,
        'GO_for_bounded_finite_basis_reproduction':numerical and representation and valid_variational and optimization and pools_ok,
        'GO_for_microhartree_BO_accuracy':high_precision and numerical,
        'claim_ceiling':'Known H2 finite-basis reproduction only, no original discovery, no sample/storage properties'}
    dump('convergence.json',convergence);dump('qa.json',qa)
    dump('error_budget.json',{'method':'All configurations within each finite orbital basis for two electrons; no correlation truncation within that space',
        'basis':'Primary remaining error relative to independent BO reference; DZ/TZ/QZ differences do not rigorously bound CBS error',
        'numerical':'CI residual, spin, tolerance change, alternative orbital representation/solver, and local minimum step-sensitivity explicitly measured',
        'physical_model':'Infinite nuclear mass nonrelativistic BO; excludes nonadiabatic/adiabatic, relativity, QED, ZPE; De is not experimental D0',
        'reference':'Published near-minimum r=1.4011 is rounded, not exact independently optimized Re; depth comparison carries this distinction',
        'BSSE':'Same-basis separate-atom endpoint; no counterpoise correction applied, nor assumed zero; exact-atomic-endpoint diagnostic also reported',
        'independence':'Literature BO numbers independent of Gaussian FCI; two FCI implementations are numerical cross-checks, not independent physical validation'})
    return 0 if qa['GO_for_bounded_finite_basis_reproduction'] else 2

def export():
    dump('resources.json',{'wall_seconds':time.monotonic()-START,'process_user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,
        'process_system_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'max_rss_KiB_linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'threadpools_final':threadpool_info(),'points_completed':len(POINTS)})
    for p in sorted(OUT.iterdir()):
        if not p.is_file():continue
        raw=p.read_bytes();data=base64.b64encode(zlib.compress(raw,9)).decode();parts=[data[i:i+8000] for i in range(0,len(data),8000)]
        for i,part in enumerate(parts):
            print('E01_FILE '+json.dumps({'name':p.name,'sha256':hashlib.sha256(raw).hexdigest(),'part':i,'parts':len(parts),'encoding':'zlib+base64','data':part}),flush=True)

if __name__=='__main__':
    signal.signal(signal.SIGALRM,deadline);signal.alarm(660)
    code=1
    try:
        with threadpool_limits(limits=1):code=main()
    except Exception as exc:
        dump('execution_error.json',{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()})
        traceback.print_exc()
    finally:
        signal.alarm(0);export()
    sys.exit(code)
