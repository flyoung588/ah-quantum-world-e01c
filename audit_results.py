"""Offline arithmetic/provenance QA of saved data only; no quantum execution."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).parent;R=ROOT/'results'
def load(p):return json.loads(Path(p).read_text())
rows=[json.loads(s) for s in (R/'all_evaluations.jsonl').read_text().splitlines()]
protocol=load(ROOT/'protocol.json');ref=load(ROOT/'reference.json')
threshold=protocol['numerical_QA']['residual_Eh_max']
manifest=load(R/'SHA256_MANIFEST.json')
hashes=all(hashlib.sha256((R/r['path']).read_bytes()).hexdigest()==r['sha256'] for r in manifest['files'])
out={'hashes_match':hashes,'quantum_computations_executed_by_audit':0,'completed_evaluations':len(rows),
    'frozen_residual_threshold_Eh':threshold,'bases':{}}
for b in protocol['basis_sequence']:
    rr=[r for r in rows if r['basis']==b]
    grid={r['R_bohr']:r for r in rr if str(r['R_bohr']) in ref['points_Eh'] and r['representation']=='overlap' and r['FCI_tol']==1e-12}
    anchor=grid[1.4011]
    out['bases'][b]={'fixed_grid_points':len(grid),'failed_residuals':sum(r['residual_norm_Eh']>threshold for r in rr),
        'anchor_energy_Eh':anchor['energy_total_Eh'],'anchor_error_Eh':anchor['energy_total_Eh']-ref['points_Eh']['1.4011'],
        'anchor_S2':anchor['S2']}
for s in load(R/'basis_summaries.json'):
    assert abs(2*s['H_atom_same_basis_Eh']-s['Emin_Eh']-s['De_basis_consistent_Eh'])<1e-13
out['De_arithmetic_matches']=True
out['execution_error']=load(R/'execution_error.json')['type']
out['overall']='NO_GO_INCOMPLETE_AND_NUMERICAL_QA'
print(json.dumps(out,indent=2))
