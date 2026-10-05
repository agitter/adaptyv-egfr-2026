"""Matched unbound minimization using the preserved classical engine."""
from pathlib import Path
import sys,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
import free_binder
from audit_refined import audit
from mm_refine import app,mm
mm.Platform.getPlatformByName('CPU').setPropertyDefaultValue('Threads','1')
cid=sys.argv[1]
if (S/'intermediate/unbound'/f'{cid}_free_audit.json').exists():raise FileExistsError('Do not silently overwrite completed unbound work')
free_binder.S3=S;free_binder.main(cid,0)
d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());p=S/'intermediate/unbound'/f'{cid}_free_relaxed.pdb'
file=app.PDBFile(str(p));chains=list(file.topology.chains());assert len(chains)==1 and len(list(chains[0].residues()))==len(d['sequence'])
a=audit(p,d,chain=chains[0].id);(S/'intermediate/unbound'/f'{cid}_free_independent.json').write_text(json.dumps(a,indent=2));print('UNBOUND STRICT',cid,a['pass'],flush=True)
