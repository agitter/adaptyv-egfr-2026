"""Equal additional refinement for the unmutated parent; no new sequence."""
from pathlib import Path
import sys,json,copy,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';C=S/'intermediate/matched_controls'
sys.path[:0]=[str(S/'code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import refine_new
from mm_refine import app,u
from audit_refined import audit

def run(cid):
 source=R/'stage3/intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');dp=R/'stage3/intermediate/designs'/(cid+'.json');d=json.loads(dp.read_text());p=app.PDBFile(str(source));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));rows=[{a.name:xyz[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'} for r in p.topology.residues() if r.chain.id=='B'];assert len(rows)==len(d['structure'])
 for r,atoms in zip(d['structure'],rows):r['atoms']=atoms
 d['stage4_control']={'role':'unmutated parent after equal additional joint refinement, not a new sequence candidate','source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'mutation_count':0};d['metrics']['scores_valid']=False
 (C/'intermediate/designs').mkdir(parents=True,exist_ok=True);(C/'intermediate/designs'/(cid+'.json')).write_text(json.dumps(d,indent=2));refine_new.S3=C;print(refine_new.run(cid,'human6ARU',450),flush=True);ref=C/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');a=audit(ref,d);(C/'intermediate/refined'/(cid+'_independent_v2.json')).write_text(json.dumps(a,indent=2));print('AUDIT',a['pass'],flush=True)
if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'C00003')
