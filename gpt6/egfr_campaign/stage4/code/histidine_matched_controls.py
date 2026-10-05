"""Equal-refinement controls for F107H probes; no new sequences.
This closes a second-level refinement confound: probes started from stage4
relaxed parents, so they need an additional unmutated round as comparator.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(S/'code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import refine_new
from mm_refine import app,u
from audit_refined import audit

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(probe):
 dprobe=json.loads((S/'intermediate/designs'/(probe+'.json')).read_text());an=dprobe['ancestry'];source=Path(an['relaxed_backbone_source']);dp=Path(an['parent_file']);assert sha(source)==an['relaxed_backbone_sha256'] and sha(dp)==an['parent_sha256'];d=json.loads(dp.read_text());cid=d['candidate_id'];out=S/'intermediate/histidine_matched_controls'/probe
 p=app.PDBFile(str(source));x=np.asarray(p.positions.value_in_unit(u.angstrom));rows=[{a.name:x[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'} for r in p.topology.residues() if r.chain.id=='B'];assert len(rows)==len(d['structure'])
 for r,a in zip(d['structure'],rows):r['atoms']=a
 assert d['sequence'][106]=='F' and dprobe['sequence'][106]=='H';assert sum(a!=b for a,b in zip(d['sequence'],dprobe['sequence']))==1
 d['stage4_control']={'role':'additional equal-refinement parent for F107H probe, not a new sequence','probe':probe,'source':str(source),'source_sha256':sha(source),'mutation_count':0};d['metrics']['scores_valid']=False
 dest=out/'intermediate/designs';dest.mkdir(parents=True,exist_ok=True);(dest/(cid+'.json')).write_text(json.dumps(d,indent=2));refine_new.S3=out;print(refine_new.run(cid,'human6ARU',450),flush=True);ref=out/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');a=audit(ref,d);(ref.parent/(cid+'_independent_v2.json')).write_text(json.dumps(a,indent=2));print('AUDIT',a['pass'],flush=True)
if __name__=='__main__':run(sys.argv[1])
