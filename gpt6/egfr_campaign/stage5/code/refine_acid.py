"""Matched all-interface-HIP structural hypothesis, independently audited.
No contact-proposal restraint is retained. Compare B00000 and chemical variants
within this same600-iteration condition; do not cherry-pick the better geometry.
"""
from pathlib import Path
import sys,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import refine_state
from audit_refined import audit
if __name__=='__main__':
 cid=sys.argv[1];species=sys.argv[2] if len(sys.argv)>2 else 'human6ARU';refine_state.S3=S;r=refine_state.run(cid,species,600);p=S/'intermediate/refined'/f'{cid}_{species}_acid_relaxed.pdb';d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text())
 if p.exists():
  a=audit(p,d);a.update(species=species,structural_hypothesis='all-interface histidines protonated');(p.parent/(p.name.replace('_relaxed.pdb','_independent.json'))).write_text(json.dumps(a,indent=2));print('INDEPENDENT',cid,species,a['pass'],flush=True)
 print(r,flush=True)
