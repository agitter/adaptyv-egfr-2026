"""Reuse historical all-atom/proton-state pipeline in isolated stage5 outputs."""
from pathlib import Path
import sys,json
R=Path('/mnt/data/egfr_campaign');S5=R/'stage5'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import refine_new
from audit_refined import audit
from mm_refine import app,u
import numpy as np
from scipy.spatial import cKDTree

def run(cid,species='human6ARU'):
 refine_new.S3=S5;r=refine_new.run(cid,species,450);d=json.loads((S5/'intermediate/designs'/f'{cid}.json').read_text());prefix=S5/'intermediate/refined'/f'{cid}_{species}';p=Path(str(prefix)+'_relaxed.pdb')
 if p.exists():
  a=audit(p,d);a.update(species=species);Path(str(prefix)+'_independent.json').write_text(json.dumps(a,indent=2));print('INDEPENDENT',cid,species,a['pass'],len(a['peptide_omega_warnings']),flush=True)
 print(r,flush=True)
if __name__=='__main__':run(*sys.argv[1:])
