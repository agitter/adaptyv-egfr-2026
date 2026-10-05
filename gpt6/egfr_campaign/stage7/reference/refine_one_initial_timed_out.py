"""Run the preserved historical MM/protonation pipeline on stage7 only."""
from pathlib import Path
import sys,json,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/x/'code') for x in ['stage6','stage5','stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
import refine_new
from audit_refined import audit
refine_new.S3=S
if __name__=='__main__':
 cid=sys.argv[1];species=sys.argv[2] if len(sys.argv)>2 else 'human6ARU';t=time.time()
 r=refine_new.run(cid,species,450);p=S/'intermediate/refined'/f'{cid}_{species}_relaxed.pdb'
 if p.exists():
  d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());a=audit(p,d);p.with_name(p.name.replace('_relaxed.pdb','_independent.json')).write_text(json.dumps(a,indent=2));print('STRICT',a['pass'],'omega_warnings',len(a['peptide_omega_warnings']),flush=True)
 print('RESULT',json.dumps(r),'seconds',time.time()-t,flush=True)
