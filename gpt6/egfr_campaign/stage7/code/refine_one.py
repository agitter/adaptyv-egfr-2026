"""Run the preserved historical MM/protonation pipeline on stage7 only."""
from pathlib import Path
import sys,json,time,faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(80,repeat=True)
print("PROCESS_START",time.time(),flush=True)
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/x/'code') for x in ['stage6','stage5','stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
import refine_new
from audit_refined import audit
refine_new.S3=S
from mm_refine import mm,app
mm.Platform.getPlatformByName("CPU").setPropertyDefaultValue("Threads","1")
_orig_h=app.Modeller.addHydrogens
def _timed_h(self,*a,**kw):
 t=time.time();print("ADD_H_START",self.topology.getNumAtoms(),flush=True)
 r=_orig_h(self,*a,**kw);print("ADD_H_DONE",time.time()-t,flush=True);return r
app.Modeller.addHydrogens=_timed_h
print("IMPORTS_DONE",time.time(),flush=True)
if __name__=='__main__':
 cid=sys.argv[1];species=sys.argv[2] if len(sys.argv)>2 else 'human6ARU';t=time.time()
 r=refine_new.run(cid,species,450);p=S/'intermediate/refined'/f'{cid}_{species}_relaxed.pdb'
 if p.exists():
  d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());a=audit(p,d);p.with_name(p.name.replace('_relaxed.pdb','_independent.json')).write_text(json.dumps(a,indent=2));print('STRICT',a['pass'],'omega_warnings',len(a['peptide_omega_warnings']),flush=True)
 print('RESULT',json.dumps(r),'seconds',time.time()-t,flush=True)
