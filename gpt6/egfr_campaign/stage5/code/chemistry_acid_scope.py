"""Explicitly approved full engineered-acid scope for E104 prototypes."""
from pathlib import Path
import sys,json,shutil,hashlib
S=Path('/mnt/data/egfr_campaign/stage5');R=S.parent
sys.path[:0]=[str(S/'code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import acid_histidine_six_sites as a
from independent_mixed_priors import run as independent_run
if __name__=='__main__':
 cid=sys.argv[1];pres=json.loads((S/'reference/six_site_scope_predeclaration.json').read_text());positions=pres['designs'][cid];root=S/'intermediate/chemistry_acid_scope';(root/'intermediate/refined').mkdir(parents=True,exist_ok=True);source=S/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';dest=root/'intermediate/refined'/source.name;shutil.copy2(source,dest);a.S3=root;a.probe_acids(cid,tuple(positions));path=root/'intermediate/acid_ensemble'/f'{cid}_acid_histidine.json';d=json.loads(path.read_text());assert len(d['sites'])==len(positions)+2 and len(d['microstates'])==3**len(d['sites']);independent_run(path,cid+'_engineered_acids')
