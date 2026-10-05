"""Repair a chain-label mismatch in postprocessing, without rerunning minimization.
OpenMM's saved single-chain free structures here have chain A. The first wrapper
asked the two-chain-complex checker for chain B and raised an assertion. Preserve
those failed wrapper logs and independently audit the actual exported chain.
"""
from pathlib import Path
import sys,json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from audit_refined import audit
from mm_refine import app
results=[]
for p in sorted((S/'intermediate/unbound').glob('*_free_relaxed.pdb')):
 cid=p.name.split('_free_relaxed')[0];dp=S/'intermediate/designs'/f'{cid}.json';ap=p.with_name(cid+'_free_audit.json')
 if not ap.exists():continue
 d=json.loads(dp.read_text());file=app.PDBFile(str(p));chains=list(file.topology.chains());assert len(chains)==1 and len(list(chains[0].residues()))==len(d['sequence'])
 result=audit(p,d,chain=chains[0].id);result.update(exported_chain_id=chains[0].id,source=str(p),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),correction='Inspect the actual single exported chain rather than assuming the complex chain B. No energy or coordinate change.')
 p.with_name(cid+'_free_independent.json').write_text(json.dumps(result,indent=2));results.append({'candidate_id':cid,'pass':result['pass'],'exported_chain':chains[0].id,'source':str(p)});print(cid,result['pass'],chains[0].id,flush=True)
(S/'reference/unbound_chain_audit_recovery.json').write_text(json.dumps({'recovered_checks':len(results),'results':results,'qualification':'Original nonzero wrapper exits preserved; this corrects postprocessing, not an observed molecular defect.'},indent=2))
