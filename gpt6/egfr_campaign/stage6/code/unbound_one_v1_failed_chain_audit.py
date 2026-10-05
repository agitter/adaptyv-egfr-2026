"""Matched local unbound minimization, not a folding/stability prediction."""
from pathlib import Path
import sys,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
import free_binder
from audit_refined import audit
cid=sys.argv[1];free_binder.S3=S;free_binder.main(cid,0)
d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());p=S/'intermediate/unbound'/f'{cid}_free_relaxed.pdb';a=audit(p,d);(S/'intermediate/unbound'/f'{cid}_free_independent.json').write_text(json.dumps(a,indent=2));print('UNBOUND STRICT',cid,a['pass'],flush=True)
