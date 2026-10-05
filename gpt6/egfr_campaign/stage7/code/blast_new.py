"""Use classical BLAST with exact-reference controls on every new sequence."""
from pathlib import Path
import json,sys,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';sys.path.insert(0,str(R/'stage6/code'))
import novelty_superset as b
b.S=S;b.OUT=S/'intermediate/blast';b.BIN=S/'runtime/bin';rows=[]
for p in sorted((S/'intermediate/designs').glob('M*.json')):
 d=json.loads(p.read_text());rows.append({'query_id':'Q'+hashlib.sha256(d['sequence'].encode()).hexdigest()[:20],'sequence':d['sequence'],'sources':[str(p)],'candidate_id':d['candidate_id']})
assert len({r['sequence'] for r in rows})==len(rows)
(S/'reference/novelty_query_superset.json').write_text(json.dumps(rows,indent=2))
for name in ['swissprot','pdb','antibodies_augmented']:b.run(name)
(S/'reference/panel_blast_complete.json').write_text(json.dumps({'complete':True,'candidate_count':len(rows),'databases':['swissprot','pdb','antibodies_augmented']}));print('ALL_SEARCHES_COMPLETE',flush=True)
