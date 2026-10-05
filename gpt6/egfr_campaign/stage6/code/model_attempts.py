"""Include failed refinement attempts even when no proton-model file was written."""
from pathlib import Path
import json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';rows=[]
for stage in ['stage2','stage3','stage4','stage5','stage6']:
 for p in sorted((R/stage/'intermediate/refined').glob('*_audit.json')):
  if '_human6ARU' in p.name:cid=p.name.split('_human6ARU')[0];species='human'
  elif '_mouseAF' in p.name:cid=p.name.split('_mouseAF')[0];species='mouse'
  else:continue
  dp=R/stage/'intermediate/designs'/f'{cid}.json'
  if not dp.exists():continue
  d=json.loads(p.read_text());design=json.loads(dp.read_text())
  if 'pass' not in d:continue
  rows.append({'key':stage+':'+cid,'species':species,'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'sequence_sha256':hashlib.sha256(design['sequence'].encode()).hexdigest(),'minimizer_geometry_pass':bool(d['pass']),'full_context_pass':d.get('postrefinement_exclusion_pass'),'basic_refinement_pass':bool(d['pass'] and d.get('postrefinement_exclusion_pass',True)),'acid_refined':'_acid_' in p.name})
(S/'reference/refinement_attempts.json').write_text(json.dumps(rows,indent=2));print('Refinement attempts',len(rows),'failed',sum(not x['basic_refinement_pass'] for x in rows))
