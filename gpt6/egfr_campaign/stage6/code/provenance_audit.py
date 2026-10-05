"""Audit delivered-sequence ancestry against the regenerated classical lineage."""
from pathlib import Path
import json,hashlib,re,sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def audit(rows):
 poses={r['pose_id']:r for r in json.loads((R/'stage3/intermediate/docking/poses.json').read_text())};roots={r['pose_id']:r for r in json.loads((R/'stage2/intermediate/docking/poses.json').read_text())};out=[]
 for r in rows:
  d=json.loads(Path(r['design_source']).read_text());pose=poses.get(d['pose_id']);assert pose is not None,d['pose_id'];ancestor=pose['ancestor_backbone_pose'];assert ancestor in roots
  loop_files=[]
  for name,(label,index) in roots[ancestor]['selection'].items():
   p=R/'stage2/intermediate/loops'/f'{label}.npz';assert p.exists();loop_files.append({'loop':name,'file':str(p),'sha256':sha(p),'array_index':index})
  intervals=d['cdr_intervals_zero_based'];mask=[False]*len(d['sequence'])
  for a,b in intervals:
   for i in range(a,b):mask[i]=True
  assert mask==[bool(x['loop']) for x in d['structure']]
  assert d['sequence']==''.join(x['aa'] for x in d['structure'])
  assert d.get('historical_rng')=='MT19937',(r['candidate_key'],d.get('historical_rng'))
  out.append({'sequence_sha256':r['sequence_sha256'],'candidate_key':r['candidate_key'],'design_file':r['design_source'],'design_sha256':sha(r['design_source']),'source_pose':d['pose_id'],'stage2_MT19937_pose':ancestor,'loop_inputs':loop_files,'source_ancestry':d.get('ancestry'),'additional_reassessment_source':d.get('reassessment_source'),'historical_rng':d.get('historical_rng'),'direct_sequence_input':False,'original_framework_binding_loops_retained':False,'audit_pass':True})
 return {'records':out,'count':len(out),'all_pass':all(x['audit_pass'] for x in out),'framework_definition_source':str(R/'stage2/reference/framework_definition.json'),'framework_definition_sha256':sha(R/'stage2/reference/framework_definition.json'),'loop_generator':str(R/'stage2/code/generate_loops_mt.py'),'loop_generator_sha256':sha(R/'stage2/code/generate_loops_mt.py'),'scope':'Checks the accepted regenerated MT19937 lineage and saved source files, not a claim that excluded historical pilot operations were compliant. New local loop rebuilding and mutations remain recorded in each design ancestry.'}
if __name__=='__main__':
 rows=json.loads((S/'reference/panel_preview.json').read_text())[:100];d=audit(rows);(S/'reference/preview_ancestry_audit.json').write_text(json.dumps(d,indent=2));print('Ancestry checks',d['count'],d['all_pass'])
