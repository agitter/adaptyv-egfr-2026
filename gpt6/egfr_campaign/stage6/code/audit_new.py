"""Apply exactly the archived-model audit to newly completed models."""
from pathlib import Path
import json,hashlib,time
from global_audit import model_audit,structural_audit,sequence_audit,sha
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';out=S/'intermediate/uniform_new';out.mkdir(exist_ok=True)
old=json.loads((S/'reference/campaign_audit.json').read_text());bycdr={}
for x in old:
 nv=x.get('local_novelty')
 if nv and x['CDRs'][2]==nv.get('cdr3_query'):
  c=x['CDRs'][2]
  if c not in bycdr or nv['max_edit_identity']>bycdr[c]['max_edit_identity']:bycdr[c]=nv
completed=set()
for batch,species in [('human_factorial','human'),('mouse_factorial','mouse')]:
 p=S/'reference'/f'{batch}_execution.json'
 if p.exists():
  for r in json.loads(p.read_text())['records']:
   if r['returncode']==0:completed.add((r['id'],species))
models=[];rows=[]
p=S/'reference/diversity_followup_execution.json'
if p.exists():
 for r in json.loads(p.read_text())['records']:
  if r['returncode']==0:
   cid,suffix=r['id'].split('_',1);completed.add((cid,'human' if suffix=='human6ARU' else 'mouse'))
p=S/'reference/mouse_completion_execution.json'
if p.exists():
 for r in json.loads(p.read_text())['records']:
  if r['returncode']==0:
   cid,suffix=r['id'].split('_',1);completed.add((cid,'mouse'))
for dp in sorted((S/'intermediate/designs').glob('[JKL]*.json')):
 design=json.loads(dp.read_text());seq=design['sequence'];nv=bycdr.get(design['cdr_sequences'][2]);g=structural_audit(design);q=sequence_audit(design)
 rows.append({'key':'stage6:'+design['candidate_id'],'sequence':seq,'design_source':str(dp),'design_sha256':sha(dp),'length':len(seq),'branch':design['branch'],'pose_id':design.get('pose_id'),'CDRs':design['cdr_sequences'],'intervals':design['cdr_intervals_zero_based'],'sequence_audit':q,'initial_geometry':g,'initial_context_pass':bool(g['geometry_prefilter_pass'] and all(g[c]['heavy_atoms_below_2A']==0 for c in ['1IVO','1NQL']) and not q['liabilities']['glycosylation_sequons']),'local_novelty_pass':bool(nv is not None and nv['max_edit_identity']<.7),'local_novelty':nv,'novelty_method':'Exact CDR3 sequence equality to an already screened campaign CDR3; source retained. Not official IMGT numbering.'})
 for species,suffix in [('human','human6ARU'),('mouse','mouseAF')]:
  if (design['candidate_id'],species) not in completed:continue
  p=S/'intermediate/refined'/f'{design["candidate_id"]}_{suffix}_proton.json'
  if not p.exists():continue
  d=json.loads(p.read_text());pro=d['proton_model'];rec={'key':'stage6:'+design['candidate_id'],'stage':'stage6','id':design['candidate_id'],'species':species,'model':str(p),'coordinate_source':d['source'],'n_sites':len(d['sites']),'minimum_contrast':pro['minimum_contrast_kcal'],'maximum_contrast':pro['maximum_contrast_kcal'],'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'branch':design['branch'],'acid_refined':False,'sites':d['sites'],'sequence':seq,'design_source':str(dp)}
  rep=[x for x in pro['scenarios'] if x['solute_dielectric']==4 and x['kappa_nm_inverse']==1.25 and x['assumed_unbound_pKa']==6.3 and x['assumed_unbound_HIE_fraction']==.5]
  if rep:rec.update(acid_proxy=rep[0]['low_pH']['conditional_binding_energy_kcal'],neutral_proxy=rep[0]['high_pH']['conditional_binding_energy_kcal'])
  dest=out/(p.stem+'.json')
  if dest.exists():res=json.loads(dest.read_text())
  else:res=model_audit(rec);dest.write_text(json.dumps(res,indent=2))
  res['result_source']=str(dest);models.append({k:v for k,v in res.items() if k not in ['strict_geometry','reference_geometry','sequence']})
  print(rec['key'],species,res['combined_geometry_pass'],round(res['independent_priors']['minimum_contrast_kcal'],5),flush=True)
(S/'reference/new_model_summary.json').write_text(json.dumps(models,indent=2));(S/'reference/new_design_audit.json').write_text(json.dumps(rows,indent=2));print('Audited',len(models),'models',len(rows),'sequences')
