from pathlib import Path
import json,collections,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';designs={};seqs={};out=[]
for stage in ['stage2','stage3','stage4','stage5']:
 root=R/stage
 for p in sorted((root/'intermediate/designs').glob('*.json')):
  d=json.loads(p.read_text())
  if not isinstance(d,dict) or not d.get('sequence'):continue
  key=stage+':'+d['candidate_id'];designs[key]=p
  seqs.setdefault(d['sequence'],[]).append(key)
 for p in sorted((root/'intermediate/refined').glob('*_proton.json')):
  d=json.loads(p.read_text());stem=p.stem
  if '_human6ARU' in stem:cid=stem.split('_human6ARU')[0];species='human'
  elif '_mouseAF' in stem:cid=stem.split('_mouseAF')[0];species='mouse'
  else:continue
  key=stage+':'+cid;dp=designs.get(key)
  if dp is None:continue
  design=json.loads(dp.read_text());pro=d['proton_model'];r={'key':key,'stage':stage,'id':cid,'species':species,'model':str(p),'coordinate_source':d['source'],'n_sites':len(d['sites']),'minimum_contrast':pro['minimum_contrast_kcal'],'maximum_contrast':pro['maximum_contrast_kcal'],'sequence_sha256':hashlib.sha256(design['sequence'].encode()).hexdigest(),'branch':design.get('branch','initial'),'acid_refined':'_acid_' in stem,'sites':d['sites'],'sequence':design['sequence'],'design_source':str(dp)}
  rep=[x for x in pro['scenarios'] if x['solute_dielectric']==4 and x['kappa_nm_inverse']==1.25 and x['assumed_unbound_pKa']==6.3 and x['assumed_unbound_HIE_fraction']==.5]
  if rep:r.update(acid_proxy=rep[0]['low_pH']['conditional_binding_energy_kcal'],neutral_proxy=rep[0]['high_pH']['conditional_binding_energy_kcal'])
  out.append(r)
(S/'reference/model_inventory.json').write_text(json.dumps(out,indent=2));(S/'reference/sequence_lineages.json').write_text(json.dumps(seqs,indent=2))
print('Unique sequences:',len(seqs),'design records:',len(designs),'core refined models:',len(out))
for x in sorted([x for x in out if x['species']=='human' and not x['acid_refined']],key=lambda x:-x['minimum_contrast']):
 print(x['key'],x['n_sites'],round(x['minimum_contrast'],3),round(x.get('acid_proxy',0),1),round(x.get('neutral_proxy',0),1),x['branch'])
