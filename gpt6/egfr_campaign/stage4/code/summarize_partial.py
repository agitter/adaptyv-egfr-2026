"""Evidence inventory; does not turn diagnostics into biological predictions."""
from pathlib import Path
import json,sys
S=Path('/mnt/data/egfr_campaign/stage4');R=S.parent

def rep(model):
    return next((r for r in model['scenarios'] if r['solute_dielectric']==1. and r['kappa_nm_inverse']==1.25 and r.get('assumed_unbound_pKa')==6.3 and r.get('assumed_unbound_HIE_fraction')==.5),None)

def main():
 out={'human':{},'mouse':{},'dynamics_corrected':{},'pilot_dynamics':{},'design_records':len(list((S/'intermediate/designs').glob('*.json')))}
 for species,key in [('human6ARU','human'),('mouseAF','mouse')]:
  for p in sorted((S/'intermediate/refined').glob(f'*_{species}_proton.json')):
   cid=p.name.split('_'+species)[0];d=json.loads(p.read_text());m=d['proton_model'];r=rep(m);ind=S/'intermediate/independent_v2'/f'{cid}_{species}_relaxed.json';geom=S/'intermediate/refined'/f'{cid}_{species}_audit.json'
   out[key][cid]={'minimum_contrast_kcal':m['minimum_contrast_kcal'],'maximum_contrast_kcal':m['maximum_contrast_kcal'],'acid_on_all54':m['all_scenarios_acid_on'],'microstates':len(d['microstates']),'scenario_count':len(m['scenarios']),'representative_acid_proxy_kcal':r['low_pH']['conditional_binding_energy_kcal'] if r else None,'representative_neutral_proxy_kcal':r['high_pH']['conditional_binding_energy_kcal'] if r else None,'independent_geometry':json.loads(ind.read_text()) if ind.exists() else None,'refinement_geometry_pass':json.loads(geom.read_text())['pass'] if geom.exists() else None}
 for p in sorted((S/'intermediate/dynamics').glob('*.json')):
  if p.name.endswith('_progress.json'):continue
  d=json.loads(p.read_text())
  if 'summary' not in d:continue
  key='dynamics_corrected' if '_v2_' in p.name else 'pilot_dynamics';out[key][p.stem]={'candidate_id':d['candidate_id'],'seed':d['seed'],'warmup_ps':d['warmup_ps'],'production_ps':d['production_ps'],**d['summary'],'final_geometry_pass':d['final_geometry_audit']['pass']}
 novelty=S/'reference/novelty_summary.json'
 if novelty.exists():out['novelty']=json.loads(novelty.read_text())
 (S/'output/progress_inventory.json').write_text(json.dumps(out,indent=2))
 for key in ['human','mouse']:
  print('\n'+key)
  for cid,r in out[key].items():print(cid,round(r['minimum_contrast_kcal'],3),round(r['representative_acid_proxy_kcal'],2) if r['representative_acid_proxy_kcal'] is not None else None,round(r['representative_neutral_proxy_kcal'],2) if r['representative_neutral_proxy_kcal'] is not None else None,r['independent_geometry']['pass'] if r['independent_geometry'] else None)
 print('\nCorrected trajectories',out['dynamics_corrected'])
if __name__=='__main__':main()
