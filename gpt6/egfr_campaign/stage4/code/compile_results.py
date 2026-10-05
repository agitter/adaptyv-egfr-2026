"""Compile evidence without converting proxies into claimed biological success."""
from pathlib import Path
import json,hashlib,collections,datetime
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
def read(p):return json.loads(p.read_text())
def short_model(p):
 d=read(p);m=d['proton_model'];r=next((z for z in m['scenarios'] if z['solute_dielectric']==1 and z['kappa_nm_inverse']==1.25 and z.get('assumed_unbound_pKa')==6.3 and z.get('assumed_unbound_HIE_fraction')==.5),None)
 return {'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'minimum_contrast_kcal':m['minimum_contrast_kcal'],'maximum_contrast_kcal':m['maximum_contrast_kcal'],'scenario_count':len(m['scenarios']),'microstate_count':len(d['microstates']),'all_acid_on':m['all_scenarios_acid_on'],'representative_acid_proxy_kcal':r['low_pH']['conditional_binding_energy_kcal'] if r else None,'representative_neutral_proxy_kcal':r['high_pH']['conditional_binding_energy_kcal'] if r else None}
def main():
 designs={p.stem:read(p) for p in (S/'intermediate/designs').glob('*.json')};assert len(designs)==63
 assert len({d['sequence'] for d in designs.values()})==63
 old={};historical=collections.defaultdict(list);counts={}
 for stage in ['stage2','stage3','stage4']:
  count=0
  for p in (R/stage/'intermediate/designs').glob('*.json'):
   d=read(p)
   if 'sequence' not in d:continue
   historical[d['sequence']].append(stage+'/'+p.stem);count+=1
   if stage in ['stage3','stage4']:old[stage+'/'+p.stem]=d['sequence']
  counts[stage]=count
 inv=read(S/'output/progress_inventory.json');dynamics={}
 for p in (S/'intermediate/dynamics').glob('*_v2_s*.json'):
  if p.name.endswith('_progress.json'):continue
  d=read(p);prod=[f for f in d['frames'] if f['phase']=='production'];early=[f for f in prod if f['time_ps']<=d['warmup_ps']+4.];row={'candidate_id':d['candidate_id'],'seed':d['seed'],'warmup_ps':d['warmup_ps'],'production_ps':d['production_ps'],**d['summary'],'first4ps_production_loop_mean_A':sum(f['loop_CA_RMSD_A'] for f in early)/len(early),'legacy_geometry_pass':d['final_geometry_audit']['pass']};dynamics[p.stem]=row
 assert len(dynamics)==11 and 'S00022_v2_s198003' in dynamics,'Final extended proline test not finished'
 controls={}
 for cid,root in [('C00003',S/'intermediate/matched_controls')]:
  controls['matched_'+cid]=short_model(root/'intermediate/refined'/(cid+'_human6ARU_proton.json'))
 for probe,parent in [('S00062','C00003'),('S00063','S00028')]:
  root=S/'intermediate/histidine_matched_controls'/probe;controls['histidine_probe_parent_'+probe]=short_model(root/'intermediate/refined'/(parent+'_human6ARU_proton.json'));a=read(root/'intermediate/refined'/(parent+'_independent_v2.json'));controls['histidine_probe_parent_'+probe]['geometry_pass']=a['pass']
 sensitivity={}
 for label,cid,root in [('matched_parent','C00003',S/'intermediate/matched_controls'),('extra_disulfide','S00028',S)]:
  sensitivity[label]={'native54':short_model(root/'intermediate/refined'/(cid+'_human6ARU_proton.json')),'capped54':short_model(root/'intermediate/capped'/(cid+'_human6ARU_capped_proton.json')),'histidine_carboxylate162':short_model(root/'intermediate/acid_ensemble'/(cid+'_acid_histidine.json'))}
 acid={cid:short_model(S/'intermediate/refined'/(cid+'_human6ARU_acid_proton.json')) for cid in ['S00017','S00022','S00028']}
 ledger={}
 for cid,d in sorted(designs.items()):
  human=inv['human'].get(cid);mouse=inv['mouse'].get(cid);decision='Not deeply evaluated; no biological or ranking promotion.'
  if human:
   if not human['independent_geometry']['pass']:decision='Exclude current minimized geometry: peptide/disulfide strain gate.'
   elif cid in ['S00062','S00063']:decision='Do not promote F107H: conservative contrast is below equal-additional-refinement parent.'
   elif cid=='S00028':decision='Retain experimental hypothesis, but no demonstrated thermal-stability advantage over matched parent; no final submission status.'
   elif cid=='S00022':decision='Retain simple proline hypothesis for comparison; extended dynamics are diagnostic, not stability proof; no final submission status.'
   else:decision='Geometry-clean physical hypothesis; incomplete mouse/stability evaluation and unresolved neutral binding.'
  ledger[cid]={'sequence_sha256':hashlib.sha256(d['sequence'].encode()).hexdigest(),'human_detailed':bool(human),'minimized_geometry_pass':human['independent_geometry']['pass'] if human else None,'mouse_detailed':bool(mouse),'mouse_minimized_geometry_pass':mouse['independent_geometry']['pass'] if mouse else None,'measured_human_affinity':False,'measured_mouse_binding':False,'measured_acid_selectivity':False,'neutral_no_detection_established':False,'decision':decision}
 f={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_records':63,'new_unique_sequences':63,'new_length_range':[125,125],'new_unique_cdr3':len({d['cdr_sequences'][2] for d in designs.values()}),'all_retained_generation_counts':counts,'stage3_plus4_records':len(old),'stage3_plus4_unique':len(set(old.values())),'historical_all_stages_unique_including_rejected':len(historical),'historical_all_stages_records_including_rejected':sum(counts.values()),'human_detailed_unique_new':len(inv['human']),'human_geometry_clean_new':sum(r['independent_geometry']['pass'] for r in inv['human'].values()),'human_all54_positive_new':sum(r['acid_on_all54'] for r in inv['human'].values()),'mouse_detailed_unique_new':len(inv['mouse']),'mouse_geometry_clean_new':sum(r['independent_geometry']['pass'] for r in inv['mouse'].values()),'human':inv['human'],'mouse':inv['mouse'],'controls':controls,'acid_structural_hypotheses':acid,'sensitivity':sensitivity,'dynamics':dynamics,'corrected_trajectory_count':len(dynamics),'corrected_total_production_ps':sum(r['production_ps'] for r in dynamics.values()),'corrected_total_simulated_ps_including_warmup':sum(r['production_ps']+r['warmup_ps'] for r in dynamics.values()),'excluded_temperature_biased_pilot_count':3,'novelty':read(S/'reference/novelty_summary.json'),'blast':read(S/'reference/blast_coverage_summary.json'),'sequence_geometry_tests':read(S/'reference/stage4_unit_tests.json')['passed'],'integrator_tests':read(S/'reference/integrator_tests.json'),'proton_polynomial_tests':read(S/'reference/proton_linkage_bound_tests.json')['test_count'],'decision_ledger':ledger,'qualification':'All scores and trajectory summaries are conditional model diagnostics; no experimental pass, folding validation, KD, kinetic rate or novelty-upload guarantee. Ranked100 pending.'}
 assert f['human_detailed_unique_new']==14 and f['human_geometry_clean_new']==9
 (S/'output/final_evidence_inventory.json').write_text(json.dumps(f,indent=2));(S/'output/decision_ledger.json').write_text(json.dumps(ledger,indent=2))
 fasta=''.join('>'+cid+' NOT_FOR_SUBMISSION unvalidated_stage4\n'+d['sequence']+'\n' for cid,d in sorted(designs.items()))
 (S/'output/stage4_INTERNAL_NOT_FOR_SUBMISSION.fasta').write_text(fasta)
 print(json.dumps({k:v for k,v in f.items() if isinstance(v,(str,int,float,list))},indent=2))
if __name__=='__main__':main()
