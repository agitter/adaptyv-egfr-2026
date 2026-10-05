"""Conservative experimental-priority ranking, not calibrated affinity prediction.
Novelty does not use whole-chain identity as an antibody rejection rule. Real
Adaptyv antibody classification/CDRH3 annotation and eligibility remain pending.

Ranking predefinition: reject observed geometry/acid-direction failures unless a
valid same-sequence geometry repairs the failed geometry. Prefer wider protonation
scope and two-species all-atom checks. Within each evidence stratum use conservative
pH contrast first, modest redundancy penalty, then fewer sequence liabilities.
No candidate is labeled as meeting experimental criteria.
"""
from pathlib import Path
import json,hashlib,collections,re,math,sys
S=Path('/mnt/data/egfr_campaign/stage6');R=S.parent

def load(path,default=None):
 return json.loads(Path(path).read_text()) if Path(path).exists() else default

def make_pool():
 rows=load(S/'reference/campaign_audit.json',[])+load(S/'reference/new_design_audit.json',[])
 models=load(S/'reference/uniform_model_summary.json',[])+load(S/'reference/new_model_summary.json',[])
 byseq=collections.defaultdict(list);bm=collections.defaultdict(list);co=collections.defaultdict(list);expanded=collections.defaultdict(list);unbound=collections.defaultdict(list)
 for r in rows:byseq[r['sequence']].append(r)
 for m in models:bm[m['sequence_sha256']].append(m)
 for c in load(S/'reference/common_coarse_summary.json',{}).get('records',[]):co[c['sequence']].append(c)
 # Earlier full engineered-acid scope; do not promote old two-acid-only checks.
 ev=load(R/'stage5/output/final_evidence_inventory.json',{})
 for cid,r in ev.get('evidence_by_candidate',{}).items():
  for m in r.get('mixed_acid_His_tests',[]):
   expanded[r['sequence']].append({'source_stage':'stage5','candidate_id':cid,'source':m.get('source',str(R/'stage5/intermediate/independent_mixed_priors')),'minimum_contrast_kcal':m['minimum_contrast_kcal'],'method':'Exhaustive saved microstate energies and independent priors','numerical_validation_pass':True,'full_record':m})
 for p in sorted((S/'intermediate/cluster_acids').glob('*_cluster.json')):
  d=load(p);seq=load(S/'intermediate/designs'/f'{d["candidate_id"]}.json')['sequence'];expanded[seq].append({'source_stage':'stage6','candidate_id':d['candidate_id'],'source':str(p),'minimum_contrast_kcal':d['proton_model']['minimum_contrast_kcal'],'method':'Third-order classical cluster interpolation with held-out exact-energy checks','numerical_validation_pass':d['numerical_validation_pass'],'held_out_max_energy_error_kcal':d['held_out_max_absolute_energy_error_kcal'],'site_count':len(d['sites'])})
 for stage in ['stage3','stage4','stage5','stage6']:
  for p in (R/stage/'intermediate/unbound').glob('*_free_audit.json'):
   cid=p.name.split('_free_audit')[0];dp=R/stage/'intermediate/designs'/f'{cid}.json'
   if not dp.exists():continue
   seq=load(dp).get('sequence');d=load(p);ip=p.with_name(cid+'_free_independent.json');ind=load(ip);unbound[seq].append({'source':str(p),'geometry_pass':d['geometry_audit']['pass'] and (ind is None or ind['pass']),'loop_displacement_A':d['displacements']['loop_CA_RMSD_A'],'potential_relaxation_kcal':d['potential_relaxation_kcal'],'qualification':'Local minimization only, not a folding/stability measurement.'})
 eligible=[];excluded=[]
 for seq,group in byseq.items():
  h=hashlib.sha256(seq.encode()).hexdigest();mods=bm[h];hm=[m for m in mods if m['species']=='human'];mm=[m for m in mods if m['species']=='mouse'];hp=[m for m in hm if m['combined_geometry_pass']];mp=[m for m in mm if m['combined_geometry_pass']];nov=[r['local_novelty'] for r in group if r.get('local_novelty_pass')];reason=[]
  if not 10<=len(seq)<=250 or not set(seq)<=set('ACDEFGHIKLMNPQRSTVWY'):reason.append('format')
  if seq.count('C')!=2:reason.append('additional_cysteine_oxidation_uncertainty')
  if re.search('N[^P][ST]',seq):reason.append('glycosylation_sequon')
  if not nov:reason.append('local_CDR3_novelty_missing_or_fail')
  if not any(r.get('initial_context_pass') for r in group) and not (hp and mp):reason.append('no_passing_full_context_geometry')
  if hm and not hp:reason.append('unrepaired_human_geometry_failure')
  if mm and not mp:reason.append('unrepaired_mouse_geometry_failure')
  if any(not m['independent_priors'].get('all_acid_on',False) for m in hp):reason.append('valid_human_model_not_acid_on_across_priors')
  ex=[e for e in expanded[seq] if e['numerical_validation_pass']]
  if any(e['minimum_contrast_kcal']<=0 for e in ex):reason.append('expanded_acid_model_not_acid_on')
  if any(not u['geometry_pass'] for u in unbound[seq]):reason.append('unbound_local_geometry_failure')
  cg=[c for c in co[seq] if c['all_acid_on']]
  if not hp and not cg:reason.append('no_positive_human_pH_model')
  if any(m.get('acid_proxy',0)>=0 for m in hp+mp):reason.append('nonfavorable_declared_interaction_proxy')
  # Prefer a representative record with a valid refined human geometry; all
  # sequence aliases and failed structures remain attached to the evidence.
  validkeys={m['key'] for m in hp};cand=sorted(group,key=lambda r:(r['key'] not in validkeys,not r.get('initial_context_pass',False),r['key']))[0]
  if reason:excluded.append({'sequence_sha256':h,'keys':[r['key'] for r in group],'reasons':sorted(set(reason))});continue
  hi=min((m['independent_priors']['minimum_contrast_kcal'] for m in hp),default=None)
  if hp and mp and ex:tier=0;contrast=min([hi]+[e['minimum_contrast_kcal'] for e in ex]);label='Expanded protonation scope and human/mouse all-atom geometry'
  elif hp and mp:tier=1;contrast=hi;label='Human/mouse all-atom geometry; His-only pH model'
  elif hp:tier=2;contrast=hi;label='Human all-atom pH/geometry; mouse geometric prefilter only'
  else:tier=3;contrast=min(c['minimum_contrast'] for c in cg);label='Exploratory reserve: coarse pH model and initial geometry only'
  # Liabilities here are caution counts, not validated expression predictions.
  li=cand['sequence_audit']['liabilities'];pen=len(li.get('NG_motifs',[]))+len(li.get('DG_motifs',[]))+len(li.get('hydrophobic_runs',[]))
  eligible.append({'sequence':seq,'sequence_sha256':h,'candidate_key':cand['key'],'all_aliases':[r['key'] for r in group],'design_source':cand['design_source'],'design_sha256':cand['design_sha256'],'length':len(seq),'molecule_class':'nanobody','CDRs':cand['CDRs'],'intervals':cand['intervals'],'pose_id':cand['pose_id'],'branch':cand['branch'],'tier':tier,'evidence_label':label,'ranking_contrast_kcal_proxy':contrast,'histidine_only_minimum_contrast_kcal':hi,'expanded_protonation':ex,'human_models':hp,'mouse_models':mp,'failed_geometry_models':[m['model'] for m in mods if not m['combined_geometry_pass']],'initial_context_sources':[r['design_source'] for r in group if r.get('initial_context_pass')],'local_novelty':max(nov,key=lambda n:n['max_edit_identity']),'local_novelty_pass':True,'official_novelty_status':'PENDING_ADAPTYV','unbound_local_checks':unbound[seq],'coarse_models':co[seq],'sequence_liabilities':li,'liability_count':pen,'biological_criteria':{'human_binding_pH6p5':'UNVERIFIED','human_no_detectable_binding_pH7p4':'UNVERIFIED','mouse_binding':'UNVERIFIED','human_affinity':'UNVERIFIED','folding_expression':'UNVERIFIED'}})
 return eligible,excluded

def rank_pool(pool,limit=100):
 selected=[];cdr_counts=collections.Counter();pose_counts=collections.Counter();remaining=list(pool)
 # Limit repetition of one CDRH3 to six entries in the entire 100-candidate pool.
 # A modest in-stratum score penalty spreads very similar hypotheses. This is
 # an explicit panel-design preference, not a calibrated probabilistic score.
 for tier in range(4):
  part=[r for r in remaining if r['tier']==tier]
  while part and len(selected)<limit:
   allowed=[r for r in part if cdr_counts[r['CDRs'][2]]<6]
   if not allowed:break
   def order(r):
    diversity=.025*cdr_counts[r['CDRs'][2]]
    if tier==3:diversity+=.01*pose_counts[r['pose_id']]
    return (r['ranking_contrast_kcal_proxy']-diversity,-r['liability_count'],r['candidate_key'])
   best=max(allowed,key=order);best=dict(best);best['selection_adjustment']=.025*cdr_counts[best['CDRs'][2]]+(.01*pose_counts[best['pose_id']] if tier==3 else 0.);best['rank']=len(selected)+1;selected.append(best);cdr_counts[best['CDRs'][2]]+=1;pose_counts[best['pose_id']]+=1;part=[r for r in part if r['sequence']!=best['sequence']]
 return selected

def main():
 pool,excluded=make_pool();ranked=rank_pool(pool,125);(S/'reference/eligible_pool.json').write_text(json.dumps(pool,indent=2));(S/'reference/exclusion_ledger.json').write_text(json.dumps(excluded,indent=2));(S/'reference/panel_preview.json').write_text(json.dumps(ranked,indent=2));print('Eligible',len(pool),'excluded',len(excluded),'preview',len(ranked),'tiers',dict(collections.Counter(r['tier'] for r in ranked)))
 for r in ranked[:25]:print(r['rank'],r['candidate_key'],'tier',r['tier'],round(r['ranking_contrast_kcal_proxy'],4),r['CDRs'][2])
if __name__=='__main__':main()
