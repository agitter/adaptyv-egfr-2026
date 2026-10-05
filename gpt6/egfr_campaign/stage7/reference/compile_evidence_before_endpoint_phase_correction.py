"""Join observed stage7 results without converting missing evidence to passes."""
from pathlib import Path
import json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 structures=json.loads((S/'reference/structure_summary.json').read_text());by={(r['candidate_id'],r['species']):r for r in structures}
 nov={r['candidate_id']:r for r in json.loads((S/'intermediate/novelty/stage7_cdr3_screen.json').read_text())}
 tpath=S/'reference/thermal_analysis.json';thermal=json.loads(tpath.read_text()) if tpath.exists() else {'replicas':[],'comparisons_without_extra_restraints':[]}
 tby={r['candidate_id']:r for r in thermal['comparisons_without_extra_restraints']};rows=[]
 endpoint=S/'reference/thermal_final_geometry.json';endpoint_rows=json.loads(endpoint.read_text())['records'] if endpoint.exists() else []
 for p in sorted((S/'intermediate/designs').glob('M*.json')):
  d=json.loads(p.read_text());cid=d['candidate_id'];cluster=S/'intermediate/cluster_acids'/f'{cid}_cluster.json';c=json.loads(cluster.read_text()) if cluster.exists() else None
  u=S/'intermediate/unbound'/f'{cid}_free_audit.json';uq=S/'intermediate/unbound'/f'{cid}_free_independent.json'
  energy={}
  for species,label in [('human','human6ARU'),('mouse','mouseAF')]:
   pp=S/'intermediate/refined'/f'{cid}_{label}_proton.json'
   if pp.exists():
    pr=json.loads(pp.read_text());ss=pr['proton_model']['scenarios'];energy[species]={'least_favorable_acid_proxy_kcal':max(x['low_pH']['conditional_binding_energy_kcal'] for x in ss),'most_favorable_acid_proxy_kcal':min(x['low_pH']['conditional_binding_energy_kcal'] for x in ss),'least_favorable_neutral_proxy_kcal':max(x['high_pH']['conditional_binding_energy_kcal'] for x in ss),'scope':'Narrow histidine model, not calibrated affinity'}
  rows.append({'candidate_id':cid,'sequence':d['sequence'],'sequence_sha256':hashlib.sha256(d['sequence'].encode()).hexdigest(),'source_file':str(p),'source_sha256':sha(p),'cdr_sequences':d['cdr_sequences'],'mutations':d['replacement_experiment']['mutations'],'unchanged_control':d['replacement_experiment']['unchanged_control'],'human':by.get((cid,'human')),'mouse':by.get((cid,'mouse')),'narrow_energy_proxies':energy,'novelty':nov[cid],'expanded_protonation':({'source':str(cluster),'source_sha256':sha(cluster),'site_count':len(c['sites']),'sites':c['sites'],'microstate_count':c['state_count'],'validation_pass':c['numerical_validation_pass'],'max_heldout_energy_error_kcal':c['held_out_max_absolute_energy_error_kcal'],'proton_model':c['proton_model']} if c else None),'unbound_minimization':json.loads(u.read_text()) if u.exists() else None,'unbound_strict_geometry':json.loads(uq.read_text()) if uq.exists() else None,'thermal_comparison':tby.get(cid),'unrestrained_thermal_replicas':[r for r in thermal['replicas'] if r['candidate_id']==cid and not r['added_chirality_restraints']],'independent_thermal_endpoints':[r for r in endpoint_rows if r['candidate_id']==cid and not r['added_chirality_restraints']],'official_novelty':'PENDING_ADAPTYV','experimental_status':'UNVERIFIED: human pH6.5 binding, human pH7.4 nondetection, mouse binding, affinity, folding/expression'})
 out={'records':rows,'counts':{'sequence_records':len(rows),'new_sequences':sum(not r['unchanged_control'] for r in rows),'human_models':sum(r['human'] is not None for r in rows),'mouse_models':sum(r['mouse'] is not None for r in rows),'both_species_geometry_pass':sum(bool(r['human'] and r['mouse'] and r['human']['all_geometry_pass'] and r['mouse']['all_geometry_pass']) for r in rows),'expanded_models':sum(r['expanded_protonation'] is not None for r in rows),'completed_primary_thermal_replicas':sum(len(r['unrestrained_thermal_replicas']) for r in rows)},'all_biological_objectives_unverified':True}
 (S/'reference/evidence_compilation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out['counts'],indent=2))
if __name__=='__main__':main()
