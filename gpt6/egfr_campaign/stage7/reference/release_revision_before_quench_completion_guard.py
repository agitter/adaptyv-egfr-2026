"""Build the revised experimental-priority pool only after all planned checks.
No experimental function or official novelty status is inferred by this export.
Retained public identifiers stay stable; rank is defined by file order.
"""
from pathlib import Path
import csv,json,hashlib,re,collections,math,importlib.util,shutil
from openpyxl import Workbook
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';OUT=S/'output'
spec=importlib.util.spec_from_file_location('fasta_formatter',R/'stage6/code/make_submission.py');formatter=importlib.util.module_from_spec(spec);spec.loader.exec_module(formatter)
NAMES={'M00003':'MolarMirage_Obsidian','M00008':'AmideAlibi_Kestrel','M00014':'PeptideParadox_Vesper','M00015':'CatalyticCaper_Quartz','M00018':'BondVoyage_Cipher','M00019':'ChiralCover_Rook','M00020':'EnzymeEnigma_Sable','M00021':'SilentCatalyst_Umbra','M00022':'IsoelectricInk_Aster'}
EXP='UNVERIFIED: human acidic binding, neutral-pH nondetection, mouse binding, affinity, folding/expression'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(path,rows):
 with path.open('w',encoding='ascii',newline='\n') as f:
  for r in rows:
   f.write('>'+r['public_name']+'\n');seq=r['sequence']
   for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')
def csv_export(path,rows):
 wb=Workbook();ws=wb.active;ws.title='Submission';ws.append(['name','sequence','molecule_class'])
 for r in rows:ws.append([r['public_name'],r['sequence'],'nanobody'])
 with path.open('w',encoding='utf-8',newline='') as f:csv.writer(f,lineterminator='\n').writerows(ws.iter_rows(values_only=True))
 wb.close()
def main():
 policy=json.loads((S/'reference/release_revision_policy.json').read_text());done=json.loads((S/'reference/followup_status.json').read_text());thermal=json.loads((S/'reference/thermal_analysis.json').read_text());evidence=json.loads((S/'reference/evidence_compilation.json').read_text());by={r['candidate_id']:r for r in evidence['records']}
 assert done['all_finished'] and len(done['completed'])==done['expected'] and all(r['returncode']==0 and not r['timed_out'] for r in done['completed']), 'Planned calculations are not all complete'
 assert thermal['all_expected_unrestrained_complete'], 'Eight complete primary thermal replicas required'
 for cid in ['M00000','M00018','M00021']:
  x=by[cid]['expanded_protonation'];assert x and x['validation_pass'] and x['microstate_count']==6561
 for p in ['input_validation.json','polynomial_validation.json']:
  assert json.loads((S/'reference'/p).read_text())['all_pass']
 oldmeta={r['public_name']:r for r in csv.DictReader((R/'stage6/output/private_codebook_and_evidence.tsv').open(),delimiter='\t')};old=[]
 for name,seq in formatter.read_fasta(R/'stage6/output/egfr_ranked_100.fasta'):
  m=oldmeta[name];assert hashlib.sha256(seq.encode()).hexdigest()==m['sequence_sha256']
  source_stage,cid=m['internal_key'].split(':');source=R/source_stage/'intermediate/designs'/f'{cid}.json';design=json.loads(source.read_text());assert seq==design['sequence']
  old.append({'public_name':name,'sequence':seq,'sequence_sha256':m['sequence_sha256'],'source_internal_key':m['internal_key'],'source_file':str(source),'source_sha256':sha(source),'previous_rank':int(m['rank']),'new_this_revision':False,'evidence_tier':int(m['evidence_tier']),'evidence_label':m['evidence_label'],'conservative_pH_proxy_kcal':float(m['conservative_pH_proxy_kcal']) if m['conservative_pH_proxy_kcal'] else None,'CDR1':m['CDR1'],'CDR2':m['CDR2'],'CDR3':m['CDR3'],'human_passing_models':int(m['human_passing_models']),'mouse_passing_models':int(m['mouse_passing_models']),'local_CDR3_max_edit_identity':float(m['local_CDR3_max_edit_identity']),'glycan_unweighted_overlap_fraction':float(m['glycan_unweighted_overlap_fraction']) if m['glycan_unweighted_overlap_fraction'] else None,'official_novelty_status':'PENDING_ADAPTYV','experimental_status':EXP,'prior_evidence':m})
 parent=by['M00000'];parent_min=parent['expanded_protonation']['proton_model']['minimum_contrast_kcal'];new=[];excluded=[]
 for cid in policy['new_sequences_considered']:
  e=by[cid];x=e['expanded_protonation'];reason=[]
  if not(e['human'] and e['mouse'] and e['human']['all_geometry_pass'] and e['mouse']['all_geometry_pass']):reason.append('missing or failed both-species structural gate')
  if not e['novelty']['below_70pct']:reason.append('local CDR3 safeguard')
  if e['human'] and e['human']['minimum_histidine_contrast']<=0:reason.append('contrary detailed human pH direction')
  if x and (not x['validation_pass'] or not x['proton_model']['all_scenarios_acid_on']):reason.append('expanded physical model contrary or numerically unvalidated')
  if e['unbound_strict_geometry'] and not e['unbound_strict_geometry']['pass']:reason.append('failed independent unbound geometry')
  if any(not t['final_geometry_audit']['pass'] for t in e['unrestrained_thermal_replicas']):reason.append('invalid final geometry in primary thermal test')
  if e.get('unrestrained_thermal_replicas') and (len(e.get('quenched_thermal_endpoints',[]))!=2 or any(not t['strict_geometry']['pass'] or not t['basic_geometry']['pass'] for t in e.get('quenched_thermal_endpoints',[]))):reason.append('failed or missing phase-matched quenched endpoint geometry')
  if reason:excluded.append({'candidate_id':cid,'reasons':reason});continue
  ph=x['proton_model']['minimum_contrast_kcal'] if x else e['human']['minimum_histidine_contrast'];tc=e['thermal_comparison'];promoted=bool(tc and tc['descriptive_thermal_rule_pass'] and x and ph>=parent_min-.10)
  r={'public_name':NAMES[cid],'sequence':e['sequence'],'sequence_sha256':e['sequence_sha256'],'source_internal_key':'stage7:'+cid,'source_file':e['source_file'],'source_sha256':e['source_sha256'],'previous_rank':None,'new_this_revision':True,'evidence_tier':0 if x else 1,'evidence_label':'Expanded protonation and both-species geometry' if x else 'Narrow histidine analysis and both-species geometry','conservative_pH_proxy_kcal':ph,'CDR1':e['cdr_sequences'][0],'CDR2':e['cdr_sequences'][1],'CDR3':e['cdr_sequences'][2],'human_passing_models':1,'mouse_passing_models':1,'local_CDR3_max_edit_identity':e['novelty']['max_edit_identity'],'glycan_unweighted_overlap_fraction':e['human']['glycan_fraction'],'official_novelty_status':'PENDING_ADAPTYV','experimental_status':EXP,'promotion_rule_pass':promoted,'thermal_comparison':tc,'new_evidence':e}
  assert len(r['sequence'])==125 and not re.search(r'N[^P][ST]',r['sequence']);new.append(r)
 promoted=sorted([r for r in new if r['promotion_rule_pass']],key=lambda r:(-r['conservative_pH_proxy_kcal'],r['source_internal_key']))
 reserves=sorted([r for r in new if not r['promotion_rule_pass']],key=lambda r:(r['evidence_tier'],-r['conservative_pH_proxy_kcal'],r['source_internal_key']))
 order=promoted+old[:20]+reserves+old[20:];selected=[];seen=set();cdr_counts=collections.Counter();omitted=[]
 for r in order:
  reason=None
  if r['sequence'] in seen:reason='duplicate sequence'
  elif cdr_counts[r['CDR3']]>=6:reason='CDR3 representation cap'
  elif len(selected)>=100:reason='lower-priority reserve beyond100'
  if reason:omitted.append({'public_name':r['public_name'],'source_internal_key':r['source_internal_key'],'reason':reason});continue
  r=dict(r);r['rank']=len(selected)+1;selected.append(r);seen.add(r['sequence']);cdr_counts[r['CDR3']]+=1
 assert len(selected)==100 and len({r['public_name'] for r in selected})==100
 for r in selected:
  assert 10<=len(r['sequence'])<=250 and set(r['sequence'])<=set('ACDEFGHIKLMNPQRSTVWY')
  assert re.fullmatch('[A-Za-z0-9_]+',r['public_name'])
 assert all(r['human_passing_models']>0 and r['mouse_passing_models']>0 and r['evidence_tier']<=1 for r in selected[:20])
 OUT.mkdir(parents=True,exist_ok=True)
 fasta(OUT/'egfr_ranked_100_v2.fasta',selected);fasta(OUT/'egfr_track3_top20_v2.fasta',selected[:20]);csv_export(OUT/'egfr_track3_top20_v2.csv',selected[:20]);csv_export(OUT/'egfr_100_v2_RESERVES_NOT_SINGLE_UPLOAD.csv',selected)
 (OUT/'final_panel_private_v2.json').write_text(json.dumps(selected,indent=2))
 fields=['rank','public_name','previous_rank','source_internal_key','new_this_revision','sequence_sha256','length','evidence_tier','evidence_label','conservative_pH_proxy_kcal','human_passing_models','mouse_passing_models','local_CDR3_max_edit_identity','CDR1','CDR2','CDR3','glycan_unweighted_overlap_fraction','promotion_rule_pass','official_novelty_status','experimental_status','source_file','source_sha256']
 with (OUT/'private_codebook_v2.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in selected:
   q={k:r.get(k,'') for k in fields};q['length']=len(r['sequence']);w.writerow(q)
 summary={'panel_size':100,'top20_size':20,'new_sequences_in_pool':sum(r['new_this_revision'] for r in selected),'new_sequences_in_top20':sum(r['new_this_revision'] for r in selected[:20]),'promoted_by_complete_comparison':[r['source_internal_key'] for r in promoted],'top20_order_changed':[(r['public_name'],r['sequence']) for r in selected[:20]]!=[(r['public_name'],r['sequence']) for r in old[:20]],'unique_CDR3':len(cdr_counts),'max_CDR3_repetition':max(cdr_counts.values()),'length_min':min(len(r['sequence']) for r in selected),'length_max':max(len(r['sequence']) for r in selected),'evidence_counts':dict(collections.Counter(r['evidence_tier'] for r in selected)),'new_entries_considered_but_excluded':excluded,'omitted_by_pool_size_or_diversity':omitted,'old_release_sha256':sha(R/'stage6/output/egfr_ranked_100.fasta'),'new_fasta_sha256':sha(OUT/'egfr_ranked_100_v2.fasta'),'new_csv_sha256':sha(OUT/'egfr_track3_top20_v2.csv'),'rank_is_file_order':True,'old_public_identifiers_preserved':True,'all_experimental_outcomes_unverified':True}
 (OUT/'revision_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
