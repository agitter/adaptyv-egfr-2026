"""Independent release checks: file formats, source ancestry, actual coordinates.
These checks do not establish folding, binding, expression or official novelty.
"""
from pathlib import Path
import csv, json, hashlib, re, subprocess, sys, collections
R=Path('/mnt/data/egfr_campaign'); S=R/'stage7'; O=S/'output'
AA=set('ACDEFGHIKLMNPQRSTVWY')
AA3=dict(zip('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),'ACDEFGHIKLMNPQRSTVWY'))
AA3.update(HID='H',HIE='H',HIP='H',ASH='D',GLH='E',CYX='C',LYN='K')
checks=[]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def check(name,condition,details=None):
 row={'test':name,'pass':bool(condition)}
 if details is not None:row['details']=details
 checks.append(row)
 if not condition:raise AssertionError(row)
def parse_fasta(p):
 rows=[]; name=None; seq=[]
 for line in p.read_text(encoding='ascii').splitlines():
  if line.startswith('>'):
   if name is not None:rows.append((name,''.join(seq)))
   name=line[1:];seq=[]
  elif line:seq.append(line)
 if name is not None:rows.append((name,''.join(seq)))
 return rows
def pdb_seq(p,chain=None):
 seen=set();a=[]
 for line in p.read_text().splitlines():
  if line.startswith('ENDMDL'):break
  if not line.startswith(('ATOM  ','HETATM')):continue
  if line[12:16].strip()!='CA' or (chain is not None and line[21]!=chain):continue
  key=(line[21],line[22:27]);res=line[17:20]
  if key in seen:continue
  if res not in AA3:raise ValueError((str(p),res))
  seen.add(key);a.append(AA3[res])
 return ''.join(a)
def main():
 rows=parse_fasta(O/'egfr_ranked_100_v2.fasta');top=parse_fasta(O/'egfr_track3_top20_v2.fasta')
 panel=json.loads((O/'final_panel_private_v2.json').read_text());summary=json.loads((O/'revision_summary.json').read_text())
 check('exactly100FASTA',len(rows)==100)
 check('unique100identifiers',len(set(n for n,s in rows))==100)
 check('unique100sequences',len(set(s for n,s in rows))==100)
 check('top20FASTAorder',top==rows[:20])
 check('private_panel_order',[(r['public_name'],r['sequence']) for r in panel]==rows)
 for fname,count,want in [('egfr_track3_top20_v2.csv',20,rows[:20]),('egfr_100_v2_RESERVES_NOT_SINGLE_UPLOAD.csv',100,rows)]:
  with (O/fname).open(newline='') as f:
   reader=csv.DictReader(f);columns=reader.fieldnames;data=list(reader)
  check(fname+'_columns',columns==['name','sequence','molecule_class'])
  check(fname+'_count',len(data)==count)
  check(fname+'_identities_order',[(r['name'],r['sequence']) for r in data]==want)
  check(fname+'_declared_class_not_official_classification',all(r['molecule_class']=='nanobody' for r in data))
 old=dict(parse_fasta(R/'stage6/output/egfr_ranked_100.fasta'))
 coordinate_checks=0
 for i,(r,(name,seq)) in enumerate(zip(panel,rows),1):
  prefix=f'rank{i:03}'
  check(prefix+'_rank',r['rank']==i)
  check(prefix+'_alphabet_length_identifier',10<=len(seq)<=250 and set(seq)<=AA and re.fullmatch(r'[A-Za-z0-9_]+',name) is not None)
  check(prefix+'_sequence_hash',hashlib.sha256(seq.encode()).hexdigest()==r['sequence_sha256'])
  source=Path(r['source_file']);d=json.loads(source.read_text())
  check(prefix+'_source_hash',sha(source)==r['source_sha256'])
  check(prefix+'_actual_design_sequence',d['sequence']==seq==''.join(x['aa'] for x in d['structure']))
  check(prefix+'_CDR_consistency',[seq[a:b] for a,b in d['cdr_intervals_zero_based']]==[r['CDR1'],r['CDR2'],r['CDR3']])
  check(prefix+'_local_CDR3_safeguard',r['local_CDR3_max_edit_identity']<.70)
  check(prefix+'_unverified_status',r['official_novelty_status']=='PENDING_ADAPTYV' and r['experimental_status'].startswith('UNVERIFIED:'))
  if not r['new_this_revision']:
   check(prefix+'_retained_public_name_identity',name in old and old[name]==seq)
  else:
   check(prefix+'_no_new_sequon',re.search(r'N[^P][ST]',seq) is None)
   ev=r['new_evidence'];cid=d['candidate_id']
   check(prefix+'_human_mouse_gates',ev['human']['all_geometry_pass'] and ev['mouse']['all_geometry_pass'])
   for label in ['human6ARU','mouseAF']:
    pp=S/'intermediate/refined'/f'{cid}_{label}_relaxed.pdb'
    check(prefix+'_'+label+'_coordinate_identity',pdb_seq(pp,'B')==seq);coordinate_checks+=1
   if r.get('promotion_rule_pass'):
    check(prefix+'_promotion_requires_complete_thermal',ev['thermal_comparison']['descriptive_thermal_rule_pass'] and len(ev['unrestrained_thermal_replicas'])==2)
    check(prefix+'_promotion_requires_quenched_geometry',len(ev['quenched_thermal_endpoints'])==2 and all(v['basic_geometry']['pass'] and v['strict_geometry']['pass'] for v in ev['quenched_thermal_endpoints']))
    parent=json.loads((S/'intermediate/cluster_acids/M00000_cluster.json').read_text())
    check(prefix+'_promotion_requires_expanded_pH_retention',ev['expanded_protonation']['validation_pass'] and r['conservative_pH_proxy_kcal']>=parent['proton_model']['minimum_contrast_kcal']-.10)
   parent_design=json.loads((R/'stage6/intermediate/designs/J00008.json').read_text())
   check(prefix+'_framework_unchanged',all(x['loop'] or x['aa']==parent_design['structure'][j]['aa'] for j,x in enumerate(d['structure'])))
  if i<=20:check(prefix+'_top20_has_both_species_evidence',r['human_passing_models']>=1 and r['mouse_passing_models']>=1 and r['evidence_tier']<=1)
 check('CDR3representation_cap',max(collections.Counter(r['CDR3'] for r in panel).values())<=6)
 check('reported_CDR3_count_matches_panel',summary['unique_CDR3']==len(set(r['CDR3'] for r in panel)))
 for p in sorted((S/'intermediate/refined').glob('M*_relaxed.pdb')):
  cid=p.name.split('_')[0];d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text())
  check('all_models_'+p.name,pdb_seq(p,'B')==d['sequence']);coordinate_checks+=1
 with (O/'private_codebook_v2.tsv').open(newline='') as f:book=list(csv.DictReader(f,delimiter='\t'))
 check('private_codebook_order',[(r['public_name'],r['sequence_sha256']) for r in book]==[(r['public_name'],r['sequence_sha256']) for r in panel])
 # Rerun deterministic selection. This writes only stage7 release outputs.
 outputs=['egfr_ranked_100_v2.fasta','egfr_track3_top20_v2.fasta','egfr_track3_top20_v2.csv','egfr_100_v2_RESERVES_NOT_SINGLE_UPLOAD.csv','final_panel_private_v2.json','private_codebook_v2.tsv','revision_summary.json']
 before={n:sha(O/n) for n in outputs}
 run=subprocess.run([sys.executable,str(S/'code/release_revision.py')],capture_output=True,text=True)
 (S/'logs/postflight_selection_rerun.log').write_text(run.stdout+'\n'+run.stderr)
 check('selection_rerun_succeeded',run.returncode==0)
 for n in outputs:check('reproducible_'+n,before[n]==sha(O/n))
 out=S/'intermediate/formatter_tests/actual_panel_top20.csv'
 run=subprocess.run([sys.executable,str(S/'code/make_submission_v2.py'),'--fasta',str(O/'egfr_ranked_100_v2.fasta'),'--output',str(out)],capture_output=True,text=True)
 check('actual_formatter_returncode',run.returncode==0)
 check('actual_formatter_agrees_with_release',out.read_bytes()==(O/'egfr_track3_top20_v2.csv').read_bytes())
 for file in ['input_validation.json','polynomial_validation.json','alignment_tests.json']:
  z=json.loads((S/'reference'/file).read_text());check(file+'_prior_tests_pass',z['all_pass'])
 z=json.loads((S/'reference/formatter_tests.json').read_text());check('formatter13tests',z['passed']==z['count']==13)
 follow=json.loads((S/'reference/followup_status.json').read_text());check('all_planned_followups_finished',follow['all_finished'] and len(follow['completed'])==follow['expected'] and all(x['returncode']==0 and not x['timed_out'] for x in follow['completed']))
 q=json.loads((S/'reference/endpoint_quench_status.json').read_text());check('all12matched_endpoint_quenches_finished',q['all_finished'] and len(q['completed'])==q['expected']==12 and all(x['returncode']==0 and not x['timed_out'] for x in q['completed']))
 history=json.loads((S/'reference/historical_sequence_uniqueness.json').read_text());check('22sequences_new_to_retained_campaign',history['new_unique_sequences']==22)
 report={'passed':len(checks),'total':len(checks),'all_pass':all(x['pass'] for x in checks),'coordinate_sequence_checks':coordinate_checks,'records':checks,'public_output_hashes':{n:sha(O/n) for n in outputs},'qualification':'File-format, ancestry, coordinate-identity and implementation tests only. No binding, folding, expression, pH-switch or official novelty result.'}
 (S/'reference/release_validation_v2.json').write_text(json.dumps(report,indent=2));print('RELEASE VALIDATION',len(checks),'PASSED; coordinate checks',coordinate_checks)
if __name__=='__main__':main()
