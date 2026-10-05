"""Independent completed-work audit before packaging. Does not validate binding."""
from pathlib import Path
import json,hashlib,os,ast,datetime
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
def read(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def workers():
 found=[]
 for entry in Path('/proc').iterdir():
  if not entry.name.isdigit() or int(entry.name)==os.getpid():continue
  try:a=[b.decode(errors='replace') for b in (entry/'cmdline').read_bytes().split(b'\0') if b]
  except (FileNotFoundError,PermissionError,ProcessLookupError):continue
  if a and ((Path(a[0]).name.startswith('python') and any('/egfr_campaign/' in s and '/code/' in s for s in a[1:])) or (Path(a[0]).name in ['blastp','makeblastdb'] and any('egfr_campaign' in s for s in a[1:]))):found.append({'pid':int(entry.name),'args':a})
 return found
def main():
 checks=[]
 def check(name,ok):
  checks.append({'check':name,'pass':bool(ok)});assert ok,name
 syntax=[]
 for p in sorted((S/'code').glob('*.py')):
  ast.parse(p.read_text());syntax.append(source(p))
 e=read(S/'output/final_evidence_inventory.json');nv=read(S/'reference/novelty_complete.json');cat=read(S/'reference/candidate_catalog.json');tests=read(S/'reference/stage5_tests.json');dyn=read(S/'reference/dynamics_analysis.json')
 check('all1618stage5checks',tests['all_pass'] and tests['passed']==1618)
 check('144independentHischecks',read(S/'reference/independent_priors_tests.json')['passed']==144)
 check('144glycanbondchecks',e['glycan_flexibility']['bond_preservation_checks']==144 and e['glycan_flexibility']['all_bond_preservation_checks_pass'])
 check('276uniquequeries',len(cat)==len({r['candidate_id'] for r in cat})==len({r['sequence'] for r in cat})==276)
 lines=(S/'output/stage5_276_NOT_FOR_SUBMISSION.fasta').read_text().splitlines();seqs={}
 for i in range(0,len(lines),2):
  h=lines[i];seq=lines[i+1];check('FASTA label '+h.split()[0],h.startswith('>') and 'NOT_FOR_SUBMISSION' in h);seqs[h[1:].split()[0]]=seq
 check('FASTA exactcatalog',seqs=={r['candidate_id']:r['sequence'] for r in cat})
 for r in cat:check('catalog sourcehash '+r['candidate_id'],sha(r['source'])==r['source_sha256'])
 for db in nv['databases']:
  check('queryset '+db['database'],db['candidate_queries']==276 and db['all_ids_verified'])
  check('exactcontrols '+db['database'],len(db['positive_controls'])==3 and all(c['exact_control_pass'] for c in db['positive_controls']))
  for s in db['sources']:
   check('BLASTsummaryhash '+Path(s['summary']).name,sha(s['summary'])==s['summary_sha256']);check('BLASTresultshash '+Path(s['results']).name,sha(s['results'])==s['results_sha256'])
 check('CDRfailure retained',nv['cdr3_rejected_ids']==['N3081002'] and nv['cdr3_local_pass']==275)
 check('8trajectories completed',dyn['complete'] and len(dyn['per_trajectory'])==8 and dyn['aggregate_production_ps']==160 and dyn['aggregate_warmup_plus_production_ps']==176)
 check('thermal integrities',all(t['geometry_pass'] for t in dyn['per_trajectory']))
 check('two declared seeds each',all(a['independent_trajectory_count']==2 and a['seeds']==[198005,198006] for a in dyn['aggregates']))
 for cid in ['B00000','H00011','H00020','H00021']:
  r=e['evidence_by_candidate'][cid];check(cid+' models',len(r['models'])==3)
  for m in r['models']:
   for key in ['source','independent_prior_source']:check(cid+' modelhash '+Path(m[key]['path']).name,sha(m[key]['path'])==m[key]['sha256'])
  for a in r['mixed_acid_His_tests']:
   check(cid+' mixedsourcehash',sha(a['source'])==a['source_sha256']);check(cid+' scenarioarrayhash',sha(a['array_file'])==a['array_sha256']);check(cid+' exactsubset',a['previous_tied_prior_maximum_error_kcal']<1e-9)
  check(cid+' unboundpass',r['unbound_local_minimization']['geometry_audit']['pass'])
  check(cid+' notexperimentallypassed',set(r['experimental_criteria'].values())=={'UNVERIFIED'})
 check('H20mousefailedgate retained',any(m['species']=='mouse' and not m['geometry_context_pass'] for m in e['evidence_by_candidate']['H00020']['models']))
 batch=[]
 for p in sorted((S/'reference').glob('*_batch_execution.json')):
  d=read(p);check('batch finished '+p.name,d['all_finished']);failed=[r for r in d['tasks'] if r['returncode']!=0]
  if failed:check('only recovered BLAST batch fails',p.name=='blast_batch_execution.json' and {r['log'] for r in failed}=={'blast_search_pdb.log','blast_search_antibodies_augmented.log'})
  batch.append({'source':source(p),'failed_historical_tasks':failed,'recovered_by_complete_novelty_proofs':bool(failed)})
 # Directly launched exhaustiveH11 worker has no launcher exit metadata. Verify
 # its complete state table, completed summary, chargechecks and independentarray.
 hp=S/'intermediate/chemistry_acid_scope/intermediate/acid_ensemble/H00011_acid_histidine.json';h=read(hp)
 check('H11all729microstates',h['microstate_count']==len(h['microstates'])==729 and len({tuple(r['state']) for r in h['microstates']})==729)
 check('H11chargeincrements',all(r['pass'] for r in h['charge_increment_tests']))
 w=workers();check('no campaign workers',not w)
 result={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'checks':checks,'checks_passed':len(checks),'source_files_syntax_checked':len(syntax),'source_files':syntax,'campaign_workers_running':w,'completed_batch_manifests':batch,'artifact_status':'NOT_FOR_SUBMISSION checkpoint; final100pending','physical_and_biological_validation':False}
 (S/'reference/postflight.json').write_text(json.dumps(result,indent=2));print('POSTFLIGHT',len(checks),'checks passed;',len(syntax),'Python sources parsed; no campaign workers')
def source(p):return {'path':str(p),'sha256':sha(p)}
if __name__=='__main__':main()
