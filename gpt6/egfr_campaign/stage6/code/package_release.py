"""Package completed campaign files and independently verify every ZIP payload."""
from pathlib import Path
import json,hashlib,zipfile,shutil,os,sys,time,datetime,shlex,subprocess
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';O=S/'output';D=Path('/mnt/data')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def check_workers():
 bad=[]
 for line in subprocess.check_output(['ps','-eo','pid,args'],text=True).splitlines()[1:]:
  try:pidtext,command=line.strip().split(None,1);pid=int(pidtext);args=shlex.split(command)
  except (ValueError,IndexError):continue
  if pid==os.getpid() or not args:continue
  if Path(args[0]).name.startswith('python') and len(args)>1 and args[1].startswith(str(R)) and '/code/' in args[1]:bad.append([pid,args])
  elif '/egfr_campaign/' in args[0] and Path(args[0]).name in ['blastp','makeblastdb']:bad.append([pid,args])
 assert not bad,bad

def main():
 start=time.time();check_workers();post=json.loads((S/'reference/final_postflight.json').read_text());assert post['checks']==post['passed'] and not post['active_workers']
 restore=json.loads((S/'reference/restoration.json').read_text());archives={Path(r['path']).name:r for r in restore['archives']}
 inputs=[]
 for name,dest in [('resources.zip','inputs/provided'),('egfr_inputs.zip','inputs/transfer1'),('egfr_refinement_inputs.zip','inputs/transfer2')]:inputs.append({'filename':name,'sha256':archives[name]['sha256'],'destination':dest})
 (S/'reference/final_checkpoint_inputs.json').write_text(json.dumps({'required_original_archives':inputs,'earlier_checkpoint_provenance':archives['egfr_stage5_checkpoint.zip'],'earlier_checkpoint_required_to_restore_this_archive':False,'runtime_and_packages':{'python':restore['python'],'platform':restore['platform'],'packages':restore['packages']},'original_input_payloads_not_duplicated':True},indent=2))
 if not (S/'reference/root_state_before_final_release.md').exists():shutil.copy2(R/'STATE.md',S/'reference/root_state_before_final_release.md')
 (R/'STATE.md').write_text('# EGFR campaign: ranked 100-sequence release completed\n\nCurrent artifacts: stage6/output/. Read stage6/output/README.md, CAMPAIGN_REPORT.md, METHODS_CERTIFICATE.md, and release_validation.json.\n\nThe 100-candidate FASTA and ordered Track-3 top20 CSV are final files from this campaign. Experimental binding, pH selectivity, mouse binding, expression/folding, and official novelty remain unverified. No submission was performed. Two Phone-a-Friend requests used, one unused.\n\nPrior stages, failures, and corrected-source versions remain preserved. Runtime and original large inputs are restored separately; see stage6/RESTORE.md.\n')
 (S/'STATE.md').write_text('# Stage 6: completed release\n\n100 unique ranked sequences; 119-127 aa; public FASTA and top20 CSV in output/. The top20 have checked human/mouse all-atom geometries, not experimentally confirmed binding.\n\nAll scientific workers and sequence searches finished. Ten original unbound wrapper exits were postprocessing errors with retained logs and successful independent chain-label recovery.\n\nDo not overwrite the archived observations while reproducing. The final checkpoint manifest covers both prior-stage payloads and new stage6 code/results.\n')
 (S/'RESTORE.md').write_text('''# Restore and inspect the final campaign checkpoint

Extract `egfr_campaign_final_checkpoint.zip` into `/mnt/data` on Linux. The public FASTA/CSV and private evidence are already present under `egfr_campaign/stage6/output`; inspection does not require running the scientific code.

Verify the archived payloads with Python's standard library:

```bash
python /mnt/data/egfr_campaign/stage6/code/restore_final.py
```

For numerical work, retain the original user-supplied `resources.zip`, `egfr_inputs.zip`, and `egfr_refinement_inputs.zip`. Their required SHA-256 hashes are recorded in `reference/final_checkpoint_inputs.json`. Restore their payloads and the supplied OpenMM wheel using:

```bash
python /mnt/data/egfr_campaign/stage6/code/restore_final.py --restore-inputs
```

The numerical environment was Linux x86_64, Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, Biopython 1.86, Numba 0.65.1, RapidFuzz 3.14.3, and the supplied OpenMM 8.4.0.post2 wheel. Modern implementations were permitted by the user. Do not substitute a modern learned potential or protein generator. Scientific scripts use `/mnt/data/egfr_campaign` paths. Set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, and OPENMM_CPU_THREADS=1. The observed run had eight CPU cores of quota and 4 GiB RAM.

BLAST indexes and normalized input FASTAs are regenerable with `stage6/code/build_databases6.py`; they and unpacked executables are not duplicated here. Final search queries/results, controls, hashes, and source database counts are retained. The supplied mouse structure is a modern input dataset; no structure prediction was run locally.

For reproduction use a separate clean working copy. The exact job plans, commands, seeds, source ancestry, and completed statuses are in `stage6/reference/`. Saved-result caches intentionally prevent silent overwriting. The observed outputs are the provenance record; a rerun can change floating-point minimization results and is not promised to be bit-for-bit identical.

The computational order for this final stage is: restoration and prior-hash checks; inventory and uniform old-model audit; matched J contact combinations and human/mouse calculations; K alternative-pose checks and L mouse completions; expanded cluster/protonation analysis and its independent tests; unbound local diagnostics with the documented postprocessing correction; local novelty query reconciliation; complete glycan screening; evidence-stratified selection; FASTA/CSV export; independent final postflight.

Read the earlier stage reports for ancestry before stage6. Failed/obsolete script variants are historical records, not an instruction to execute them against the final outputs. `restore_stage6.py` is the original pre-stage6 restoration script; use `restore_final.py` for this final snapshot because the root state file has intentionally changed.

The checkpoint archive has a complete included-payload manifest and a separate list of omitted regenerable stage6 files. The original input archives are external dependencies, not missing computed results. The independent archive verification is delivered alongside the ZIP. The excluded modern-RNG pilot is preserved and explicitly disclaimed in METHODS_CERTIFICATE.md.
''')
 (S/'reference/refresh_timeout_record.json').write_text(json.dumps({'operation':'Intermediate audit/glycan/selection refresh under a 45-second container command limit','outcome':'Command timed out after intermediate files were produced; no biological conclusion taken from completion status.','recovery':'Completed the full refresh under a 120-second command limit, then independently validated final source sequences, gates, query coverage, and idempotent selection.','final_evidence':str(S/'reference/final_postflight.json')},indent=2))
 # Preserve every file in the complete stage5 checkpoint, plus new stage6 data.
 base=D/'egfr_stage5_checkpoint.zip'
 with zipfile.ZipFile(base) as z:paths={x.filename[len('egfr_campaign/'):] for x in z.infolist() if not x.is_dir() and x.filename.startswith('egfr_campaign/')}
 excluded=[];special={'stage6/reference/final_checkpoint_included_manifest.json','stage6/reference/final_checkpoint_excluded_manifest.json'}
 for p in S.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(R).as_posix();sub=p.relative_to(S).as_posix();reason=None
  if rel in special:continue
  if '__pycache__' in p.parts or p.suffix in ['.pyc','.nbc','.nbi']:reason='Regenerable interpreter/JIT cache'
  elif sub.startswith('runtime/'):reason='Unpacked software restored from supplied archives'
  elif sub.startswith('intermediate/blast/') and (p.name.startswith(('swissprot.','pdb.','antibodies_augmented.')) or p.name=='antibody_sources.jsonl'):reason='Regenerable reference FASTA, BLAST index, or reference normalization'
  if reason:excluded.append({'path':rel,'bytes':p.stat().st_size,'reason':reason})
  else:paths.add(rel)
 included=[]
 for rel in sorted(paths):
  p=R/rel;assert p.is_file(),rel;included.append({'path':rel,'bytes':p.stat().st_size,'sha256':sha(p)})
 im=S/'reference/final_checkpoint_included_manifest.json';em=S/'reference/final_checkpoint_excluded_manifest.json';im.write_text(json.dumps(included,indent=2));em.write_text(json.dumps(excluded,indent=2))
 target=D/'egfr_campaign_final_checkpoint.zip';temporary=D/'egfr_campaign_final_checkpoint.zip.tmp'
 with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for r in included:z.write(R/r['path'],'egfr_campaign/'+r['path'])
  for p in [im,em]:z.write(p,'egfr_campaign/'+p.relative_to(R).as_posix())
 os.replace(temporary,target);print('ARCHIVED',len(included),'payloads',target.stat().st_size,'bytes',flush=True)
 # Independently hash every decompressed ZIP payload, not just archive CRCs.
 failures=[]
 with zipfile.ZipFile(target) as z:
  names=set(z.namelist());expected={'egfr_campaign/'+r['path'] for r in included}|{'egfr_campaign/'+p.relative_to(R).as_posix() for p in [im,em]};assert names==expected
  stored=json.loads(z.read('egfr_campaign/stage6/reference/final_checkpoint_included_manifest.json'));assert stored==included
  for row in stored:
   name='egfr_campaign/'+row['path'];h=hashlib.sha256();n=0
   with z.open(name) as f:
    for b in iter(lambda:f.read(1024*1024),b''):h.update(b);n+=len(b)
   if n!=row['bytes'] or h.hexdigest()!=row['sha256']:failures.append(row['path'])
 assert not failures,failures
 # Convenient standalone files have the same hashes as the archived originals.
 copies={'egfr_ranked_100.fasta':'egfr_ranked_100.fasta','egfr_track3_top20.csv':'egfr_track3_top20.csv','egfr_track3_top20.fasta':'egfr_track3_top20.fasta','CAMPAIGN_REPORT.md':'egfr_campaign_report.md','METHODS_CERTIFICATE.md':'egfr_methods_certificate.md','private_codebook_and_evidence.tsv':'egfr_private_codebook.tsv','release_validation.json':'egfr_release_validation.json'}
 for source,dest in copies.items():shutil.copy2(O/source,D/dest);assert sha(O/source)==sha(D/dest)
 small=D/'egfr_submission_files.zip';public=['egfr_ranked_100.fasta','egfr_track3_top20.csv','egfr_track3_top20.fasta','egfr_100_RESERVE_POOL_NOT_SINGLE_UPLOAD.csv','README.md','make_submission.py']
 with zipfile.ZipFile(small,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for name in public:z.write(O/name,name)
 with zipfile.ZipFile(small) as z:
  assert set(z.namelist())==set(public)
  for name in public:assert hashlib.sha256(z.read(name)).hexdigest()==sha(O/name)
 check_workers()
 verification={'archive':str(target),'archive_bytes':target.stat().st_size,'archive_sha256':sha(target),'payload_files':len(included),'ZIP_members':len(included)+2,'expanded_payload_bytes':sum(r['bytes'] for r in included),'all_payloads_SHA256_and_size_verified':not failures,'failures':failures,'omitted_regenerable_stage6_files':len(excluded),'small_submission_archive':{'path':str(small),'bytes':small.stat().st_size,'sha256':sha(small),'all_members_verified':True},'standalone_files':{dest:{'bytes':(D/dest).stat().st_size,'sha256':sha(D/dest)} for dest in copies.values()},'seconds':time.time()-start,'no_outstanding_scientific_workers':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 (D/'egfr_campaign_final_verification.json').write_text(json.dumps(verification,indent=2));print(json.dumps({k:v for k,v in verification.items() if k!='standalone_files'},indent=2),flush=True)
if __name__=='__main__':main()
