"""Package all preserved scientific files, with hashes for exclusions too.
Run only after postflight; original inputs, runtimes and rebuilt databases remain
recoverable from their original uploaded archives. No fonts are distributed.
"""
from pathlib import Path
import hashlib,json,zipfile,datetime,os
R=Path('/mnt/data/egfr_campaign');S=R/'stage5';DEST=Path('/mnt/data/egfr_stage5_checkpoint.zip')
NAMES=['checkpoint_included_manifest.json','checkpoint_excluded_manifest.json','checkpoint_archive_inputs.json','checkpoint_summary_packaging.json','checkpoint_independent_verification.json']
SELF={'stage5/reference/'+x for x in NAMES}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def include(p):
 r=p.relative_to(R);s=str(r)
 if '__pycache__' in r.parts or p.suffix in ['.pyc','.nbi','.nbc']:return False,'regenerable cache'
 if p.suffix.lower() in ['.ttf','.otf','.woff','.woff2']:return False,'font file not distributed'
 if any(s.startswith(st+'/runtime/') for st in ['stage2','stage3','stage4','stage5']):return False,'runtime restored from original verified archives/wheels'
 if s.startswith('stage2/intermediate/blast/'):return False,'regenerable prior BLAST index'
 if s in ['stage2/intermediate/swissprot.fasta','stage2/intermediate/pdb.fasta']:return False,'decompressed original sequence input'
 if s.startswith('stage2/intermediate/novelty/') and (p.name in ['antibodies.fasta','antibodies_augmented.fasta','antibody_sources.jsonl','cdrh3_candidates.jsonl','cdrh3_augmented.jsonl','cdrh3_sources.jsonl'] or p.stat().st_size>5*1024*1024):return False,'regenerable reference normalization from supplied inputs'
 if any(s.startswith(st+'/intermediate/blast/') for st in ['stage4','stage5']) and (p.suffix in ['.phr','.pin','.psq','.pjs','.pog','.psd','.psi','.pot','.pto','.ptf'] or p.name in ['swissprot.fasta','pdb.fasta','antibodies_augmented.fasta','antibody_sources.jsonl']):return False,'regenerable reference database; scientific query/result files retained'
 if s in ['stage4/intermediate/novelty/cdrh3_rebuilt.jsonl','stage5/intermediate/novelty/cdrh3_rebuilt.jsonl']:return False,'regenerable CDR3 reference; all query nearest-neighbor results retained'
 if s.startswith('inputs/'):
  if p.name.endswith('.provenance.json') or p.name in ['download_manifest.json','manifest.json','README.txt']:return True,'input provenance metadata'
  return False,'immutable original input in uploaded archives'
 return True,'campaign code, scientific intermediate, control, log or result'
def main():
 post=json.loads((S/'reference/postflight.json').read_text());assert not post['campaign_workers_running']
 workers=[]
 for e in Path('/proc').iterdir():
  if not e.name.isdigit() or int(e.name)==os.getpid():continue
  try:args=[a.decode(errors='replace') for a in (e/'cmdline').read_bytes().split(b'\0') if a]
  except (FileNotFoundError,PermissionError,ProcessLookupError):continue
  if not args:continue
  name=Path(args[0]).name
  if (name.startswith('python') and any('/egfr_campaign/' in a and '/code/' in a for a in args[1:])) or (name in ['blastp','makeblastdb'] and any('egfr_campaign' in a for a in args[1:])):workers.append({'pid':e.name,'args':args})
 assert not workers,workers
 inc=[];exc=[]
 for p in sorted(R.rglob('*')):
  if not p.is_file() or str(p.relative_to(R)) in SELF:continue
  yes,why=include(p);row={'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p),'reason':why};(inc if yes else exc).append(row)
 archives=[]
 for name in ['resources.zip','egfr_inputs.zip','egfr_refinement_inputs.zip','egfr_campaign_checkpoint.zip','egfr_request2_audit.zip','egfr_stage2_checkpoint.zip','egfr_stage3_checkpoint.zip','egfr_stage4_checkpoint.zip']:
  p=Path('/mnt/data')/name
  if p.exists():archives.append({'name':name,'bytes':p.stat().st_size,'sha256':sha(p),'role':'original input or restoration archive'})
 for name,data in zip(NAMES[:3],[inc,exc,archives]):(S/'reference'/name).write_text(json.dumps(data,indent=2))
 with zipfile.ZipFile(DEST,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
  for row in inc:z.write(R/row['path'],arcname='egfr_campaign/'+row['path'])
  for name in NAMES[:3]:p=S/'reference'/name;z.write(p,arcname='egfr_campaign/'+str(p.relative_to(R)))
 with zipfile.ZipFile(DEST) as z:
  bad=z.testzip();assert bad is None,bad;assert len(set(z.namelist()))==len(inc)+3
 result={'archive':str(DEST),'bytes':DEST.stat().st_size,'MiB':round(DEST.stat().st_size/1024**2,2),'sha256':sha(DEST),'included_files':len(inc)+3,'manifested_scientific_payload_files':len(inc),'manifested_excluded_files':len(exc),'included_uncompressed_bytes':sum(a['bytes'] for a in inc),'zip_integrity_test':'pass','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'not_final_submission':True,'no_running_campaign_workers':True}
 (S/'reference/checkpoint_summary_packaging.json').write_text(json.dumps(result,indent=2));Path('/mnt/data/egfr_stage5_checkpoint.sha256').write_text(result['sha256']+'  '+DEST.name+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
