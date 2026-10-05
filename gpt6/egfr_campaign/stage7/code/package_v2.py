"""Archive the prior scientific record plus v2 and verify every payload hash."""
from pathlib import Path
import json,hashlib,zipfile,os,shlex,subprocess,shutil,datetime,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';O=S/'output';D=Path('/mnt/data')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(4*1024*1024),b''):h.update(x)
 return h.hexdigest()
def workers():
 bad=[]
 for line in subprocess.check_output(['ps','-eo','pid,args'],text=True).splitlines()[1:]:
  try:p,cmd=line.strip().split(None,1);pid=int(p);a=shlex.split(cmd)
  except (ValueError,IndexError):continue
  if pid==os.getpid() or not a:continue
  if Path(a[0]).name.startswith('python') and any(x.startswith(str(R)) and '/code/' in x for x in a[1:]):bad.append([pid,a])
  elif str(R) in a[0] and Path(a[0]).name in ['blastp','makeblastdb']:bad.append([pid,a])
 return bad

def main():
 start=time.time();assert not workers(),workers()
 post=json.loads((S/'reference/release_validation_v2.json').read_text());assert post['all_pass']
 for n in ['followup_status.json','endpoint_quench_status.json']:
  d=json.loads((S/'reference'/n).read_text());assert d['all_finished'] and all(r['returncode']==0 and not r['timed_out'] for r in d['completed'])
 # These historic payloads must remain unchanged, even when earlier conclusions differ.
 old=json.loads((R/'stage6/reference/final_checkpoint_included_manifest.json').read_text())
 for row in old:
  p=R/row['path'];assert p.is_file() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],row['path']
 (S/'reference/old_payloads_verified_at_v2_release.json').write_text(json.dumps({'payloads':len(old),'all_size_and_SHA256_match':True,'historical_state_not_overwritten':True},indent=2))
 (S/'reference/no_active_campaign_workers.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'active_scientific_workers':workers(),'two_phone_a_friend_requests_used':True,'one_unused':True},indent=2))
 base=D/'egfr_campaign_final_checkpoint.zip'
 with zipfile.ZipFile(base) as z:
  paths={i.filename[len('egfr_campaign/'):] for i in z.infolist() if not i.is_dir() and i.filename.startswith('egfr_campaign/')}
 paths.add('CURRENT_RELEASE.md');omitted=[];special={'stage7/reference/checkpoint_v2_included_manifest.json','stage7/reference/checkpoint_v2_omitted_manifest.json'}
 for p in S.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(R).as_posix();sub=p.relative_to(S).as_posix();reason=None
  if rel in special:continue
  if '__pycache__' in p.parts or p.suffix in {'.pyc','.nbc','.nbi'}:reason='Regenerable interpreter/JIT cache'
  elif sub.startswith('runtime/'):reason='Regenerable software from the original supplied archive'
  elif sub.startswith('intermediate/blast/') and (p.name.startswith(('swissprot.','pdb.','antibodies_augmented.')) or p.name=='antibody_sources.jsonl'):reason='Regenerable sequence database, normalization or BLAST index; actual query and result files are included'
  elif sub=='intermediate/novelty/cdrh3_rebuilt.jsonl':reason='Regenerable CDRH3 reference normalization; screening result, source hash and code are included'
  if reason:omitted.append({'path':rel,'bytes':p.stat().st_size,'sha256':sha(p),'reason':reason})
  else:paths.add(rel)
 included=[]
 for rel in sorted(paths):
  p=R/rel;assert p.is_file(),rel;included.append({'path':rel,'bytes':p.stat().st_size,'sha256':sha(p)})
 im=S/'reference/checkpoint_v2_included_manifest.json';om=S/'reference/checkpoint_v2_omitted_manifest.json';im.write_text(json.dumps(included,indent=2));om.write_text(json.dumps(omitted,indent=2))
 target=D/'egfr_campaign_v2_checkpoint.zip';temp=target.with_suffix('.zip.tmp')
 with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for row in included:z.write(R/row['path'],'egfr_campaign/'+row['path'])
  for p in [im,om]:z.write(p,'egfr_campaign/'+p.relative_to(R).as_posix())
 os.replace(temp,target)
 print('ARCHIVE_WRITTEN',target.stat().st_size,'bytes',len(included),'manifested payloads',flush=True)
 # Independent decompression-based verification, not merely a CRC or filesystem check.
 failures=[]
 with zipfile.ZipFile(target) as z:
  assert len(z.namelist())==len(set(z.namelist()))
  expect={'egfr_campaign/'+r['path'] for r in included}|{'egfr_campaign/'+p.relative_to(R).as_posix() for p in [im,om]}
  assert set(z.namelist())==expect
  stored=json.loads(z.read('egfr_campaign/stage7/reference/checkpoint_v2_included_manifest.json'));assert stored==included
  for row in stored:
   n=0;h=hashlib.sha256()
   with z.open('egfr_campaign/'+row['path']) as f:
    for x in iter(lambda:f.read(1024*1024),b''):h.update(x);n+=len(x)
   if n!=row['bytes'] or h.hexdigest()!=row['sha256']:failures.append(row['path'])
 assert not failures,failures
 copies={'egfr_ranked_100_v2.fasta':'egfr_ranked_100_v2.fasta','egfr_track3_top20_v2.csv':'egfr_track3_top20_v2.csv','egfr_track3_top20_v2.fasta':'egfr_track3_top20_v2.fasta','REVISION_REPORT.md':'egfr_revision_v2_report.md','METHODS_CERTIFICATE_V2.md':'egfr_methods_certificate_v2.md','private_codebook_v2.tsv':'egfr_private_codebook_v2.tsv','revision_summary.json':'egfr_revision_v2_summary.json'}
 for a,b in copies.items():shutil.copy2(O/a,D/b);assert sha(O/a)==sha(D/b)
 shutil.copy2(S/'reference/release_validation_v2.json',D/'egfr_release_validation_v2.json')
 small=D/'egfr_submission_files_v2.zip';public=['egfr_ranked_100_v2.fasta','egfr_track3_top20_v2.csv','egfr_track3_top20_v2.fasta','egfr_100_v2_RESERVES_NOT_SINGLE_UPLOAD.csv','README.md','make_submission_v2.py']
 with zipfile.ZipFile(small,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for name in public:z.write(O/name,name)
 with zipfile.ZipFile(small) as z:
  assert set(z.namelist())==set(public)
  for name in public:assert hashlib.sha256(z.read(name)).hexdigest()==sha(O/name)
 assert not workers(),workers()
 v={'archive':str(target),'archive_bytes':target.stat().st_size,'archive_sha256':sha(target),'manifested_payload_count':len(included),'zip_member_count':len(included)+2,'expanded_payload_bytes':sum(r['bytes'] for r in included),'all_payload_sizes_and_SHA256_verified':True,'failures':failures,'omitted_regenerable_files':len(omitted),'unchanged_prior_payloads':len(old),'submission_zip':{'path':str(small),'bytes':small.stat().st_size,'sha256':sha(small),'all_members_verified':True},'standalone_files':{n:{'bytes':(D/n).stat().st_size,'sha256':sha(D/n)} for n in copies.values()},'no_campaign_workers_running':True,'seconds':time.time()-start,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 (D/'egfr_campaign_v2_verification.json').write_text(json.dumps(v,indent=2));print(json.dumps({k:x for k,x in v.items() if k!='standalone_files'},indent=2),flush=True)
if __name__=='__main__':main()
