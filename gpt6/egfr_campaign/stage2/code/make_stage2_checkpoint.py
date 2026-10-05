"""Compact campaign snapshot with provenance for included AND excluded files.

Run only after campaign worker processes have completed or been stopped.
Large immutable inputs and rebuildable runtime/index files are represented by
hashes and retained in the original user-uploaded archives, not duplicated here.
"""
from pathlib import Path
import hashlib,json,zipfile,sys,os,datetime
R=Path('/mnt/data/egfr_campaign'); S=R/'stage2'
DEST=Path('/mnt/data/egfr_stage2_checkpoint.zip')
SELF_GENERATED={'stage2/reference/checkpoint_included_manifest.json','stage2/reference/checkpoint_excluded_manifest.json','stage2/reference/checkpoint_archive_inputs.json','stage2/reference/checkpoint_summary_packaging.json'}

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def include(p):
    r=p.relative_to(R);s=str(r)
    if '__pycache__' in r.parts or p.suffix in ['.pyc','.nbi','.nbc']:return False,'regenerable cache'
    if s.startswith('stage2/runtime/'):return False,'runtime restored from verified original wheel/archive'
    if s.startswith('stage2/intermediate/blast/'):return False,'regenerable BLAST index'
    if s in ['stage2/intermediate/swissprot.fasta','stage2/intermediate/pdb.fasta']:return False,'decompressed original sequence input'
    if s.startswith('stage2/intermediate/novelty/') and p.name in ['antibodies.fasta','antibodies_augmented.fasta','antibody_sources.jsonl','cdrh3_candidates.jsonl','cdrh3_augmented.jsonl','cdrh3_sources.jsonl']:
        return False,'reference normalization output; regenerable from original antibody inputs'
    if s.startswith('stage2/intermediate/novelty/') and p.stat().st_size>5*1024*1024:return False,'large regenerable normalized reference data'
    if s.startswith('inputs/'):
        if p.name.endswith('.provenance.json') or p.name in ['download_manifest.json','manifest.json','README.txt']:return True,'input provenance metadata'
        return False,'original input retained in uploaded resources/transfer archives'
    return True,'campaign source or small intermediate/result'

def main():
    inc=[];exc=[]
    for p in sorted(R.rglob('*')):
        if not p.is_file() or str(p.relative_to(R)) in SELF_GENERATED:continue
        yes,reason=include(p)
        row={'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p),'reason':reason}
        (inc if yes else exc).append(row)
    archives=[]
    for name in ['resources.zip','egfr_inputs.zip','egfr_refinement_inputs.zip','egfr_campaign_checkpoint.zip','egfr_request2_audit.zip']:
        p=Path('/mnt/data')/name
        if p.exists():archives.append({'name':name,'bytes':p.stat().st_size,'sha256':sha(p),'role':'original input/recovery archive'})
    for name,data in [('checkpoint_included_manifest.json',inc),('checkpoint_excluded_manifest.json',exc),('checkpoint_archive_inputs.json',archives)]:
        (S/'reference'/name).write_text(json.dumps(data,indent=2))
    with zipfile.ZipFile(DEST,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for row in inc:
            p=R/row['path'];z.write(p,arcname='egfr_campaign/'+row['path'])
        for name in ['checkpoint_included_manifest.json','checkpoint_excluded_manifest.json','checkpoint_archive_inputs.json']:
            p=S/'reference'/name;z.write(p,arcname='egfr_campaign/'+str(p.relative_to(R)))
    with zipfile.ZipFile(DEST) as z:
        bad=z.testzip();assert bad is None,bad
        archived={a.filename for a in z.infolist()};assert len(archived)==len(inc)+3
    result={'archive':str(DEST),'bytes':DEST.stat().st_size,'MiB':round(DEST.stat().st_size/1024**2,2),
      'sha256':sha(DEST),'included_files':len(inc)+3,'manifested_excluded_files':len(exc),
      'included_uncompressed_bytes':sum(a['bytes'] for a in inc),'zip_integrity_test':'pass',
      'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'missing_prior_stage1_archive_disclosed':True,'not_final_submission':True}
    (S/'reference/checkpoint_summary_packaging.json').write_text(json.dumps(result,indent=2))
    Path('/mnt/data/egfr_stage2_checkpoint.sha256').write_text(result['sha256']+'  '+DEST.name+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
