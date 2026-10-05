"""Verify the final checkpoint and optionally restore user-supplied inputs/runtime.
Inspection needs only Python's standard library. Scientific reruns require the
recorded modern numerical packages and the pinned Linux OpenMM implementation.
"""
from pathlib import Path
import argparse,hashlib,json,zipfile
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def extract(src,dest):
 dest=dest.resolve()
 with zipfile.ZipFile(src) as z:
  for i in z.infolist():
   if not (dest/i.filename).resolve().is_relative_to(dest):raise ValueError('Unsafe archive path: '+i.filename)
  z.extractall(dest)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--restore-inputs',action='store_true');p.add_argument('--archive-directory',type=Path,default=Path('/mnt/data'));a=p.parse_args()
 rows=json.loads((S/'reference/final_checkpoint_included_manifest.json').read_text());fail=[]
 for r in rows:
  f=R/r['path']
  if not f.is_file() or f.stat().st_size!=r['bytes'] or sha(f)!=r['sha256']:fail.append(r['path'])
 if fail:raise ValueError('Checkpoint verification failed: '+repr(fail[:20]))
 print('Verified',len(rows),'checkpoint payloads')
 if not a.restore_inputs:return
 meta=json.loads((S/'reference/final_checkpoint_inputs.json').read_text())
 for r in meta['required_original_archives']:
  src=a.archive_directory/r['filename']
  if not src.is_file() or sha(src)!=r['sha256']:raise ValueError('Missing or mismatched original archive: '+str(src))
  extract(src,R/r['destination']);print('Restored',r['filename'])
 wheels=list((R/'inputs/transfer2').rglob('openmm-*.whl'));assert len(wheels)==1
 extract(wheels[0],R/'stage2/runtime/python');print('Restored the supplied OpenMM wheel; no downloaded program was run during extraction.')
if __name__=='__main__':main()
