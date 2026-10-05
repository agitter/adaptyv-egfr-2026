"""Verify the complete v2 payload; optionally restore original input archives.
Verification uses only Python's standard library and never runs scientific jobs.
"""
from pathlib import Path
import argparse,json,hashlib,subprocess,sys
R=Path(__file__).resolve().parents[2];S=R/'stage7'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(4*1024*1024),b''):h.update(x)
 return h.hexdigest()
def main():
 a=argparse.ArgumentParser(description=__doc__);a.add_argument('--restore-inputs',action='store_true');a.add_argument('--archive-directory',type=Path,default=Path('/mnt/data'));args=a.parse_args()
 rows=json.loads((S/'reference/checkpoint_v2_included_manifest.json').read_text());fail=[]
 for row in rows:
  p=(R/row['path']).resolve()
  if not p.is_relative_to(R) or not p.is_file() or p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:fail.append(row['path'])
 if fail:raise ValueError('Verification failed: '+repr(fail[:20]))
 print('Verified',len(rows),'v2 payloads')
 if args.restore_inputs:
  if R!=Path('/mnt/data/egfr_campaign'):raise ValueError('Scientific restore scripts use /mnt/data/egfr_campaign; extract there before restoring inputs')
  subprocess.run([sys.executable,str(R/'stage6/code/restore_final.py'),'--restore-inputs','--archive-directory',str(args.archive_directory)],check=True)
if __name__=='__main__':main()
