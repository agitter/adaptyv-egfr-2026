from pathlib import Path
import zipfile, hashlib, json, sys, platform, importlib.metadata
ROOT=Path('/mnt/data');R=ROOT/'egfr_campaign';S=R/'stage6'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 return h.hexdigest()
def extract(src,dest):
 dest.mkdir(exist_ok=True,parents=True)
 with zipfile.ZipFile(src) as z:
  for info in z.infolist():
   if not (dest/info.filename).resolve().is_relative_to(dest.resolve()):raise ValueError(info.filename)
  z.extractall(dest)
  return {'path':str(src),'sha256':sha(src),'entries':len(z.infolist()),'expanded_bytes':sum(i.file_size for i in z.infolist())}
archives=[extract(ROOT/'egfr_stage5_checkpoint.zip',ROOT)]
for d in ['code','logs','reference','intermediate','output']:(S/d).mkdir(parents=True,exist_ok=True)
manifest=json.loads((R/'stage5/reference/checkpoint_included_manifest.json').read_text());failed=[]
for row in manifest:
 p=R/row['path']
 if not p.is_file() or p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:failed.append(row['path'])
assert not failed,failed
for name,dest in [('resources.zip',R/'inputs/provided'),('egfr_inputs.zip',R/'inputs/transfer1'),('egfr_refinement_inputs.zip',R/'inputs/transfer2')]:archives.append(extract(ROOT/name,dest))
wheels=list((R/'inputs/transfer2').rglob('openmm*.whl'));assert len(wheels)==1
archives.append(extract(wheels[0],R/'stage2/runtime/python'))
packages={}
for p in ['numpy','scipy','biopython','numba','rapidfuzz']:
 try:packages[p]=importlib.metadata.version(p)
 except importlib.metadata.PackageNotFoundError:packages[p]=None
record={'verified_stage5_payloads':len(manifest),'failures':failed,'archives':archives,'python':sys.version,'platform':platform.platform(),'packages':packages}
(S/'reference/restoration.json').write_text(json.dumps(record,indent=2));(S/'code/restore_stage6.py').write_text(Path(__file__).read_text());print(json.dumps(record,indent=2))
