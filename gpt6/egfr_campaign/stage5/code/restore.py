from pathlib import Path
import zipfile, json, hashlib, platform, sys, importlib.metadata
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
manifest=json.loads((R/'stage4/reference/checkpoint_included_manifest.json').read_text())
fail=[]
for r in manifest:
 p=R/r['path']
 if not p.is_file() or p.stat().st_size!=r['bytes'] or sha(p)!=r['sha256']:fail.append(r['path'])
assert not fail,fail
records=[]
for name,dest in [('resources.zip',R/'inputs/provided'),('egfr_inputs.zip',R/'inputs/transfer1'),('egfr_refinement_inputs.zip',R/'inputs/transfer2')]:
 src=Path('/mnt/data')/name
 with zipfile.ZipFile(src) as z:
  for info in z.infolist():
   p=(dest/info.filename).resolve()
   if not p.is_relative_to(dest.resolve()):raise ValueError(info.filename)
  z.extractall(dest)
  records.append(dict(path=str(src),sha256=sha(src),expanded_bytes=sum(i.file_size for i in z.infolist())))
wheel=list((R/'inputs/transfer2').rglob('OpenMM*.whl'))+list((R/'inputs/transfer2').rglob('openmm*.whl'))
assert len(wheel)==1,wheel
runtime=R/'stage2/runtime/python';runtime.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(wheel[0]) as z:z.extractall(runtime)
records.append(dict(path=str(wheel[0]),sha256=sha(wheel[0]),operation='unpack wheel; modern implementation of pre2011 physical methods'))
meta={'verified_stage4_payload_files':len(manifest),'failed_stage4_files':fail,'inputs':records,'python':sys.version,'platform':platform.platform(),'packages':{x:importlib.metadata.version(x) for x in ['numpy','scipy','biopython','numba','rapidfuzz']}}
(S/'reference/restoration.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))
