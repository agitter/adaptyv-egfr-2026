import zipfile,pathlib,hashlib,json,datetime,os,tarfile,subprocess,sys
root=pathlib.Path('/mnt/data/egfr_campaign')
records=[]
for fn,sub in [('resources.zip','provided'),('egfr_inputs.zip','transfer1')]:
 src=root.parent/fn; dst=root/'inputs'/sub; dst.mkdir(exist_ok=True)
 with zipfile.ZipFile(src) as z:
  expanded=sum(i.file_size for i in z.infolist())
  if expanded>4_000_000_000: raise ValueError('Unexpected expanded size')
  for i in z.infolist():
   target=(dst/i.filename).resolve()
   if not target.is_relative_to(dst.resolve()): raise ValueError('Unsafe member')
  z.extractall(dst)
 records.append({'file':str(src),'sha256':hashlib.file_digest(open(src,'rb'),'sha256').hexdigest(),'expanded_bytes':expanded})
for p in sorted((root/'inputs').rglob('*')):
 if p.is_file():records.append({'file':str(p.relative_to(root)),'size':p.stat().st_size,'sha256':hashlib.file_digest(open(p,'rb'),'sha256').hexdigest()})
(root/'reference'/'input_manifest.json').write_text(json.dumps(records,indent=2))
p=root/'inputs'/'transfer1'/'legacy_fetch'
manifest=json.loads((p/'manifest_v2.json').read_text())
print('MANIFEST TYPE:', type(manifest).__name__)
print(json.dumps(manifest,indent=2)[:3200])
for fn in ['1IVO.pdb','1NQL.pdb','1YY9.pdb','2A3D.pdb','1FSD.pdb','1L2Y.pdb','3DWT.pdb','3EAK.pdb','3EBA.pdb']:
 text=(p/fn).read_text()
 print('\n'+fn+'\n'+'\n'.join(line for line in text.splitlines() if line.startswith(('HEADER','TITLE','COMPND','JRNL        AUTH','JRNL        TITL','JRNL        REF ')))[:2300])
print('Python:',sys.version)
print('CPU:',os.cpu_count())
print(open('/proc/meminfo').read()[:140])
