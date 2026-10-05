"""Package the current campaign, excluding reproducible large data and fonts."""
from pathlib import Path
import hashlib,json,time,zipfile
R=Path('/mnt/data/egfr_campaign')
Z=Path('/mnt/data/egfr_campaign_checkpoint.zip')
font_ext={'.ttf','.otf','.woff','.woff2','.eot','.pfb','.pfa','.ttc'}
rows=[]
for p in sorted(R.rglob('*')):
 if not p.is_file() or p.name in {'checkpoint_manifest.json','checkpoint_summary.json'}:continue
 rel=p.relative_to(R);parts=rel.parts
 include=True;reason='included'
 if p.suffix.lower() in font_ext:include=False;reason='font file excluded'
 elif '__pycache__' in parts or p.suffix in {'.pyc','.nbc','.nbi'}:include=False;reason='regenerable compiled cache'
 elif parts[0] in {'databases','tools'}:include=False;reason='rebuild from transferred software/data'
 elif p.name.endswith(('.tar.gz','.fasta.gz')):include=False;reason='retained in original egfr_inputs.zip'
 elif p.name=='pdb_seqres.txt.gz':include=False;reason='retained in original egfr_inputs.zip'
 elif p.name in {'swissprot.fasta','pdb.fasta'}:include=False;reason='decompress from transferred dataset'
 elif len(parts)>2 and parts[0]=='inputs' and parts[1]=='provided' and p.suffix not in {'.pdf','.cif','.md','.html'}:
  include=False;reason='saved web-page asset reproducible from resources.zip'
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 rows.append(dict(path=rel.as_posix(),bytes=p.stat().st_size,sha256=h.hexdigest(),included=include,reason=reason))
manifest=R/'reference/checkpoint_manifest.json'
manifest.write_text(json.dumps(dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),files=rows),indent=2))
with zipfile.ZipFile(Z,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as out:
 for rec in rows:
  if rec['included']:out.write(R/rec['path'],'egfr_campaign/'+rec['path'])
 out.write(manifest,'egfr_campaign/reference/checkpoint_manifest.json')
with zipfile.ZipFile(Z) as z:
 assert z.testzip() is None
 assert not any(Path(n).suffix.lower() in font_ext for n in z.namelist())
summary=dict(path=str(Z),bytes=Z.stat().st_size,included_files=sum(x['included'] for x in rows),
 excluded_files=sum(not x['included'] for x in rows),sha256=hashlib.sha256(Z.read_bytes()).hexdigest())
(R/'reference/checkpoint_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
