"""Recover from supplied inputs; do not treat missing earlier outputs as verified."""
from pathlib import Path
import json,hashlib,sys,platform,zipfile,tarfile,subprocess,time,shutil,gzip
R=Path('/mnt/data/egfr_campaign');S=R/'stage2';T=R/'inputs/transfer2/egfr_refinement_fetch'
rep={'date':'2026-10-01','python':sys.version,'platform':platform.platform(),'missing_checkpoint':'egfr_stage1_checkpoint.zip not present in current runtime','recovered_archives':['egfr_campaign_checkpoint.zip','egfr_request2_audit.zip','egfr_inputs.zip','egfr_refinement_inputs.zip'],'checks':[]}
for p in sorted(T.glob('*.provenance.json')):
 d=json.loads(p.read_text());fn=p.with_name(p.name.replace('.provenance.json',''));exp=d.get('sha256');actual=hashlib.sha256(fn.read_bytes()).hexdigest() if fn.exists() else None
 rep['checks'].append({'file':fn.name,'expected':exp,'actual':actual,'pass':exp==actual})
if not all(x['pass'] for x in rep['checks']):print('Some provenance schemas need inspection',rep['checks'])
wheel=next(T.glob('openmm-*.whl')); dest=S/'runtime/python'; dest.mkdir(exist_ok=True)
with zipfile.ZipFile(wheel) as z:z.extractall(dest)
print('OpenMM extracted')
arc=R/'inputs/transfer1/legacy_fetch/ncbi-blast-2.17.0+-x64-linux.tar.gz'
with tarfile.open(arc) as t:
 for m in t.getmembers():
  if m.isfile() and m.name.split('/')[-1] in ['blastp','makeblastdb','blastdbcmd']:
   m.name=m.name.split('/')[-1];t.extract(m,path=S/'runtime/bin',filter='data')
for src,name in [('uniprot_sprot.fasta.gz','swissprot'),('pdb_seqres.txt.gz','pdb')]:
 out=S/'intermediate'/f'{name}.fasta'
 with gzip.open(R/'inputs/transfer1/legacy_fetch'/src,'rb') as f,open(out,'wb') as w:shutil.copyfileobj(f,w)
 print('Restored',name,out.stat().st_size,flush=True)
rep['numpy_version']=__import__('numpy').__version__
(S/'reference/recovery.json').write_text(json.dumps(rep,indent=2))
print('Recovery completed')
