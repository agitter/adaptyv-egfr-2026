"""Verify no campaign workers, current-source syntax and run completion."""
from pathlib import Path
import os,ast,json,datetime,importlib.metadata,sys,platform
R=Path('/mnt/data/egfr_campaign');S=R/'stage4';workers=[]
for e in Path('/proc').iterdir():
 if not e.name.isdigit() or int(e.name)==os.getpid():continue
 try:args=[a.decode(errors='replace') for a in (e/'cmdline').read_bytes().split(b'\0') if a]
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if not args:continue
 name=Path(args[0]).name
 if (name.startswith('python') and any('/egfr_campaign/' in a and '/code/' in a for a in args[1:])) or (name in ['blastp','makeblastdb'] and any('egfr_campaign' in a for a in args[1:])):workers.append({'pid':int(e.name),'args':args})
assert not workers,workers
syntax=[]
for p in sorted((S/'code').glob('*.py')):ast.parse(p.read_text(),filename=str(p));syntax.append(p.name)
completed={}
for p in (S/'reference').glob('*selection_status.json'):
 d=json.loads(p.read_text());rows=d.get('completed',[]);expected=len(d.get('jobs',[]));assert len(rows)==expected and all(r['returncode']==0 for r in rows),p;completed[p.name]={'expected':expected,'completed_successfully':len(rows)}
f=json.loads((S/'output/final_evidence_inventory.json').read_text());assert f['corrected_trajectory_count']==11
p=S/'intermediate/dynamics/S00022_v2_s198003.json';assert p.exists() and json.loads(p.read_text())['production_ps']==20
versions={}
for name in ['numpy','scipy','biopython','rapidfuzz','numba','openmm']:
 try:versions[name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:versions[name]='not in base metadata; campaign-local runtime recorded separately'
mem=Path('/sys/fs/cgroup/memory.events');memory={a:int(b) for a,b in [l.split() for l in mem.read_text().splitlines()]} if mem.exists() else {}
record={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'campaign_workers_running':workers,'all_stage4_current_python_sources_parse':True,'python_sources_parsed':len(syntax),'source_names':syntax,'completed_task_batches':completed,'final_single_proline_task_completed':True,'memory_events':memory,'python':sys.version,'platform':platform.platform(),'package_metadata':versions,'phone_a_friend_used':2,'checkpoint_only_not_final_submission':True}
(S/'reference/postflight.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
