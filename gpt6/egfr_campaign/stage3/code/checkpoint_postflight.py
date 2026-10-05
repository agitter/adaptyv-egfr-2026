"""Record worker completion, memory counters, and Python syntax checks."""
from pathlib import Path
import os,ast,json,datetime
R=Path('/mnt/data/egfr_campaign');S=R/'stage3';workers=[]
for entry in Path('/proc').iterdir():
    if not entry.name.isdigit() or int(entry.name)==os.getpid():continue
    try:args=[a.decode(errors='replace') for a in (entry/'cmdline').read_bytes().split(b'\0') if a]
    except (FileNotFoundError,PermissionError,ProcessLookupError):continue
    if not args:continue
    name=Path(args[0]).name
    if (name.startswith('python') and any('/egfr_campaign/stage3/code/' in a for a in args[1:])) or (name in ['blastp','makeblastdb'] and any('egfr_campaign' in a for a in args[1:])):workers.append({'pid':int(entry.name),'args':args})
if workers:raise RuntimeError('Active scientific workers: '+repr(workers))
syntax=[]
for p in sorted((S/'code').glob('*.py')):ast.parse(p.read_text(),filename=str(p));syntax.append(p.name)
mem=Path('/sys/fs/cgroup/memory.events');memory={a:int(b) for a,b in [x.split() for x in mem.read_text().splitlines()]} if mem.exists() else {}
record={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'campaign_workers_running':workers,'all_stage3_python_sources_parse':True,'python_sources_parsed':len(syntax),'source_names':syntax,'memory_events':memory,'phone_a_friend_used':2,'checkpoint_only_not_final_submission':True}
(S/'reference/postflight.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
