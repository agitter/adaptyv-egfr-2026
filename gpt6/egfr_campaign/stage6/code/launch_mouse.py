from pathlib import Path
import json,subprocess
S=Path('/mnt/data/egfr_campaign/stage6');files=sorted((S/'intermediate/designs').glob('J*.json'));plan={'name':'mouse_factorial','concurrency':2,'selection':'All predeclared common-backbone sequences, including unchanged controls; no human-score cherry-picking.','jobs':[{'id':p.stem,'command':['python',str(S/'code/refine_one.py'),p.stem,'mouseAF'],'timeout':600} for p in files]};path=S/'reference/mouse_factorial.json';path.write_text(json.dumps(plan,indent=2))
with (S/'logs/mouse_factorial_manager.log').open('w') as f:p=subprocess.Popen(['python',str(S/'code/run_batch.py'),str(path)],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
(S/'reference/mouse_factorial_process.json').write_text(json.dumps({'pid':p.pid,'plan':str(path)}));print('Mouse followup started',len(files),p.pid)
