"""Bounded subprocess execution with durable status, no deferred deliverables."""
import json,sys,subprocess,os,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
S=Path('/mnt/data/egfr_campaign/stage4');config=Path(sys.argv[1]);jobs=json.loads(config.read_text());workers=int(sys.argv[2]) if len(sys.argv)>2 else 3;state={'config':str(config),'jobs':jobs,'completed':[]};dest=S/'reference'/(config.stem+'_status.json')
def worker(job):
    start=time.time();log=S/'logs'/(job['name']+'.log')
    with log.open('w') as f:
        try:proc=subprocess.run([sys.executable,*job['args']],stdout=f,stderr=subprocess.STDOUT,env=os.environ.copy(),timeout=job.get('timeout',1500));status=proc.returncode
        except subprocess.TimeoutExpired:status='timeout'
    return {'name':job['name'],'returncode':status,'seconds':time.time()-start,'log':str(log)}
with ThreadPoolExecutor(max_workers=workers) as pool:
    futures=[pool.submit(worker,j) for j in jobs]
    for f in as_completed(futures):
        r=f.result();state['completed'].append(r);dest.write_text(json.dumps(state,indent=2));print(r,flush=True)
print('ALL DONE',len(jobs),flush=True)
