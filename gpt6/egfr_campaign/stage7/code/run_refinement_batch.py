"""Bounded local jobs, all explicitly tracked. No unattended continuation."""
from pathlib import Path
import sys,json,subprocess,time,os,concurrent.futures
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
def one(task):
 cid,species=task;log=S/'logs'/f'{cid}_{species}_batch.log';cmd=[sys.executable,str(S/'code/refine_one.py'),cid,species];start=time.time();env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1')
 try:
  with log.open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=600)
  out={'candidate_id':cid,'species':species,'returncode':p.returncode,'seconds':time.time()-start,'command':cmd,'log':str(log)}
 except subprocess.TimeoutExpired:out={'candidate_id':cid,'species':species,'returncode':'timeout','seconds':time.time()-start,'command':cmd,'log':str(log)}
 print('FINISHED',json.dumps(out),flush=True);return out
if __name__=='__main__':
 mode=sys.argv[1];species='human6ARU' if mode=='human' else 'mouseAF';ids=sys.argv[2:] or [f'M{i:05d}' for i in range(1,23)];tasks=[(cid,species) for cid in ids];plan={'tasks':tasks,'workers':3,'per_job_timeout_seconds':600,'thread_environment':{'OPENMM_CPU_THREADS':1,'OMP_NUM_THREADS':1,'OPENBLAS_NUM_THREADS':1}}
 (S/'reference'/f'{mode}_batch_plan.json').write_text(json.dumps(plan,indent=2));done=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  for future in concurrent.futures.as_completed([pool.submit(one,t) for t in tasks]):
   done.append(future.result());(S/'reference'/f'{mode}_batch_status.json').write_text(json.dumps({'completed':done,'expected':len(tasks),'all_finished':len(done)==len(tasks)},indent=2))
 print('BATCH_COMPLETE',len(done),flush=True)
