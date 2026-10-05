"""Bounded local batch runner. Every subprocess is joined before completion.
Status and logs are retained; no successful result is inferred from launch alone.
"""
from pathlib import Path
import concurrent.futures,subprocess,os,json,time,sys
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
def main():
 plan_path=Path(sys.argv[1]);plan=json.loads(plan_path.read_text());status=S/'reference'/(plan['name']+'_execution.json');start=time.time();env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1');records=[]
 def run(job):
  t=time.time();log=S/'logs'/(plan['name']+'_'+job['id']+'.log');cmd=job['command']
  try:
   with log.open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env,timeout=job.get('timeout',600))
   return {'id':job['id'],'command':cmd,'returncode':p.returncode,'seconds':time.time()-t,'log':str(log)}
  except subprocess.TimeoutExpired:return {'id':job['id'],'command':cmd,'returncode':124,'seconds':time.time()-t,'log':str(log),'error':'timeout; subprocess terminated'}
 def save(done=False):status.write_text(json.dumps({'plan':str(plan_path),'jobs':len(plan['jobs']),'completed':len(records),'finished':done,'seconds':time.time()-start,'records':records},indent=2))
 save()
 with concurrent.futures.ThreadPoolExecutor(max_workers=plan.get('concurrency',3)) as pool:
  futures=[pool.submit(run,j) for j in plan['jobs']]
  for f in concurrent.futures.as_completed(futures):
   r=f.result();records.append(r);save();print(r['id'],r['returncode'],round(r['seconds'],1),flush=True)
 save(True);print('BATCH COMPLETE',plan['name'],len(records),flush=True)
if __name__=='__main__':main()
