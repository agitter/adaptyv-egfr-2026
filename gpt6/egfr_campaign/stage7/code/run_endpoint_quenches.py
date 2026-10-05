"""Join twelve endpoint quenches while respecting eight numerical CPU workers.
The queue exits after all inputs from the ongoing, in-session campaign finish.
"""
from pathlib import Path
import sys,os,subprocess,time,json,shlex
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1')

def numerical_workers():
 count=0
 for line in subprocess.check_output(['ps','-eo','args'],text=True).splitlines():
  try:args=shlex.split(line)
  except ValueError:continue
  if not args or not Path(args[0]).name.startswith('python'):continue
  scripts=[x for x in args[1:] if x.startswith(str(R)) and '/code/' in x]
  if not scripts:continue
  name=Path(scripts[0]).name
  if name in {'dynamics_one.py','dynamics_unrestrained_one.py','cluster_one.py','refine_one.py','unbound_one.py','quench_endpoint.py'}:count+=1
 return count

def main():
 plan=[]
 for extra,cids in [(True,['M00000','NATIVE_3EAK']),(False,['M00000','M00018','M00021','NATIVE_3EAK'])]:
  for cid in cids:
   for seed in [198005,198006]:
    sub=S/'intermediate/dynamics' if extra else S/'intermediate/thermal_unrestrained'/f'{cid}_s{seed}'/'intermediate/dynamics'
    source=sub/f'{cid}_v2_s{seed}.json';tag=f'{cid}_s{seed}_'+('extra' if extra else 'noextra')
    plan.append({'tag':tag,'source':str(source),'candidate_id':cid,'seed':seed,'added_restraints':extra})
 (S/'reference/endpoint_quench_plan.json').write_text(json.dumps(plan,indent=2));remaining=list(plan);active=[];done=[];begin=time.monotonic()
 while remaining or active:
  for rec in active[:]:
   job,p,start,log=rec
   if p.poll() is None and time.monotonic()-start<600:continue
   timed=p.poll() is None
   if timed:p.kill()
   code=p.wait();log.close();row={**job,'returncode':code,'timed_out':timed,'seconds':time.monotonic()-start};done.append(row);active.remove(rec);print('DONE',row,flush=True)
  ready=[j for j in remaining if Path(j['source']).is_file()]
  while ready and numerical_workers()<8:
   job=ready.pop(0);remaining.remove(job);log=(S/'logs'/('quench_'+job['tag']+'.log')).open('w');cmd=[sys.executable,str(S/'code/quench_endpoint.py'),job['source'],job['tag']];p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env);active.append((job,p,time.monotonic(),log));print('START',job['tag'],p.pid,flush=True)
  (S/'reference/endpoint_quench_status.json').write_text(json.dumps({'expected':len(plan),'completed':done,'active':[j['tag'] for j,p,t,l in active],'waiting':[j['tag'] for j in remaining],'all_finished':not remaining and not active},indent=2))
  if time.monotonic()-begin>4000:raise TimeoutError('Endpoint queue exceeded bounded session run')
  if remaining or active:time.sleep(3)
 assert len(done)==12 and all(j['returncode']==0 and not j['timed_out'] for j in done)
 print('ALL_ENDPOINT_QUENCHES_FINISHED',flush=True)
if __name__=='__main__':main()
