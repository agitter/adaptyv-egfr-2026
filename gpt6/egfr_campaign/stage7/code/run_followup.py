"""Resource-bounded, fully recorded local follow-up jobs. No external service."""
from pathlib import Path
import sys,os,json,subprocess,time,datetime
S=Path('/mnt/data/egfr_campaign/stage7');R=S.parent
NUMERICAL={'refine_one.py','cluster_one.py','dynamics_one.py','dynamics_unrestrained_one.py','unbound_one.py'}
def count_jobs():
 text=subprocess.check_output(['ps','-eo','pid=,comm=,args='],text=True)
 found=[]
 for line in text.splitlines():
  cells=line.strip().split(None,2)
  if len(cells)<3 or not cells[1].startswith('python') or str(S/'code') not in cells[2]:continue
  if any(('/'+s+' ') in cells[2] or cells[2].endswith('/'+s) for s in NUMERICAL):found.append(int(cells[0]))
 return found

def anon_bytes():
 rows=dict(x.split() for x in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines());return int(rows['anon'])
def write(path,obj):
 t=path.with_suffix('.tmp');t.write_text(json.dumps(obj,indent=2));t.replace(path)

def main():
 tasks=[]
 for cid in ['M00000','M00018','M00021','M00015','M00003','M00008','M00014','M00019','M00020','M00022']:
  tasks.append({'tag':cid+'_mouse','script':'refine_one.py','args':[cid,'mouseAF'],'timeout':700})
 for cid in ['M00000','M00018','M00021']:
  tasks.append({'tag':cid+'_unbound','script':'unbound_one.py','args':[cid],'timeout':600})
 for cid in ['M00018','M00021']:
  tasks.append({'tag':cid+'_cluster','script':'cluster_one.py','args':[cid],'timeout':2200})
 for seed in [198005,198006]:
  for cid in ['M00000','M00018','M00021','NATIVE_3EAK']:
   tasks.append({'tag':f'{cid}_s{seed}_noextra','script':'dynamics_unrestrained_one.py','args':[cid,str(seed)],'timeout':2400})
 ref=S/'reference/followup_plan.json';status=S/'reference/followup_status.json'
 if ref.exists():raise FileExistsError('Refuse to overwrite the followup plan')
 write(ref,{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'tasks':tasks,'max_concurrent_numeric_processes_including_existing':8,'max_anon_before_launch_bytes':int(3.45*1024**3),'individual_threads':1,'precondition':'Human batch finishes first; checkpoint evidence retained.'})
 while True:
  f=S/'reference/human_batch_status.json'
  try:d=json.loads(f.read_text());ready=d.get('all_finished',False)
  except (FileNotFoundError,json.JSONDecodeError):ready=False
  if ready:break
  time.sleep(1)
 active={};done=[];pending=list(tasks)
 while pending or active:
  for tag,item in list(active.items()):
   p=item['process'];elapsed=time.monotonic()-item['start']
   if p.poll() is None and elapsed>item['task']['timeout']:
    p.kill();p.wait();item['timed_out']=True
   if p.poll() is not None:
    item['log_handle'].close();result={**item['task'],'returncode':p.returncode,'timed_out':item.get('timed_out',False),'seconds':elapsed,'command':item['command'],'log':str(item['log'])};done.append(result);print('DONE',json.dumps(result),flush=True);del active[tag]
  while pending and len(count_jobs())<8 and anon_bytes()<3.45*1024**3:
   task=pending.pop(0);cmd=[sys.executable,str(S/'code'/task['script']),*task['args']];log=S/'logs'/('followup_'+task['tag']+'.log');handle=log.open('w');p=subprocess.Popen(cmd,stdout=handle,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1'))
   active[task['tag']]={'process':p,'start':time.monotonic(),'task':task,'log':log,'log_handle':handle,'command':cmd};print('START',task['tag'],p.pid,flush=True);time.sleep(.2)
  write(status,{'completed':done,'running':[{'tag':k,'pid':v['process'].pid,'elapsed_seconds':time.monotonic()-v['start']} for k,v in active.items()],'pending':[x['tag'] for x in pending],'expected':len(tasks),'all_finished':not pending and not active})
  time.sleep(1)
 print('FOLLOWUP_FINISHED',len(done),flush=True)
if __name__=='__main__':main()
