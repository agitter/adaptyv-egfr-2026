"""Finite local batch; records every worker exit."""
from pathlib import Path
import subprocess,os,sys,time,json
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
config=json.loads(Path(sys.argv[1]).read_text());parallel=int(sys.argv[2]) if len(sys.argv)>2 else 3
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1')
active=[];results=[];todo=list(config);start=time.time()
while todo or active:
 while todo and len(active)<parallel:
  r=todo.pop(0);log=S/'logs'/r['log'];f=log.open('w');p=subprocess.Popen([sys.executable,*r['args']],env=env,stdout=f,stderr=subprocess.STDOUT);active.append((p,f,r,time.time()));print('START',r['log'],p.pid,flush=True)
 for a in list(active):
  p,f,r,t=a;c=p.poll()
  if c is not None:
   f.close();active.remove(a);results.append({**r,'returncode':c,'seconds':time.time()-t});print('END',r['log'],c,round(time.time()-t,1),flush=True)
 if active:time.sleep(1)
(S/'reference'/(Path(sys.argv[1]).stem+'_execution.json')).write_text(json.dumps({'tasks':results,'seconds':time.time()-start,'all_finished':True},indent=2));sys.exit(int(any(r['returncode'] for r in results)))
