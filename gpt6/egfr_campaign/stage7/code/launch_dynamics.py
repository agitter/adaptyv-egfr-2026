from pathlib import Path
import sys,json,subprocess,os,concurrent.futures,time
S=Path('/mnt/data/egfr_campaign/stage7')
def run(t):
 cid,seed=t;cmd=[sys.executable,str(S/'code/dynamics_one.py'),cid,str(seed)];f=S/'logs'/f'dynamics_{cid}_{seed}.log';tm=time.time()
 try:
  with f.open('w') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=1800,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1'))
  rc=p.returncode
 except subprocess.TimeoutExpired:rc='timeout'
 d={'id':cid,'seed':seed,'returncode':rc,'seconds':time.time()-tm,'command':cmd,'log':str(f)};print(json.dumps(d),flush=True);return d
if __name__=='__main__':
 group=sys.argv[1];ids=sys.argv[2:];tasks=[(cid,seed) for cid in ids for seed in [198005,198006]];status=S/'reference'/f'dynamics_{group}_execution.json';(S/'reference'/f'dynamics_{group}_plan.json').write_text(json.dumps({'tasks':tasks,'workers':4,'threads_each':1,'production_ps':20,'warmup_ps':2,'source':'preserved stage4 dynamics.py','definition':'Replica seeds fixed before observing trajectories. Report both complete trajectories, not selected early segments. RMSD is not folding free energy.'},indent=2));done=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for f in concurrent.futures.as_completed([ex.submit(run,t) for t in tasks]):done.append(f.result());status.write_text(json.dumps({'completed':done,'expected':len(tasks),'all_finished':len(done)==len(tasks)},indent=2))
 print('DYNAMICS_GROUP_COMPLETE',group,flush=True)
