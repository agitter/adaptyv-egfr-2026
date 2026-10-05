"""Run a predeclared matched longer stress test; wait only for local control."""
from pathlib import Path
import sys,time
S=Path('/mnt/data/egfr_campaign/stage4')
sys.path.insert(0,str(S/'code'))
from dynamics import run
if __name__=='__main__':
 cid=sys.argv[1]
 if cid.startswith('MATCHED_'):
  cid0=cid.removeprefix('MATCHED_');needed=S/'intermediate/matched_controls/intermediate/refined'/(cid0+'_independent_v2.json');start=time.time()
  while not needed.exists():
   if time.time()-start>600:raise TimeoutError('Matched-parent refinement did not complete')
   time.sleep(2)
 run(cid,seed=198003,production_ps=20.,warm_ps=2.,threads=2)
