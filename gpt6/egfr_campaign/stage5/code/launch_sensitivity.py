"""Run declared sensitivity cases locally and record all worker exits."""
from pathlib import Path
import json,subprocess,sys,time
S=Path('/mnt/data/egfr_campaign/stage5')
# Wait for the finite mouse batch to release its slot; this driver completes
# during the current campaign response and must not be left running.
p=S/'reference/mouse_batch_execution.json'
start=time.time()
while not p.exists():
 if time.time()-start>900:raise TimeoutError('Mouse batch did not finish')
 time.sleep(2)
ids=['B00000','N2070600','N2091100'];tasks=[{'log':f'sensitivity_{cid}_{mode}.log','args':[str(S/'code/sensitivity.py'),cid,mode]} for cid in ids for mode in ['cap','acid']]
(S/'reference/sensitivity_selection.json').write_text(json.dumps({'ids':ids,'modes':['native-coordinate target ACE/NME cap','explicit protonation of both retained designed carboxylates'],'same_acid_residues_via_parent_position_map':True,'not_a_pKa_prediction':True},indent=2))
config=S/'reference/sensitivity_batch.json';config.write_text(json.dumps(tasks,indent=2));subprocess.run([sys.executable,str(S/'code/launch.py'),str(config),'1'],check=True)
