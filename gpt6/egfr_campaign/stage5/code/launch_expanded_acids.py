from pathlib import Path
import json,time,subprocess,sys
S=Path('/mnt/data/egfr_campaign/stage5');start=time.time()
while not (S/'reference/sensitivity_batch_execution.json').exists():
 if time.time()-start>1800:raise TimeoutError('Prior sensitivity batch did not finish')
 time.sleep(2)
ids=['B00000','N2070600','N2091100','T3070400'];tasks=[{'log':f'three_acids_{cid}.log','args':[str(S/'code/sensitivity_three_acids.py'),cid]} for cid in ids]
(S/'reference/expanded_acid_scope.json').write_text(json.dumps({'reason':'Actual closest target-H358 contact is binder parent E106, not E103. Include E106 protonation instead of assuming it remains charged at both pH values. Preserve two-acid results separately.','parent_residue_ids':[54,103,106],'sites_total':5,'microstates_expected':243,'ids':ids,'no_pKa_prediction':True},indent=2));config=S/'reference/expanded_acid_batch.json';config.write_text(json.dumps(tasks,indent=2));subprocess.run([sys.executable,str(S/'code/launch.py'),str(config),'1'],check=True)
