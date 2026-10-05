from pathlib import Path
import sys,json,time,hashlib
import numpy as np
from cluster_validate import expand
from mixed_tensor import summarize
R=Path('/mnt/data/egfr_campaign');S=R/'stage6';records=[]
for p in sorted((R/'stage5/intermediate').glob('**/*acid_histidine.json')):
 d=json.loads(p.read_text());ms=d['microstates'];st=np.asarray([m['state'] for m in ms]);n=st.shape[1]
 if n<5:continue
 E=np.asarray([[e['delta_kcal_proxy'] for e in m['energy']['states']] for m in ms]);exact,C=summarize(E,d['sites']);orders=[]
 for order in [2,3]:
  samples={tuple(s):e for s,e in zip(st,E) if np.count_nonzero(s)<=order};pred=expand(samples,st,order);pr,P=summarize(pred,d['sites']);mask=np.count_nonzero(st,axis=1)>order
  r={'order':order,'basis_samples':len(samples),'held_out_states':int(mask.sum()),'max_energy_error_kcal':float(abs(pred-E)[mask].max()),'RMS_energy_error_kcal':float(np.sqrt(np.mean((pred[mask]-E[mask])**2))),'max_contrast_error_kcal':float(abs(P-C).max()),'minimum_predicted_contrast_kcal':pr['minimum_contrast_kcal'],'minimum_exact_contrast_kcal':exact['minimum_contrast_kcal']};orders.append(r);print(p.parent.name,p.stem,'order',order,r,flush=True)
 records.append({'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'site_count':n,'exact_summary':exact,'cluster_validation':orders})
(S/'reference/cluster_validation_expanded.json').write_text(json.dumps(records,indent=2))
