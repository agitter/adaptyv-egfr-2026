"""Compare matched designs at the SAME assumed protonation parameters.
The scenario grid is sensitivity analysis, not independent measurements or a
probability distribution over biochemical outcomes.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';sys.path.insert(0,str(R/'stage6/code'))
from mixed_tensor import energies_at_ph

def load(cid):
 p=S/'intermediate/cluster_acids'/f'{cid}_cluster.json';d=json.loads(p.read_text());npz=p.with_name(f'{cid}_energy_tensor.npz');z=np.load(npz)
 key='energies_kcal';E=z[key];sites=d['sites'];lo=energies_at_ph(E,sites,6.5);hi=energies_at_ph(E,sites,7.4)
 assert abs(float((hi-lo).min())-d['proton_model']['minimum_contrast_kcal'])<1e-8
 return d,hi-lo,lo,hi

def main():
 parent,base,base_lo,base_hi=load('M00000');rows=[]
 for cid in ['M00018','M00021']:
  d,x,lo,hi=load(cid)
  def signature(s):return (s['chain'],str(s['pdb_residue_id']),s['kind'])
  assert [signature(s) for s in d['sites']]==[signature(s) for s in parent['sites']]
  delta=x-base
  rows.append({'candidate_id':cid,'scenario_count':int(delta.size),'minimum_contrast_kcal':float(x.min()),'matched_parent_minimum_contrast_kcal':float(base.min()),'difference_of_minima_kcal':float(x.min()-base.min()),'paired_parameter_contrast_change_min_kcal':float(delta.min()),'paired_parameter_contrast_change_max_kcal':float(delta.max()),'paired_parameter_contrast_change_median_kcal':float(np.median(delta)),'number_of_parameter_settings_with_larger_contrast':int((delta>0).sum()),'acid_interaction_proxy_change_range_kcal':[float((lo-base_lo).min()),float((lo-base_lo).max())],'neutral_interaction_proxy_change_range_kcal':[float((hi-base_hi).min()),float((hi-base_hi).max())],'scope':'Fixed-coordinate classical conditional model with matched priors; no predicted Kd, binding ratio or nondetection claim.'})
 out={'parent':'M00000, unchanged J00008 with matched additional refinement','comparisons':rows,'scenario_grid_is_not_probabilistic':True};(S/'reference/paired_protonation_comparison.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
