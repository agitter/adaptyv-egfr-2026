"""Deterministic independent reference cases for the RMSD analysis."""
from pathlib import Path
import sys,json,math
import numpy as np
S=Path('/mnt/data/egfr_campaign/stage7');sys.path.insert(0,str(S/'code'))
from analyze_thermal import fit_metrics

def main():
 ref=np.array([[0.,0,0],[1.,0,0],[0,2.,0],[0,0,3.],[1,1,1],[2.,2.,1],[0.,3.,2],[2.,0.,1]])
 mask=np.array([False,False,False,False,True,True,True,True]);tests=[]
 def ck(name,value):tests.append({'test':name,'pass':bool(value)});assert value,name
 theta=.73;rot=np.array([[math.cos(theta),-math.sin(theta),0],[math.sin(theta),math.cos(theta),0],[0,0,1.]])
 x=ref@rot+[12.,-4.,3.];m=fit_metrics(x,ref,mask)
 ck('proper_rigid_rotation_and_translation_zero_loop_RMSD',m['loop']<1e-12)
 ck('proper_rigid_rotation_and_translation_zero_framework_RMSD',m['framework']<1e-12)
 x=ref.copy();x[mask,2]+=2.3;m=fit_metrics(x,ref,mask)
 ck('known_loop_only_displacement',abs(m['loop']-2.3)<1e-12)
 ck('loop_displacement_does_not_drive_framework_fit',m['framework']<1e-12)
 x=ref.copy();x[4]+=[1.,2.,3.];m=fit_metrics(x,ref,mask)
 ck('correct_RMSD_normalization',abs(m['loop']-math.sqrt(14/4))<1e-12)
 x=ref.copy();x[:,2]*=-1;m=fit_metrics(x,ref,mask)
 ck('reflection_is_not_silently_treated_as_a_proper_rotation',m['framework']>.1)
 out={'tests':tests,'passed':sum(r['pass'] for r in tests),'all_pass':all(r['pass'] for r in tests),'qualification':'Mathematical implementation checks, not protein stability validation.'}
 (S/'reference/alignment_tests.json').write_text(json.dumps(out,indent=2));print('Alignment tests',out['passed'],'passed')
if __name__=='__main__':main()
