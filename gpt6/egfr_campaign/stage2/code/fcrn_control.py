"""Pre-2011 acid-on reference; used only to test a scoring model, never as a seed."""
from pathlib import Path
import json,copy
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from proton_polynomial import evaluate_design
from classical_design import S,R
p=R/'inputs/transfer2/egfr_refinement_fetch/1I1A.pdb'
st=PDBParser(QUIET=True).get_structure('1I1A',str(p));chains={}
for chain in st[0]:
 rows=[]
 for r in chain:
  aa=seq1(r.resname)
  if r.id[0]!=' ' or aa=='X' or 'CA' not in r:continue
  rows.append({'aa':aa,'id':r.id[1],'atoms':{a.name:a.coord.tolist() for a in r if a.element!='H'}})
 chains[chain.id]=rows
binder=chains['C']+chains['D'];target=chains['A']+chains['B']
res=evaluate_design(binder,target,True)
res['source']='PDB 1I1A (2001), chains C+D Fc; A+B FcRn/beta2 microglobulin'
res['role']='qualitative acid-on control, no calibration of absolute affinity or screening thresholds'
res['nontitrating_control']='Setting all histidine protonation energy differences to zero gives zero pH contrast; covered by polynomial tests.'
(S/'reference/fcrn_control.json').write_text(json.dumps(res,indent=2))
print('FcRn model acid-on across scenarios',res['all_scenarios_acid_on'],'contrast interval',res['worst_contrast_kcal_surrogate'],res['best_contrast_kcal_surrogate']);print('sites',res['site_count'],res['site_details'])
