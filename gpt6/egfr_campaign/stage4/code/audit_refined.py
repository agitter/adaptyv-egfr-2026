"""Independent stereochemistry, peptide and disulfide strain audit.
Checks geometry, not folding or correct oxidation in a synthesis experiment.
"""
from pathlib import Path
import sys,json,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u
from legacy_geometry import dihedral

def angle(a,b,c):
 v=np.asarray(a)-b;w=np.asarray(c)-b;return math.degrees(math.acos(float(np.clip(np.dot(v,w)/np.linalg.norm(v)/np.linalg.norm(w),-1,1))))
def audit(path,design,chain='B'):
 p=app.PDBFile(str(path));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));rows=[{a.name:xyz[a.index] for a in r.atoms()} for r in p.topology.residues() if r.chain.id==chain];assert len(rows)==len(design['sequence']);inversions=[];omega=[];cystines=[];strain=[]
 for i,r in enumerate(rows):
  if 'CB' in r and np.linalg.det(np.array([r[n]-r['CA'] for n in ['N','C','CB']]))<=0:inversions.append(i+1)
  if i:
   theta=math.degrees(dihedral(rows[i-1]['CA'],rows[i-1]['C'],r['N'],r['CA']));deviation=abs(180-abs(theta))
   if deviation>25:omega.append({'previous_position':i,'next_position':i+1,'omega_deg':theta,'next_is_proline':design['sequence'][i]=='P'})
 for a,b in p.topology.bonds():
  if a.name==b.name=='SG' and a.residue.chain.id==b.residue.chain.id==chain:
   i,j=int(a.residue.id),int(b.residue.id);ai,aj=rows[i-1],rows[j-1];dd=float(np.linalg.norm(ai['SG']-aj['SG']));chi=math.degrees(dihedral(ai['CB'],ai['SG'],aj['SG'],aj['CB']));aa=[angle(ai['CB'],ai['SG'],aj['SG']),angle(ai['SG'],aj['SG'],aj['CB'])];bad=not(1.9<=dd<=2.2 and abs(abs(chi)-90)<=35 and all(80<=x<=130 for x in aa));row={'positions':sorted([i,j]),'SG_distance_A':dd,'chi3_deg':chi,'CB_S_S_angles_deg':aa,'strain_gate_pass':not bad};cystines.append(row)
   if bad:strain.append(row)
 present={tuple(x['positions']) for x in cystines};intended=design.get('intended_extra_disulfide_pairs',[]);extra_ok=all(tuple(sorted(x)) in present for x in intended)
 return {'candidate_id':design['candidate_id'],'absolute_alpha_inversions':inversions,'peptide_omega_warnings':omega,'actual_cystines':cystines,'disulfide_strain_warnings':strain,'intended_extra_pairs_present':extra_ok,'pass':not inversions and not omega and extra_ok and not strain,'thresholds':{'omega_max_deviation_from_trans_deg':25,'SS_A':[1.9,2.2],'absolute_chi3_deviation_from90_max_deg':35,'CB_S_S_angles_deg':[80,130]},'qualification':'Heuristic geometry gates; not a folding, oxidation-yield or binding prediction.'}
def main():
 records=[]
 for root,pattern in [(S,'*_relaxed.pdb'),(R/'stage3','C0000[39]_human6ARU_relaxed.pdb')]:
  for p in sorted((root/'intermediate/refined').glob(pattern)):
   cid=p.name.split('_')[0];dp=root/'intermediate/designs'/(cid+'.json')
   if not dp.exists():continue
   d=json.loads(dp.read_text());r=audit(p,d);r['source']=str(p);out=S/'intermediate/independent_v2';out.mkdir(exist_ok=True);(out/(p.stem+'.json')).write_text(json.dumps(r,indent=2));records.append(r);print(p.stem,r['pass'],'omega',len(r['peptide_omega_warnings']),'SSstrain',len(r['disulfide_strain_warnings']),flush=True)
 (S/'reference/independent_geometry_v2_summary.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
