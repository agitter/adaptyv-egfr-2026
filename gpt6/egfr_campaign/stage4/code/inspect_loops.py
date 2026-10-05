"""Inspect classical phi/psi geometry; not a stability predictor."""
from pathlib import Path
import sys,json,math,itertools
import numpy as np
R=Path('/mnt/data/egfr_campaign'); S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u
from legacy_geometry import dihedral
from classical_design import conformer,LIB,BB
from scipy.spatial import cKDTree
out={}
for cid in ['C00003','C00009']:
 d=json.loads((R/'stage3/intermediate/designs'/f'{cid}.json').read_text())
 p=app.PDBFile(str(R/'stage3/intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb'))
 xyz=np.asarray(p.positions.value_in_unit(u.angstrom));rs=[r for r in p.topology.residues() if r.chain.id=='B']
 rows=[{a.name:xyz[a.index] for a in r.atoms() if a.element!=app.element.hydrogen} for r in rs]
 tt=np.array([xyz[a.index] for a in p.topology.atoms() if a.residue.chain.id=='A' and a.element!=app.element.hydrogen]);tree=cKDTree(tt)
 report=[]
 for i,row in enumerate(rows):
  if not d['structure'][i]['loop']:continue
  phi=math.degrees(dihedral(rows[i-1]['C'],row['N'],row['CA'],row['C']))
  psi=math.degrees(dihedral(row['N'],row['CA'],row['C'],rows[i+1]['N']))
  target_min=float(tree.query(np.array(list(row.values())))[0].min())
  r={'position':i+1,'aa':d['sequence'][i],'phi':phi,'psi':psi,'target_min_A':target_min,'pro_allowed_phi':-90<=phi<=-40,'gly_to_ala_phi':phi<0}
  report.append(r)
  print(cid,i+1,r['aa'],round(phi),round(psi),'target',round(target_min,1))
 # Contact-independent geometric Cys pairs; native disulfide positions excluded.
 ds=[]
 for i in range(len(rows)):
  if not d['structure'][i]['loop']:continue
  for j in range(i+2,len(rows)):
   if d['sequence'][i] in 'CP' or d['sequence'][j] in 'CP':continue
   if i+1 in [52,54,100,103] or j+1 in [52,54,100,103]:continue
   if np.linalg.norm(rows[i]['CA']-rows[j]['CA'])>7.2:continue
   best=None
   for ri,rj in itertools.product(range(len(LIB['C']['xyz'])),repeat=2):
    ai=conformer('C',ri,np.array([rows[i][n] for n in BB]));aj=conformer('C',rj,np.array([rows[j][n] for n in BB]))
    sep=float(np.linalg.norm(ai['SG']-aj['SG']));chi=math.degrees(dihedral(ai['CB'],ai['SG'],aj['SG'],aj['CB']))
    score=abs(sep-2.03)+abs(abs(chi)-90)/90
    if best is None or score<best['score']:best={'positions':[i+1,j+1],'from':d['sequence'][i]+d['sequence'][j],'rotamers':[ri,rj],'SG_distance_A':sep,'chi3_deg':chi,'score':score}
   if best and best['score']<1.2:ds.append(best)
 out[cid]={'loop_geometry':report,'disulfide_geometries':sorted(ds,key=lambda r:r['score'])}
 print('CYS',sorted(ds,key=lambda r:r['score'])[:20])
(S/'reference/loop_inspection.json').write_text(json.dumps(out,indent=2))
