from pathlib import Path
import sys,json,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u
from legacy_geometry import dihedral
p=R/'stage4/intermediate/matched_controls/intermediate/refined/C00003_human6ARU_relaxed.pdb'
x=app.PDBFile(str(p));co=np.asarray(x.positions.value_in_unit(u.angstrom));rows=[{a.name:co[a.index] for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'} for r in x.topology.residues() if r.chain.id=='B'];target={int(r.id)+333:{a.name:co[a.index] for a in r.atoms() if a.element!=app.element.hydrogen} for r in x.topology.residues() if r.chain.id=='A'}
d=json.loads((R/'stage3/intermediate/designs/C00003.json').read_text());hs=[370,383,409,433]
print('CDRs',d['cdr_intervals_zero_based'])
for i in range(len(rows)):
 if not d['structure'][i]['loop']:continue
 a=rows[i];phi=math.degrees(dihedral(rows[i-1]['C'],a['N'],a['CA'],a['C']));psi=math.degrees(dihedral(a['N'],a['CA'],a['C'],rows[i+1]['N']))
 hc={h:round(min(np.linalg.norm(v-target[h][n]) for v in a.values() for n in ['ND1','NE2']),2) for h in hs if h in target and 'ND1' in target[h]}
 print(i+1,d['sequence'][i],round(phi),round(psi),hc)
for a,b in [(54,66),(103,116),(26,37),(50,66)]:
 print('anchors',a,b,'CA distance',np.linalg.norm(rows[a-1]['CA']-rows[b-1]['CA']))
print('Matched parent exists',p)
(S/'reference/parent_his_contacts.json').write_text(json.dumps({'source':str(p),'cdrs':d['cdr_intervals_zero_based'],'anchors':[[54,66],[103,116]],'protected_contacts':['D54','Y100','E103'],'proposed_rebuild':'noncontact loop-return segments, not wholesale removal of pH contacts'},indent=2))
