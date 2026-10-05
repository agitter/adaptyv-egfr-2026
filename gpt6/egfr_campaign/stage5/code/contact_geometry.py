"""Continuous contact geometry, not a declaration of hydrogen bonding.
Includes both histidine nitrogens so neutral-tautomer alternatives are visible.
"""
from pathlib import Path
import sys,json,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u

def run(path):
 cid=path.name.split('_')[0];d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());p=app.PDBFile(str(path));x=np.asarray(p.positions.value_in_unit(u.angstrom));rows={};target={}
 for r in p.topology.residues():
  if r.chain.id=='B':rows[int(r.id)]={a.name:x[a.index] for a in r.atoms()}
  if r.chain.id=='A' and r.name in ['HIS','HID','HIE','HIP']:target[int(r.id)+333]={a.name:x[a.index] for a in r.atoms()}
 record=[]
 for hp,hr in target.items():
  if hp not in [358,383]:continue
  neighbors=[]
  for row in d['structure']:
   if not row['loop'] or row['aa'] not in 'DEY':continue
   br=rows[row['position']];names=[n for n in br if n in ['OD1','OD2','OE1','OE2','OH']]
   for n in names:
    for hn in ['ND1','NE2']:
     distance=float(np.linalg.norm(br[n]-hr[hn]))
     if distance<=7.:neighbors.append({'binder_position':row['position'],'parent_position':row.get('parent_position',row['position']),'binder_aa':row['aa'],'acceptor_atom':n,'histidine_N':hn,'heavy_distance_A':distance})
  record.append({'target_human_position':hp,'nearby_candidate_oxygens':sorted(neighbors,key=lambda z:z['heavy_distance_A'])})
 out=S/'intermediate/contacts';out.mkdir(exist_ok=True);r={'candidate_id':cid,'source':str(path),'contacts':record,'qualification':'Distances alone do not establish proton-specific hydrogen bonds; original neutral tautomers remain competitors.'};(out/(path.stem+'.json')).write_text(json.dumps(r,indent=2));return r
if __name__=='__main__':
 for cid in ['B00000','N2070600','N2091100']:
  path=S/'intermediate/refined'/f'{cid}_human6ARU_relaxed.pdb';r=run(path)
  for h in r['contacts']:print(cid,h['target_human_position'],h['nearby_candidate_oxygens'][:4],flush=True)
