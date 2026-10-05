"""Geometric probe for additional acid-on proton uptake, not a sequence generator.
Classical finite rotamers and directional donor/acceptor geometry only.
"""
from pathlib import Path
import json,sys,math
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(S/'code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from classical_design import LIB,BB,conformer,atom_data,pair_score
from mm_refine import app,u
from legacy_geometry import dihedral
from Bio.SeqUtils import seq1

def main():
 allout=[]
 for root,cid in [(S,'S00028'),(S/'intermediate/matched_controls','C00003')]:
  source=root/'intermediate/refined'/(cid+'_human6ARU_relaxed.pdb');d=json.loads((root/'intermediate/designs'/(cid+'.json')).read_text());p=app.PDBFile(str(source));x=np.asarray(p.positions.value_in_unit(u.angstrom));res=list(p.topology.residues());rows=[{'aa':seq1(r.name),'atoms':{a.name:x[a.index] for a in r.atoms() if a.element!=app.element.hydrogen}} for r in res if r.chain.id=='B'];target=[{'aa':seq1(r.name),'human_pos':int(r.id)+333,'atoms':{a.name:x[a.index] for a in r.atoms() if a.element!=app.element.hydrogen}} for r in res if r.chain.id=='A'];tenv=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);gly=np.array([g['xyz'] for g in json.loads((R/'stage2/intermediate/target_glycans.json').read_text())]);gt=cKDTree(gly);oxy=[]
  for r in target:
   for n in (['OD1','OD2'] if r['aa']=='D' else ['OE1','OE2'] if r['aa']=='E' else []):oxy.append((r['human_pos'],n,r['atoms'][n]))
  targetN=np.array([r[2] for r in oxy]);options=[]
  for i,row in enumerate(rows):
   if not d['structure'][i]['loop'] or i+1 in [54,100,103] or row['aa'] in 'CP':continue
   ph=math.degrees(dihedral(rows[i-1]['atoms']['C'],row['atoms']['N'],row['atoms']['CA'],row['atoms']['C']))
   if row['aa']=='G' and ph>0:continue
   env=np.concatenate([atom_data(r['aa'],r['atoms']) for j,r in enumerate(rows) if j!=i]);bb=np.array([row['atoms'][n] for n in BB]);at=[]
   for k in range(len(LIB['H']['xyz'])):
    a=conformer('H',k,bb);data=atom_data('H',a,LIB['H']['names']);rp,vw,hb=pair_score(data,env);tr,tv,th=pair_score(data,tenv)
    if rp>12 or tr>8 or gt.query(data[:,:3])[0].min()<2.5:continue
    contacts=[]
    for n,adj in [('ND1',['CG','CE1']),('NE2',['CD2','CE1'])]:
     vector=a[n]-(a[adj[0]]+a[adj[1]])/2;vector/=np.linalg.norm(vector)
     dist=np.linalg.norm(targetN-a[n],axis=1)
     for j in np.where(dist<4.)[0]:
      cosine=float(vector@(targetN[j]-a[n])/dist[j]);contacts.append({'His_atom':n,'target_position':oxy[j][0],'target_atom':oxy[j][1],'NO_distance_A':float(dist[j]),'donor_alignment_cosine':cosine})
    good=[r for r in contacts if r['NO_distance_A']<3.6 and r['donor_alignment_cosine']>.6]
    if not good:continue
    both=len({r['His_atom'] for r in good})==2;two=len({r['target_position'] for r in good})>=2
    item={'parent':cid,'position':i+1,'from':row['aa'],'rotamer':k,'rotamer_source':LIB['H']['sources'][k],'parent_phi_deg':ph,'packing_score':float(rp+vw+hb+tr+tv+th),'good_contacts':good,'both_His_nitrogens':both,'two_distinct_acid_residues':two};at.append(item)
   if at:options.extend(sorted(at,key=lambda z:(not(z['both_His_nitrogens'] and z['two_distinct_acid_residues']),z['packing_score']))[:2])
  allout.extend(options);print(cid,len(options),'dual contacts',sum(r['both_His_nitrogens'] and r['two_distinct_acid_residues'] for r in options),flush=True)
  for r in options:print(r['position'],r['from'],r['both_His_nitrogens'],r['two_distinct_acid_residues'],[(g['target_position'],g['His_atom'],round(g['NO_distance_A'],2)) for g in r['good_contacts']],flush=True)
 (S/'reference/histidine_opportunities.json').write_text(json.dumps({'role':'proposal only; no new sequences created or affinity claim','options':allout},indent=2))
if __name__=='__main__':main()
