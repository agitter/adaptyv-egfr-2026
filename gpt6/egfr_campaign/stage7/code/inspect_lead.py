"""Classical torsion/steric audit of the prior released lead, before design.
No learned model, inferred affinity, or alteration of previous results.
"""
from pathlib import Path
import sys,json,hashlib,math,copy
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/x/'code') for x in ['stage6','stage5','stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u
from classical_design import LIB,BB,conformer,atom_data,pair_score
from legacy_geometry import dihedral
from Bio.SeqUtils import seq1

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_model(p):
 m=app.PDBFile(str(p));x=np.asarray(m.positions.value_in_unit(u.angstrom));br=[r for r in m.topology.residues() if r.chain.id=='B']
 if not br and len(list(m.topology.chains()))==1:br=list(m.topology.residues())
 assert len(br)==125,(str(p),len(br))
 rows=[{'position':i+1,'aa':seq1(r.name,custom_map={'HID':'H','HIE':'H','HIP':'H'}),'atoms':{a.name:x[a.index] for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}} for i,r in enumerate(br)]
 other=[{'position':int(r.id)+333,'aa':seq1(r.name,custom_map={'HID':'H','HIE':'H','HIP':'H'}),'atoms':{a.name:x[a.index] for a in r.atoms() if a.element!=app.element.hydrogen}} for r in m.topology.residues() if r not in br]
 return rows,other

def main():
 src=R/'stage6/intermediate/designs/J00008.json';d=json.loads(src.read_text());models={}
 for k,p in [('human',R/'stage6/intermediate/refined/J00008_human6ARU_relaxed.pdb'),('mouse',R/'stage6/intermediate/refined/J00008_mouseAF_relaxed.pdb'),('free',R/'stage6/intermediate/unbound/J00008_free_relaxed.pdb')]:
  rows,t=load_model(p);assert ''.join(r['aa'] for r in rows)==d['sequence'];models[k]=(rows,t)
 rec=[]
 for i,row in enumerate(d['structure']):
  if not row['loop']:continue
  vals={}
  for model,(rows,tar) in models.items():
   r=rows[i];a=r['atoms'];ph=math.degrees(dihedral(rows[i-1]['atoms']['C'],a['N'],a['CA'],a['C']));ps=math.degrees(dihedral(a['N'],a['CA'],a['C'],rows[i+1]['atoms']['N']));v={'phi':ph,'psi':ps}
   if r['aa']=='G':
    ala=conformer('A',0,np.asarray([a[n] for n in BB]));dat=atom_data('A',ala,LIB['A']['names']);env=np.concatenate([atom_data(rr['aa'],rr['atoms']) for j,rr in enumerate(rows) if j!=i]);v['A_nonlocal_cb_min_A']=float(cKDTree([z for j,rr in enumerate(rows) if abs(i-j)>1 for z in rr['atoms'].values()]).query(ala['CB'])[0]);v['A_intramol_score']=list(pair_score(dat,env));v['A_target_score']=list(pair_score(dat,np.concatenate([atom_data(rr['aa'],rr['atoms']) for rr in tar]))) if tar else [0.,0.,0.]
   if tar:
    pts=np.concatenate([list(rr['atoms'].values()) for rr in tar]);v['target_min_A']=float(cKDTree(pts).query(np.asarray(list(a.values())))[0].min())
   vals[model]=v
  rec.append({'position':i+1,'aa':row['aa'],'models':vals})
 out={'source':str(src),'source_sha256':sha(src),'sequence':d['sequence'],'records':rec,'hypothesis':'Selected Gly to Ala/Ser substitutions can reduce unfolded conformational freedom; this is not a measured stability benefit. All pH-contact residues initially preserved.'}
 (S/'reference/lead_torsions_and_cb_clearance.json').write_text(json.dumps(out,indent=2))
 for r in rec:
  if r['aa']=='G':
   print(r['position'],'G', 'phi',[round(r['models'][k]['phi'],1) for k in models], 'psi',[round(r['models'][k]['psi'],1) for k in models], 'Amin',[round(r['models'][k]['A_nonlocal_cb_min_A'],2) for k in models],'pack',[round(r['models'][k]['A_intramol_score'][0],2) for k in models],flush=True)
 print('CDRs',d['cdr_sequences']);print('SOURCES',len(models));print('LIB sizes',{a:len(LIB[a]['xyz']) for a in 'AGSNTVHP'})
if __name__=='__main__':main()
