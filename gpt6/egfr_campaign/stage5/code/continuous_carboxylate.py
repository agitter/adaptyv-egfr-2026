"""Powell1964 side-chain chi optimization for a two-acceptor histidine hypothesis.
A geometry restraint is used only to propose conformers. It is absent in the
subsequent Amber/OBC interface refinement and proton-state calculation.
"""
from pathlib import Path
import sys,json,copy,math,time,itertools
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
S=Path('/mnt/data/egfr_campaign/stage5');R=S.parent;sys.path[:0]=[str(S/'code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from build_return_loops import load_parent,PARENT,SOURCE,sha
from classical_design import LIB,BB,conformer,atom_data,pair_score,charge_center,screened_charge
from mm_refine import app,u
from Bio.SeqUtils import seq1

def rotate(a,left,right,names,theta):
 a={k:np.array(v) for k,v in a.items()};base=a[left];axis=a[right]-base;axis/=np.linalg.norm(axis);co=math.cos(theta);si=math.sin(theta)
 for name in names:
  v=a[name]-base;a[name]=base+co*v+si*np.cross(axis,v)+(1-co)*(v@axis)*axis
 return a

def modify(a,x):
 a=rotate(a,'CA','CB',['CG','CD','OE1','OE2'],x[0]);a=rotate(a,'CB','CG',['CD','OE1','OE2'],x[1]);return rotate(a,'CG','CD',['OE1','OE2'],x[2])

def main():
 parent,rows,tx=load_parent();p=app.PDBFile(str(SOURCE));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));target=[];hist=None
 for r in p.topology.residues():
  if r.chain.id!='A':continue
  a={a.name:xyz[a.index] for a in r.atoms() if a.element!=app.element.hydrogen};target.append({'aa':seq1(r.name,custom_map={'HID':'H','HIE':'H','HIP':'H'}),'atoms':a})
  if int(r.id)+333==358:hist=a
 assert hist is not None
 te=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);mo=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());me=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mo]);gly=cKDTree([r['xyz'] for r in json.loads((R/'stage2/intermediate/target_glycans.json').read_text())]);hn=hist['NE2'];direction=hn-(hist['CD2']+hist['CE1'])/2;direction/=np.linalg.norm(direction);report=[]
 for cid,precursor in [('H00020','H00014'),('H00021','H00016')]:
  # Build neutral E103Q, but retain Q105: precursor supplies only Q103 and optional E54 geometry.
  base=copy.deepcopy(rows);old=json.loads((S/'intermediate/designs'/f'{precursor}.json').read_text())
  for pos in ([103] if cid=='H00020' else [54,103]):base[pos-1]=copy.deepcopy(old['structure'][pos-1])
  fixed=np.concatenate([atom_data(r['aa'],r['atoms']) for r in base if r['position']!=104]);charges=[(r['aa'],charge_center(r['aa'],r['atoms'])) for r in base if r['position']!=104 and r['aa'] in 'DEKR'];bb=np.array([base[103]['atoms'][n] for n in BB]);best=None;trials=[]
  for rn in range(len(LIB['E']['xyz'])):
   init=conformer('E',rn,bb)
   def score(x,full=False):
    a=modify(init,x);dat=atom_data('E',a,LIB['E']['names']);rep,vw,hb=pair_score(dat,fixed);tr=pair_score(dat,te);mr=pair_score(dat,me);contacts=[]
    for o in ['OE1','OE2']:
     v=a[o]-hn;dist=float(np.linalg.norm(v));v/=max(dist,1e-8);av=a[o]-a['CD'];av/=np.linalg.norm(av);orient=max(0.,float(direction@v))*max(0.,float(av@(-v)));contacts.append((12*(dist-2.85)**2+6*(1-orient)**2,dist,orient,o))
    contact=min(contacts);gmin=float(gly.query(dat[:,:3])[0].min());wall=100*max(0,2.5-gmin)**2;cc=charge_center('E',a);electro=sum(float(screened_charge(-1.,-1. if aa in 'DE' else 1.,np.linalg.norm(cc-c),40.)) for aa,c in charges);val=rep+vw+hb+max(sum(tr),sum(mr))+electro+contact[0]+wall
    if full:return val,a,rep,tr[0],mr[0],gmin,contact
    return val
   opt=minimize(score,np.zeros(3),method='Powell',options={'maxiter':70,'xtol':1e-4,'ftol':1e-5});value,a,rep,tr,mr,gmin,contact=score(opt.x,True);valid=rep<15 and max(tr,mr)<12 and gmin>=2.5 and 2.55<contact[1]<3.3 and contact[2]>.35;trials.append({'rotamer':rn,'angles_rad':opt.x.tolist(),'value':float(value),'binder_repulsion':rep,'human_repulsion':tr,'mouse_repulsion':mr,'glycan_min_A':gmin,'NE2_contact_distance_A':contact[1],'orientation':contact[2],'pass_proposal_gate':valid})
   if valid and (best is None or value<best[0]):best=(value,a,rn,opt.x,contact)
  record={'candidate_id':cid,'precursor':precursor,'trials':trials,'new_sequence_generated':best is not None};report.append(record)
  if best is None:print(cid,'NO acceptable continuous proposal',flush=True);continue
  value,a,rn,angles,contact=best;base[103].update(aa='E',atoms={k:np.asarray(v).tolist() for k,v in a.items()});seq=''.join(r['aa'] for r in base);d=copy.deepcopy(parent);edits=[{'parent_position':r['position'],'from':parent['sequence'][r['position']-1],'to':r['aa']} for r in base if r['aa']!=parent['sequence'][r['position']-1]]
  d.update(candidate_id=cid,sequence=seq,structure=base,cdr_sequences=[seq[a:b] for a,b in d['cdr_intervals_zero_based']],campaign_stage=5,branch='own_design_continuous_carboxylate_proposal',seed=None,mutation_selection='deterministic empirical rotamers plus Powell1964 chi optimization',ancestry={'parent_candidate':'C00003','parent_type':'own untested de novo design','parent_file':str(PARENT),'parent_sha256':sha(PARENT),'relaxed_backbone_source':str(SOURCE),'relaxed_backbone_sha256':sha(SOURCE),'neutral_carboxylate_geometry_source':str(S/'intermediate/designs'/f'{precursor}.json')},metrics={'scores_valid':False,'note':'Geometric proposal only; explicit neutral/protonated competition not yet evaluated.'},contact_redesign={'edits':edits,'geometric_search_objective':float(value),'continuous_chi':{'rotamer':rn,'source':LIB['E']['sources'][rn],'angles_radians':angles.tolist()},'His358_NE2_distance_A':contact[1],'orientation':contact[2],'proposal_restraint_removed_before_refinement':True})
  path=S/'intermediate/designs'/f'{cid}.json';assert not path.exists();path.write_text(json.dumps(d));print(cid,edits,'NE2 contact',contact[1],contact[2],flush=True)
 (S/'reference/continuous_carboxylate_trials.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':main()
