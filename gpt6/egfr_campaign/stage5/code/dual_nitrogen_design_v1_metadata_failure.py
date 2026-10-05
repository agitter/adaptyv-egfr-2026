"""Finite classical carboxylate redesign aimed at both histidine nitrogens.
Angular donor/acceptor geometry and steric packing are search heuristics only.
Explicit HIE/HID/HIP calculations must subsequently test proton specificity.
Known binders are not seeds: parent is our own untested de novo C00003.
"""
from pathlib import Path
import sys,json,copy,math,itertools,re,time,hashlib
import numpy as np
from scipy.spatial import cKDTree
S=Path('/mnt/data/egfr_campaign/stage5');R=S.parent
sys.path[:0]=[str(S/'code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from build_return_loops import load_parent,PARENT,SOURCE,sha
from classical_design import LIB,BB,conformer,atom_data,pair_score
from mm_refine import app,u
from legacy_geometry import dihedral

def geom_score(rows,his):
 out=[]
 for hp,ha in his.items():
  best=[]
  for hn,neighbors in [('ND1',['CG','CE1']),('NE2',['CD2','CE1'])]:
   h=ha[hn];dv=h-np.mean([ha[n] for n in neighbors],0);dv/=np.linalg.norm(dv);hits=[]
   for row in rows:
    if not row['loop'] or row['aa'] not in 'DE':continue
    a={k:np.asarray(v) for k,v in row['atoms'].items()};cn='CG' if row['aa']=='D' else 'CD';ns=['OD1','OD2'] if row['aa']=='D' else ['OE1','OE2']
    for name in ns:
     v=a[name]-h;dist=np.linalg.norm(v);v/=max(dist,1e-8);av=a[name]-a[cn];av/=np.linalg.norm(av);orient=max(0.,float(dv@v))*max(0.,float(av@(-v)));q=math.exp(-((dist-2.85)/.55)**2)*orient;hits.append((q,row['position'],name,float(dist),orient))
   q=max(hits,default=(0.,None,None,None,None));best.append(q)
  score=sum(q[0] for q in best)+2*min(q[0] for q in best)
  out.append({'human_histidine':hp,'N_contacts':dict(zip(['ND1','NE2'],best)),'two_nitrogen_geometric_score':score})
 return sum(r['two_nitrogen_geometric_score']*(2. if r['human_histidine']==383 else 1.) for r in out),out

def main():
 parent,rows,tx=load_parent();p=app.PDBFile(str(SOURCE));x=np.asarray(p.positions.value_in_unit(u.angstrom));target=[];his={}
 from Bio.SeqUtils import seq1
 for res in p.topology.residues():
  if res.chain.id!='A':continue
  aa=seq1(res.name,custom_map={'HID':'H','HIE':'H','HIP':'H'});a={a.name:x[a.index] for a in res.atoms() if a.element!=app.element.hydrogen};target.append({'aa':aa,'atoms':a})
  if int(res.id)+333 in [358,383]:his[int(res.id)+333]=a
 assert set(his)=={358,383}
 te=np.concatenate([atom_data(r['aa'],r['atoms']) for r in target]);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());me=np.concatenate([atom_data(r['aa'],r['atoms']) for r in mouse]);gt=cKDTree([g['xyz'] for g in json.loads((R/'stage2/intermediate/target_glycans.json').read_text())])
 edits=[{54:'E'},{52:'D',54:'E'},{52:'E',54:'E'},{53:'D',54:'E'},{53:'E',54:'E'},{54:'E',103:'Q'},{52:'E',54:'E',103:'Q'},{52:'D',54:'E',103:'Q'},{103:'D'},{104:'D'},{104:'E'},{105:'D'},{105:'E'},{103:'Q',105:'E'},{54:'E',105:'E'},{54:'E',103:'Q',105:'E'},{103:'D',105:'E'},{106:'D'},{54:'E',106:'D'}]
 history={}
 for stage in ['stage2','stage3','stage4','stage5']:
  for q in (R/stage/'intermediate/designs').glob('*.json'):
   z=json.loads(q.read_text());history.setdefault(z['sequence'],[]).append(str(q))
 report=[];serial=0
 for edit in edits:
  fixed=np.concatenate([atom_data(r['aa'],r['atoms']) for r in rows if r['position'] not in edit]);opts={};failed=False
  for pos,aa in edit.items():
   bb=np.array([rows[pos-1]['atoms'][n] for n in BB]);phi=math.degrees(dihedral(np.array(rows[pos-2]['atoms']['C']),*bb[:3]))
   if phi>0:failed=True;break
   candidates=[]
   for rn in range(len(LIB[aa]['xyz'])):
    a=conformer(aa,rn,bb);dat=atom_data(aa,a,LIB[aa]['names']);rp,vw,hb=pair_score(dat,fixed);tt=pair_score(dat,te);mt=pair_score(dat,me)
    if rp>15 or max(tt[0],mt[0])>12 or (len(dat) and gt.query(dat[:,:3])[0].min()<2.5):continue
    rr=copy.deepcopy(rows);rr[pos-1].update(aa=aa,atoms={k:np.asarray(v).tolist() for k,v in a.items()});geom,contacts=geom_score(rr,his);packing=rp+vw+hb+max(sum(tt),sum(mt));val=packing-4.*geom
    candidates.append((val,packing,rn,a,dat))
   opts[pos]=sorted(candidates,key=lambda z:z[0])[:4]
   if not opts[pos]:failed=True;break
  if failed:report.append({'edits':edit,'status':'no stereochemical/steric rotamer'});continue
  best=None
  for combo in itertools.product(*opts.values()):
   rr=copy.deepcopy(rows);packing=sum(c[1] for c in combo);severe=False
   for (pos,aa),o in zip(edit.items(),combo):rr[pos-1].update(aa=aa,atoms={k:np.asarray(v).tolist() for k,v in o[3].items()})
   for left,right in itertools.combinations(combo,2):
    z=pair_score(left[4],right[4]);packing+=sum(z)
    if z[0]>15:severe=True
   if severe:continue
   geom,contacts=geom_score(rr,his);val=packing-4*geom
   if best is None or val<best[0]:best=(val,rr,combo,contacts,geom)
  if best is None:report.append({'edits':edit,'status':'new sidechain pair clash'});continue
  value,rr,combo,contacts,geom=best;seq=''.join(r['aa'] for r in rr)
  if seq in history:report.append({'edits':edit,'status':'existing sequence not counted again','existing':history[seq]});continue
  serial+=1;cid=f'H{serial:05d}';d=copy.deepcopy(parent);d.update(candidate_id=cid,sequence=seq,structure=rr,cdr_sequences=[seq[a:b] for a,b in parent['cdr_intervals_zero_based']],seed=None,branch='own_design_dual_histidine_nitrogen_contact_hypothesis',campaign_stage=5,mutation_selection='deterministic finite rotamer enumeration',metrics={'scores_valid':False,'note':'All inherited scores invalidated; geometric search objective is not pH-binding evidence.'},ancestry={'parent_candidate':'C00003','parent_type':'own untested de novo design','parent_file':str(PARENT),'parent_sha256':sha(PARENT),'relaxed_backbone_source':str(SOURCE),'relaxed_backbone_sha256':sha(SOURCE)},contact_redesign={'edits':[{'parent_position':pos,'from':parent['sequence'][pos-1],'to':aa,'rotamer':o[2],'rotamer_source':LIB[aa]['sources'][o[2]]} for (pos,aa),o in zip(edit.items(),combo)],'geometric_search_objective':float(value),'contacts':contacts,'hypothesis':'Longer D54E and additional carboxylates may reach a second His N; actual neutral-tautomer competition requires explicit testing.'})
  (S/'intermediate/designs'/f'{cid}.json').write_text(json.dumps(d));history[seq]=[cid];report.append({'candidate_id':cid,'edits':edit,'score':value,'dual_N_score':geom,'status':'generated unvalidated hypothesis'});print(cid,edit,round(value,2),round(geom,3),flush=True)
 (S/'reference/dual_nitrogen_design.json').write_text(json.dumps(report,indent=2));print('GENERATED',serial,flush=True)
if __name__=='__main__':main()
