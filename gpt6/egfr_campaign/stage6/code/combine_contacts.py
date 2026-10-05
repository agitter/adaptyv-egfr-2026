"""Classical finite-rotamer factorial experiment on our own de novo parent.
Shared backbone and unchanged-sequence controls; no established binder loops.
Search geometry is never treated as measured affinity or pH-switch evidence.
"""
from pathlib import Path
import sys,json,copy,math,itertools,hashlib,re,time
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(R/x/'code') for x in ['stage5','stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
from classical_design import LIB,BB,conformer,atom_data,pair_score
from dual_nitrogen_design import geom_score
from mm_refine import app,u
from Bio.SeqUtils import seq1
from legacy_geometry import dihedral

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 start=time.time();source=R/'stage5/intermediate/refined/B00000_human6ARU_relaxed.pdb';parentfile=R/'stage5/intermediate/designs/B00000.json';parent=json.loads(parentfile.read_text());p=app.PDBFile(str(source));x=np.asarray(p.positions.value_in_unit(u.angstrom));rows=copy.deepcopy(parent['structure']);target=[];his={};br=[r for r in p.topology.residues() if r.chain.id=='B']
 for row,res in zip(rows,br):row['atoms']={a.name:x[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
 for res in p.topology.residues():
  if res.chain.id!='A':continue
  aa=seq1(res.name,custom_map={'HID':'H','HIE':'H','HIP':'H'});a={a.name:x[a.index] for a in res.atoms() if a.element!=app.element.hydrogen};target.append({'aa':aa,'atoms':a})
  if int(res.id)+333 in [358,383]:his[int(res.id)+333]=a
 assert set(his)=={358,383};te=np.concatenate([atom_data(t['aa'],t['atoms']) for t in target]);mouse=json.loads((R/'intermediate/mouseAF_aligned.json').read_text());me=np.concatenate([atom_data(t['aa'],t['atoms']) for t in mouse]);gt=cKDTree([g['xyz'] for g in json.loads((R/'stage2/intermediate/target_glycans.json').read_text())])
 controls=[{}, {52:'E'}, {52:'E',53:'D'}, {104:'E'}, {54:'E',103:'Q',104:'E'}]
 experimental=[{52:'E',104:'E'},{53:'D',104:'E'},{53:'E',104:'E'},{52:'E',53:'D',104:'E'},{52:'E',53:'E',104:'E'},{52:'D',104:'E'},{52:'E',103:'Q',104:'E'},{53:'E',103:'Q',104:'E'},{52:'E',53:'D',103:'Q',104:'E'},{52:'E',54:'E',104:'E'},{52:'E',54:'E',103:'Q',104:'E'},{53:'D',54:'E',103:'Q',104:'E'},{52:'D',54:'E',103:'Q',104:'E'},{52:'E',53:'D',54:'E',103:'Q',104:'E'},{52:'E',53:'E',54:'E',103:'Q',104:'E'},{52:'D',53:'E',103:'Q',104:'E'}]
 plan={'parent_sequence_source':str(parentfile),'parent_sha256':sha(parentfile),'shared_backbone_source':str(source),'shared_backbone_sha256':sha(source),'controls':controls,'experimental':experimental,'method':'Deterministic finite empirical rotamer search; factorial combination of earlier physical hypotheses. All energetic metrics recomputed; modern implementations of historical algorithms only.','refinement_iterations':450,'unchanged_controls_receive_same_refinement':True,'selection':'Evaluate every sterically feasible generated sequence with human MM; prioritize followup by strict geometry plus conservative independent-prior contrast, not search score.'};(S/'reference/combination_predeclaration.json').write_text(json.dumps(plan,indent=2))
 history=json.loads((S/'reference/sequence_lineages.json').read_text());report=[];out=S/'intermediate/designs';out.mkdir(exist_ok=True)
 for serial,edit in enumerate(controls+experimental):
  cid=f'J{serial:05d}';fixed=np.concatenate([atom_data(r['aa'],r['atoms']) for r in rows if r['position'] not in edit]);options={};failed=None
  for pos,aa in edit.items():
   assert rows[pos-1]['loop'];bb=np.array([rows[pos-1]['atoms'][n] for n in BB]);phi=math.degrees(dihedral(np.array(rows[pos-2]['atoms']['C']),*bb[:3]))
   if phi>0:failed='positive phi for acidic substitution';break
   candidates=[]
   for rn in range(len(LIB[aa]['xyz'])):
    a=conformer(aa,rn,bb);dat=atom_data(aa,a,LIB[aa]['names']);pck=pair_score(dat,fixed);hh=pair_score(dat,te);mm=pair_score(dat,me)
    if pck[0]>15 or max(hh[0],mm[0])>12 or gt.query(dat[:,:3])[0].min()<2.5:continue
    rr=copy.deepcopy(rows);rr[pos-1].update(aa=aa,atoms={k:np.asarray(v).tolist() for k,v in a.items()});geom,_=geom_score(rr,his);packing=sum(pck)+max(sum(hh),sum(mm));candidates.append((packing-4*geom,packing,rn,a,dat))
   options[pos]=sorted(candidates,key=lambda z:z[0])[:2]
   if not options[pos]:failed='no feasible rotamer';break
  if failed:report.append({'candidate_id':cid,'edits':edit,'status':failed});continue
  best=None
  for combo in itertools.product(*options.values()):
   rr=copy.deepcopy(rows);packing=sum(c[1] for c in combo);bad=False
   for (pos,aa),o in zip(edit.items(),combo):rr[pos-1].update(aa=aa,atoms={k:np.asarray(v).tolist() for k,v in o[3].items()})
   for a,b in itertools.combinations(combo,2):
    z=pair_score(a[4],b[4]);packing+=sum(z);bad=bad or z[0]>15
   if bad:continue
   geom,contacts=geom_score(rr,his);score=packing-4*geom
   if best is None or score<best[0]:best=(score,rr,combo,contacts)
  if best is None:report.append({'candidate_id':cid,'edits':edit,'status':'mutated sidechains clash'});continue
  score,rr,combo,contacts=best;seq=''.join(r['aa'] for r in rr);d=copy.deepcopy(parent)
  d.update(candidate_id=cid,sequence=seq,structure=rr,cdr_sequences=[seq[a:b] for a,b in parent['cdr_intervals_zero_based']],campaign_stage=6,branch='shared_backbone_factorial_contact_test',seed=None,mutation_selection='deterministic finite rotamer enumeration',metrics={'scores_valid':False,'note':'Inherited scores invalid. Search score is not binding evidence.'},ancestry={'parent_type':'own untested de novo design','parent_file':str(parentfile),'parent_sha256':sha(parentfile),'shared_backbone_file':str(source),'shared_backbone_sha256':sha(source)},combination={'edits':[{'position':pos,'from':parent['sequence'][pos-1],'to':aa,'rotamer':o[2],'rotamer_source':LIB[aa]['sources'][o[2]]} for (pos,aa),o in zip(edit.items(),combo)],'search_score':float(score),'contacts':contacts,'predeclared_control':serial<len(controls),'previous_same_sequence_records':history.get(seq,[])})
  d['sidechain_sources']=[r for r in d.get('sidechain_sources',[]) if r.get('position') not in edit]+d['combination']['edits'];assert not re.search(r'N[^P][ST]',seq);assert len(seq)==125 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY');assert seq==''.join(r['aa'] for r in rr)
  (out/(cid+'.json')).write_text(json.dumps(d));report.append({'candidate_id':cid,'edits':edit,'status':'generated','is_new_sequence':seq not in history,'previous_same_sequence':history.get(seq,[]),'sequence':seq});print(cid,'new',seq not in history,'edits',edit,flush=True)
 (S/'reference/combination_generation.json').write_text(json.dumps({'records':report,'seconds':time.time()-start},indent=2));print('Generation finished',round(time.time()-start,2),flush=True)
if __name__=='__main__':main()
