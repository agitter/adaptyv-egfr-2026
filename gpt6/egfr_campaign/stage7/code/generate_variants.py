"""Matched-backbone finite-rotamer Gly replacement experiment (1987 rationale).
Ala/Ser substitutions avoid positive-phi Gly and retain pH-contact identities.
Search minimizes worst local packing over human, mouse and unbound coordinates.
"""
from pathlib import Path
import sys,json,copy,itertools,hashlib,re,time
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';sys.path.insert(0,str(S/'code'))
from inspect_lead import load_model,sha
from classical_design import LIB,BB,conformer,atom_data,pair_score

def main():
 start=time.time();parentfile=R/'stage6/intermediate/designs/J00008.json';parent=json.loads(parentfile.read_text());paths={'human':R/'stage6/intermediate/refined/J00008_human6ARU_relaxed.pdb','mouse':R/'stage6/intermediate/refined/J00008_mouseAF_relaxed.pdb','free':R/'stage6/intermediate/unbound/J00008_free_relaxed.pdb'};models={k:load_model(p) for k,p in paths.items()}
 edits=[{}]+[{p:'A'} for p in [35,51,99,102,111,113]]+[{p:'S'} for p in [35,99,111,113]]+[{p:'A' for p in group} for group in [(35,111),(51,111),(102,111),(99,102,111),(35,51,111),(35,102,111),(51,102,111),(35,51,102,111),(35,51,99,102,111),(35,51,99,102,111,113)]]
 edits += [{35:'A',51:'A',102:'A',111:'S'},{35:'A',51:'A',102:'A',111:'A',113:'S'}]
 plan={'parent_file':str(parentfile),'parent_sha256':sha(parentfile),'model_sources':{k:{'path':str(p),'sha256':sha(p)} for k,p in paths.items()},'edits':edits,'methods':'Geometry-compatible Gly replacements; deterministic empirical rotamer enumeration; worst-case packing over human/mouse/free starting structures; unchanged-sequence matched-refinement control. Pre-2011 methods.','rationale_source':'https://doi.org/10.1073/pnas.84.19.6663','gates':{'positive_phi_Gly_preserved':True,'all_intended_pH_contact_identities_preserved':True,'new_cysteines_or_prolines':False,'new_glycosylation_sequons':False,'minimization_iterations_each':450,'relaxed_geometry_omega_tolerance_deg':25},'promotion_plan':'No claim of measured improvement. Require both-species strict geometry and retained acid-on direction. Expanded pH comparison for finalists. Thermal controls needed before claiming stress-test stabilization. Failures and unchanged parent retained.'}
 (S/'reference/gly_replacement_predeclaration.json').write_text(json.dumps(plan,indent=2))
 history={}
 for st in range(2,7):
  for p in (R/f'stage{st}/intermediate/designs').glob('*.json'):
   d=json.loads(p.read_text())
   if isinstance(d,dict) and isinstance(d.get('sequence'),str):history.setdefault(d['sequence'],[]).append(str(p.relative_to(R)))
 report=[]
 for ix,edit in enumerate(edits):
  cid=f'M{ix:05d}';options={};failed=None
  for pos,aa in edit.items():
   assert parent['sequence'][pos-1]=='G';assert parent['structure'][pos-1]['loop'];rank=[]
   for rn in range(len(LIB[aa]['xyz'])):
    scores=[]
    for key,(rows,target) in models.items():
     bb=np.asarray([rows[pos-1]['atoms'][n] for n in BB]);a=conformer(aa,rn,bb);dat=atom_data(aa,a,LIB[aa]['names']);env=np.concatenate([atom_data(r['aa'],r['atoms']) for j,r in enumerate(rows) if j!=pos-1]);iv=pair_score(dat,env);tv=pair_score(dat,np.concatenate([atom_data(r['aa'],r['atoms']) for r in target])) if target else (0.,0.,0.)
     scores.append({'model':key,'intra':list(iv),'inter':list(tv)})
    if max(r['intra'][0] for r in scores)>12 or max(r['inter'][0] for r in scores)>12:continue
    rank.append((max(sum(r['intra'])+sum(r['inter']) for r in scores),rn,scores))
   options[pos]=sorted(rank,key=lambda z:z[0])[:2]
   if not rank:failed='no rotamer passes predeclared steric cutoff';break
  if failed:report.append({'candidate_id':cid,'edits':edit,'status':failed});continue
  best=None
  for combo in itertools.product(*options.values()):
   total=max([0.]+[sum(c[0] for c in combo)]);bad=False
   for key,(rr,target) in models.items():
    atoms=[]
    for (pos,aa),(_,rn,sc) in zip(edit.items(),combo):
     a=conformer(aa,rn,np.asarray([rr[pos-1]['atoms'][n] for n in BB]));atoms.append(atom_data(aa,a,LIB[aa]['names']))
    for a,b in itertools.combinations(atoms,2):
     ps=pair_score(a,b);bad |= ps[0]>12;total+=sum(ps)
   if not bad and (best is None or total<best[0]):best=(total,combo)
  if best is None:report.append({'candidate_id':cid,'status':'sidechains clash as a set','edits':edit});continue
  score,combo=best;rows=copy.deepcopy(parent['structure'])
  for row,source in zip(rows,models['human'][0]):row['atoms']={n:np.asarray(v).tolist() for n,v in source['atoms'].items()}
  mutations=[]
  for (pos,aa),(_,rn,sc) in zip(edit.items(),combo):
   a=conformer(aa,rn,np.asarray([rows[pos-1]['atoms'][n] for n in BB]));rows[pos-1].update(aa=aa,atoms={n:v.tolist() for n,v in a.items()});mutations.append({'position':pos,'from':'G','to':aa,'rotamer':rn,'rotamer_source':LIB[aa]['sources'][rn],'multistructure_packing':sc})
  seq=''.join(r['aa'] for r in rows);assert set(seq)<=set('ACDEFGHIKLMNPQRSTVWY') and len(seq)==125;assert not re.search('N[^P][ST]',seq)
  d=copy.deepcopy(parent);d.update(candidate_id=cid,sequence=seq,structure=rows,cdr_sequences=[seq[a:b] for a,b in parent['cdr_intervals_zero_based']],campaign_stage=7,branch='geometry_compatible_gly_replacement',seed=None,metrics={'scores_valid':False,'note':'Inherited search/evaluation values invalid. All reported molecular energies recomputed.'},ancestry={'parent_file':str(parentfile),'parent_sha256':sha(parentfile),'parent_type':'own untested de novo computational design','backbone_file':str(paths['human']),'backbone_sha256':sha(paths['human'])},replacement_experiment={'mutations':mutations,'unchanged_control':not edit,'previous_same_sequence':history.get(seq,[]),'search_score_not_binding_energy':float(score)})
  d['sidechain_sources']=[r for r in d.get('sidechain_sources',[]) if r.get('position') not in edit]+mutations
  dest=S/'intermediate/designs'/f'{cid}.json';dest.write_text(json.dumps(d));report.append({'candidate_id':cid,'edits':edit,'status':'generated','sequence':seq,'new_sequence':seq not in history,'old_matches':history.get(seq,[])})
  print(cid,'new',seq not in history,'edit',edit,'search',round(score,3),flush=True)
 (S/'reference/generation.json').write_text(json.dumps({'records':report,'seconds':time.time()-start},indent=2));print('Finished',len(report),'records',flush=True)
if __name__=='__main__':main()
