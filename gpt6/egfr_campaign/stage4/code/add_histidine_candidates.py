"""Targeted, explicitly falsifiable F107H probe on two own-design backgrounds.
A histidine neutral tautomer can already donate here; no acid-on claim is made.
"""
from pathlib import Path
import sys,json,copy,hashlib
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from classical_design import LIB,BB,conformer
from mm_refine import app,u

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 opportunities=json.loads((S/'reference/histidine_opportunities.json').read_text())['options'];made=[]
 for cid,parent,root in [('S00062','C00003',S/'intermediate/matched_controls'),('S00063','S00028',S)]:
  src=root/'intermediate/designs'/(parent+'.json');pfile=root/'intermediate/refined'/(parent+'_human6ARU_relaxed.pdb');d=json.loads(src.read_text());p=app.PDBFile(str(pfile));x=np.asarray(p.positions.value_in_unit(u.angstrom));res=[r for r in p.topology.residues() if r.chain.id=='B'];rows=copy.deepcopy(d['structure'])
  for row,r in zip(rows,res):row['atoms']={a.name:x[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
  option=next(z for z in opportunities if z['parent']==parent and z['position']==107);k=option['rotamer'];a=conformer('H',k,np.array([rows[106]['atoms'][n] for n in BB]));rows[106]['aa']='H';rows[106]['atoms']={n:v.tolist() for n,v in a.items()};seq=list(d['sequence']);assert seq[106]=='F';seq[106]='H';seq=''.join(seq);d.update(candidate_id=cid,sequence=seq,structure=rows,branch='additional_proton_uptake_falsification_probe',campaign_stage=4,status='untested computational candidate');d['cdr_sequences']=[seq[a:b] for a,b in d['cdr_intervals_zero_based']];d['metrics']={'scores_valid':False,'note':'Parent scores are not current results. All neutral His tautomers must be explicitly tested.'};d['ancestry']={'parent_candidate':parent,'parent_type':'own untested de novo design; equal-refinement control when applicable','parent_file':str(src),'parent_sha256':sha(src),'relaxed_backbone_source':str(pfile),'relaxed_backbone_sha256':sha(pfile)};d['charge_addition_probe']={'edit':'F107H','proposed_target_carboxylate_human_position':344,'geometry':option,'interpretation':'Only one histidine donor has a geometric partner; neutral-tautomer bypass is a known risk, not ignored.'};d['sidechain_sources']=[z for z in d['sidechain_sources'] if z['position']!=107]+[{'position':107,'aa':'H','rotamer':k,'source':LIB['H']['sources'][k]}];out=S/'intermediate/designs'/(cid+'.json');assert not out.exists();out.write_text(json.dumps(d,indent=2));made.append({'candidate_id':cid,'sequence':seq,'parent':parent,'out_sha256':sha(out)});print(cid,parent,'F107H',flush=True)
 (S/'reference/additional_proton_probe_designs.json').write_text(json.dumps(made,indent=2))
if __name__=='__main__':main()
