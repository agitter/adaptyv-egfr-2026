"""Reuse audited classical interface refinement; independently inspect cystines."""
from pathlib import Path
import sys,json,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
import refine_new
from mm_refine import app,u
from legacy_geometry import dihedral

def run(cid,species='human6ARU'):
    refine_new.S3=S;result=refine_new.run(cid,species,450);print(result,flush=True)
    path=S/'intermediate/refined'/f'{cid}_{species}_relaxed.pdb'
    if not path.exists():return
    d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());p=app.PDBFile(str(path));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));rows=[{a.name:xyz[a.index] for a in r.atoms()} for r in p.topology.residues() if r.chain.id=='B'];inversions=[];omega=[];cystines=[]
    for i,r in enumerate(rows):
        if 'CB' in r and np.linalg.det(np.array([r[n]-r['CA'] for n in ['N','C','CB']]))<=0:inversions.append(i+1)
        if i:
            theta=math.degrees(dihedral(rows[i-1]['CA'],rows[i-1]['C'],r['N'],r['CA']));deviation=abs(180-abs(theta))
            if deviation>25:omega.append({'previous_position':i,'next_position':i+1,'omega_deg':theta,'next_is_proline':d['sequence'][i]=='P'})
    for a,b in p.topology.bonds():
        if a.name==b.name=='SG' and a.residue.chain.id==b.residue.chain.id=='B':
            i,j=int(a.residue.id),int(b.residue.id);ai,aj=rows[i-1],rows[j-1];cystines.append({'positions':sorted([i,j]),'SG_distance_A':float(np.linalg.norm(ai['SG']-aj['SG'])),'chi3_deg':math.degrees(dihedral(ai['CB'],ai['SG'],aj['SG'],aj['CB']))})
    intended=d.get('intended_extra_disulfide_pairs',[]);present={tuple(x['positions']) for x in cystines};extra_ok=all(tuple(sorted(x)) in present for x in intended)
    audit={'candidate_id':cid,'species':species,'absolute_alpha_inversions':inversions,'peptide_omega_warnings':omega,'actual_cystines':cystines,'intended_extra_pairs_present':extra_ok,'pass':not inversions and not omega and extra_ok}
    (S/'intermediate/refined'/f'{cid}_{species}_independent.json').write_text(json.dumps(audit,indent=2));print('INDEPENDENT',cid,audit,flush=True)
if __name__=='__main__':run(*sys.argv[1:])
