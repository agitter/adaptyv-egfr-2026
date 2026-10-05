"""Independent full-target/glycan/alternative-conformation and loop geometry checks.
Distances, signed volumes and Ramachandran regions are classical heuristics,
not a fold predictor. Do not merge reference and relaxed coordinate frames.
"""
from pathlib import Path
import sys,json,copy,math,hashlib
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u
from legacy_geometry import dihedral,rama_quality
from evaluate_pool import structural_audit,TARGETS

def run(path):
 cid=path.name.split('_')[0];d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());p=app.PDBFile(str(path));x=np.asarray(p.positions.value_in_unit(u.angstrom));binder=[];target={};species='mouseAF' if '_mouseAF_' in path.name else 'human6ARU'
 for res in p.topology.residues():
  a={a.name:x[a.index].tolist() for a in res.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
  if res.chain.id=='B':binder.append(a)
  if res.chain.id=='A':target[int(res.id)+333]=a
 assert len(binder)==len(d['structure'])
 for row,a in zip(d['structure'],binder):row['atoms']=a
 baseline=structural_audit(d);whole=np.array([p for r in TARGETS[species] for p in target.get(r['human_pos'],r['atoms']).values()]);bx=np.array([p for row in binder for p in row.values()]);dist=cKDTree(whole).query(bx)[0]
 angles=[]
 for i,row in enumerate(d['structure']):
  if not row['loop']:continue
  a=binder[i];ph=dihedral(np.asarray(binder[i-1]['C']),np.asarray(a['N']),np.asarray(a['CA']),np.asarray(a['C']));ps=dihedral(np.asarray(a['N']),np.asarray(a['CA']),np.asarray(a['C']),np.asarray(binder[i+1]['N']));q,pos=rama_quality(np.array([[ph,ps]]));phi=math.degrees(ph);psi=math.degrees(ps)
  angles.append({'position':i+1,'aa':row['aa'],'phi_deg':phi,'psi_deg':psi,'generic_region_distance':float(q[0]),'nonG_positive_phi':row['aa']!='G' and phi>0,'proline_outside_phi_window':row['aa']=='P' and not -100<phi<-35})
 # Positive phi is a diagnostic, not automatic rejection of all non-glycine.
 geom=json.loads(path.with_name(path.name.replace('_relaxed.pdb','_independent.json')).read_text())
 same={'species':species,'minimum_heavy_distance_A':float(dist.min()),'heavy_atoms_below2A':int(sum(dist<2.))}
 gate=geom['pass'] and same['heavy_atoms_below2A']==0 and baseline['resolved_glycans']['atoms_below_2A']==0 and baseline['1IVO']['heavy_atoms_below_2A']==0 and baseline['1NQL']['heavy_atoms_below_2A']==0
 r={'candidate_id':cid,'species':species,'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'same_species_full_context_with_relaxed_crop':same,'static_reference_comparisons':baseline,'loop_torsions':angles,'independent_geometry_pass':geom['pass'],'combined_geometry_context_pass':gate,'qualification':'Limited to supplied conformations and resolved glycans; no folding or binding guarantee.'}
 out=S/'intermediate/context';out.mkdir(exist_ok=True);(out/(path.stem+'.json')).write_text(json.dumps(r,indent=2));print(cid,species,gate,'same-clash',same['heavy_atoms_below2A'],'alt',[baseline[k]['heavy_atoms_below_2A'] for k in ['1IVO','1NQL']],'gly',baseline['resolved_glycans']['atoms_below_2A'],flush=True);return r
if __name__=='__main__':
 records=[run(p) for p in sorted((S/'intermediate/refined').glob('*_relaxed.pdb')) if p.with_name(p.name.replace('_relaxed.pdb','_independent.json')).exists()]
 (S/'reference/context_summary.json').write_text(json.dumps(records,indent=2))
