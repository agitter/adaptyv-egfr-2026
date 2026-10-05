"""Independent final coordinates in full and alternate target environments.
Classical distances and inherited least-squares-aligned target input data only.
"""
from pathlib import Path
import sys,json,copy,hashlib
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from mm_refine import app,u
from evaluate_pool import structural_audit,TARGETS

def run(root,path):
 cid=path.name.split('_')[0];dp=root/'intermediate/designs'/(cid+'.json');d=json.loads(dp.read_text());p=app.PDBFile(str(path));x=np.asarray(p.positions.value_in_unit(u.angstrom));binder=[];relaxed={}
 for r in p.topology.residues():
  a={a.name:x[a.index].tolist() for a in r.atoms() if a.element!=app.element.hydrogen and a.name!='OXT'}
  if r.chain.id=='B':binder.append(a)
  elif r.chain.id=='A':relaxed[int(r.id)+333]=a
 assert len(binder)==len(d['structure'])
 for r,a in zip(d['structure'],binder):r['atoms']=a
 baseline=structural_audit(d);species='mouseAF' if '_mouseAF_' in path.name else 'human6ARU';whole=[]
 for r in TARGETS[species]:whole.extend(relaxed.get(r['human_pos'],r['atoms']).values())
 bx=np.array([a for r in binder for a in r.values()]);dist=cKDTree(whole).query(bx)[0]
 record={'candidate_id':cid,'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'static_reference_coordinate_comparisons':baseline,'same_species_full_context_with_refined_crop':{'species':species,'minimum_heavy_distance_A':float(dist.min()),'heavy_atoms_below2A':int(sum(dist<2.))},'qualification':'Severe overlap filter only. Alternative receptor coordinates remain fixed; no conformational ensemble, glycans beyond resolved inputs, or induced-fit claim.'}
 out=S/'intermediate/final_context';out.mkdir(exist_ok=True);tag=path.stem if root==S else root.name+'_'+path.stem;(out/(tag+'.json')).write_text(json.dumps(record,indent=2));return record

def main():
 records=[]
 roots=[S,S/'intermediate/matched_controls']+[S/'intermediate/histidine_matched_controls'/q for q in ['S00062','S00063']]
 for root in roots:
  for p in sorted((root/'intermediate/refined').glob('*_relaxed.pdb')):
   rec=run(root,p);records.append(rec);b=rec['static_reference_coordinate_comparisons'];print(root.name,p.stem,'same-species',rec['same_species_full_context_with_refined_crop']['heavy_atoms_below2A'],'1IVO',b['1IVO']['heavy_atoms_below_2A'],'1NQL',b['1NQL']['heavy_atoms_below_2A'],'gly',b['resolved_glycans']['atoms_below_2A'],flush=True)
 (S/'reference/final_context_summary.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
