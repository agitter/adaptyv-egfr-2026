"""Apply the same archived, unweighted glycan torsion grid to the candidate pool.
This is a geometric stress test, not a population or binding-probability estimate.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u
from select_panel import make_pool
source=R/'stage5/intermediate/glycan_flexibility/admissible_resolved_glycans.npz';grids=np.load(source);rows,_=make_pool(apply_glycan_filter=False);results=[];cache={}
for row in rows:
 paths=[m['coordinate_source'] for m in row['human_models']];models=[]
 if not paths:paths=[row['design_source']]
 for path in paths:
  if path in cache:r=cache[path]
  else:
   if path.endswith('.pdb'):
    p=app.PDBFile(path);x=np.asarray(p.positions.value_in_unit(u.angstrom));bx=x[[a.index for a in p.topology.atoms() if a.residue.chain.id=='B' and a.element!=app.element.hydrogen]]
   else:d=json.loads(Path(path).read_text());bx=np.asarray([v for res in d['structure'] for v in res['atoms'].values()])
   tree=cKDTree(bx);byroot=[]
   for key in grids.files:
    distances=[float(tree.query(g)[0].min()) for g in grids[key]];clashes=sum(v<2 for v in distances);byroot.append({'root_human_position':int(key),'admissible_grid_models':len(distances),'overlapping_grid_models':clashes,'overlap_fraction_unweighted':clashes/len(distances),'minimum_distances_A':distances})
   r={'source':path,'source_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'by_root':byroot,'maximum_unweighted_overlap_fraction':max(z['overlap_fraction_unweighted'] for z in byroot)};cache[path]=r
  models.append(r)
 worst=max(r['maximum_unweighted_overlap_fraction'] for r in models);results.append({'sequence_sha256':row['sequence_sha256'],'candidate_key':row['candidate_key'],'models':models,'worst_unweighted_overlap_fraction':worst,'conservative_panel_filter_pass':worst<=.5})
report={'grid_source':str(source),'grid_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'maximum_overlap_fraction_for_panel':.5,'qualification':'A conservative panel-selection heuristic on an unweighted finite grid. Fractions are not physical glycan populations or binding probabilities. Other glycans and induced fit are not modeled.','records':results};(S/'reference/glycan_pool_screen.json').write_text(json.dumps(report,indent=2));print('Glycan grid screened',len(results),'retained',sum(r['conservative_panel_filter_pass'] for r in results))
for r in results:
 if r['candidate_key'].startswith('stage6:J') or r['candidate_key'] in ['stage5:H00011','stage5:H00021','stage5:N2091100','stage5:N2070600']:print(r['candidate_key'],r['worst_unweighted_overlap_fraction'],r['conservative_panel_filter_pass'])
