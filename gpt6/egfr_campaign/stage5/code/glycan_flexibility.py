"""Geometric glycan torsion stress test using ONLY resolved supplied glycans.
Classical rigid rotations preserve covalent geometry. Protein/glycan clashes
filter the grid before any binder comparison. These are not Boltzmann-weighted
populations, a glycan forcefield, or claims about unresolved carbohydrate atoms.
"""
from pathlib import Path
import sys,json,math,itertools,hashlib
import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u

def rotation(x,origin,axis,degrees):
 axis=np.asarray(axis,float);axis/=np.linalg.norm(axis);v=np.asarray(x)-origin;a=math.radians(degrees)
 return origin+v*math.cos(a)+np.cross(axis,v)*math.sin(a)+np.outer(v@axis,axis)*(1-math.cos(a))

def main():
 gly=json.loads((R/'stage2/intermediate/target_glycans.json').read_text());target=json.loads((R/'intermediate/human6ARU_aligned.json').read_text());byid={r['id']:r for r in target};cif=MMCIF2Dict(str(R/'inputs/provided/6ARU.cif'));keys=[k for k in cif if k.startswith('_struct_conn.')];connections=[{k.split('.',1)[1]:cif[k][i] for k in keys} for i in range(len(cif[keys[0]]))];attachments=[]
 for c in connections:
  if c['conn_type_id']=='covale' and c['ptnr1_auth_asym_id']=='A' and c['ptnr1_auth_comp_id']=='ASN' and c['ptnr2_auth_comp_id']=='NAG':attachments.append((int(c['ptnr1_auth_seq_id']),c['ptnr2_auth_asym_id'],int(c['ptnr2_auth_seq_id'])))
 records=[];models={};tests=[]
 for pdbid,chain,rid in attachments:
  residue=byid[pdbid];hp=residue['human_pos'];indices=[i for i,g in enumerate(gly) if g['chain']==chain and (chain in ['D','E'] or g['resid']==rid)];assert indices
  gg=[gly[i] for i in indices];xyz=np.array([g['xyz'] for g in gg]);lookup={(g['resid'],g['atom']):i for i,g in enumerate(gg)};c1=lookup[(rid,'C1')];origin=np.array(residue['atoms']['ND2']);axis=xyz[c1]-origin
  prot=np.array([v for r in target for n,v in r['atoms'].items() if not(r['id']==pdbid and n=='ND2')]);other=np.array([g['xyz'] for i,g in enumerate(gly) if i not in indices]);obstacle=np.vstack([prot,other]);tree=cKDTree(obstacle)
  # Heavy-atom covalent graph inferred within each resolved sugar, plus CIF links.
  graph={i:set() for i in range(len(gg))}
  for i,j in itertools.combinations(range(len(gg)),2):
   if gg[i]['resid']==gg[j]['resid'] and np.linalg.norm(xyz[i]-xyz[j])<1.85:graph[i].add(j);graph[j].add(i)
  for c in connections:
   if c['conn_type_id']!='covale' or c['ptnr1_auth_asym_id']!=chain or c['ptnr2_auth_asym_id']!=chain:continue
   try:i=lookup[(int(c['ptnr1_auth_seq_id']),c['ptnr1_label_atom_id'])];j=lookup[(int(c['ptnr2_auth_seq_id']),c['ptnr2_label_atom_id'])]
   except (KeyError,ValueError):continue
   graph[i].add(j);graph[j].add(i)
  nonbonded=[]
  for i in range(len(gg)):
   reached={i};front={i}
   for _ in range(3):front={q for p in front for q in graph[p]}-reached;reached|=front
   nonbonded.extend((i,j) for j in range(i+1,len(gg)) if j not in reached)
  pairs=np.array(nonbonded,dtype=int).reshape(-1,2);bondpairs=np.array([(i,j) for i in graph for j in graph[i] if i<j],dtype=int);bondlength=np.linalg.norm(xyz[bondpairs[:,0]]-xyz[bondpairs[:,1]],axis=1)
  branch=np.array([i for i,g in enumerate(gg) if g['resid']>rid]);has_link=(rid,'O4') in lookup and (rid+1,'C1') in lookup and len(branch)>0;kept=[];grid=[]
  for phi in range(-60,61,15):
   for psi in (range(-45,46,15) if has_link else [0]):
    x=rotation(xyz,origin,axis,phi)
    if has_link:
     l=lookup[(rid,'O4')];r=lookup[(rid+1,'C1')];x[branch]=rotation(x[branch],x[l],x[r]-x[l],psi)
    outside=float(tree.query(x)[0].min());selfmin=float(np.linalg.norm(x[pairs[:,0]]-x[pairs[:,1]],axis=1).min()) if len(pairs) else 99.;bond_error=float(np.max(abs(np.linalg.norm(x[bondpairs[:,0]]-x[bondpairs[:,1]],axis=1)-bondlength)));root_error=abs(float(np.linalg.norm(x[c1]-origin)-np.linalg.norm(xyz[c1]-origin)));tests.append(bond_error<1e-8 and root_error<1e-8)
    ok=outside>=2.0 and selfmin>=2.0;row={'phi_perturbation_deg':phi,'first_internal_link_perturbation_deg':psi,'reference_environment_min_A':outside,'glycan_self_nonbonded_min_A':selfmin,'bond_length_max_error_A':bond_error,'root_bond_error_A':root_error,'steric_admissible':ok};grid.append(row)
    if ok:kept.append(x)
  zero=next(r for r in grid if r['phi_perturbation_deg']==r['first_internal_link_perturbation_deg']==0)
  # A failing native reference makes this grid uninterpretable; do not conceal it.
  records.append({'human_position':hp,'reference_PDB_position':pdbid,'chain':chain,'first_sugar_resid':rid,'resolved_atom_count':len(gg),'sampled_models':len(grid),'sterically_admissible_models':len(kept),'unperturbed_model_pass':zero['steric_admissible'],'grid':grid});models[hp]=np.array(kept)
 out=S/'intermediate/glycan_flexibility';out.mkdir(exist_ok=True);np.savez_compressed(out/'admissible_resolved_glycans.npz',**{str(k):v for k,v in models.items()});assert all(tests)
 binders=[]
 for path in sorted((S/'intermediate/refined').glob('*_human6ARU_relaxed.pdb')):
  cid=path.name.split('_')[0];p=app.PDBFile(str(path));x=np.asarray(p.positions.value_in_unit(u.angstrom));bx=x[[a.index for a in p.topology.atoms() if a.residue.chain.id=='B' and a.element!=app.element.hydrogen]];bt=cKDTree(bx);rr=[]
  for hp,mm in models.items():
   distances=[float(bt.query(g)[0].min()) for g in mm];rr.append({'human_glycan_position':hp,'admissible_models':len(mm),'models_with_binder_overlap_below2A':sum(d<2.0 for d in distances),'models_below2p5A':sum(d<2.5 for d in distances),'minimum_distances_A':distances,'minimum_over_grid_A':min(distances) if distances else None})
  binders.append({'candidate_id':cid,'source':str(path),'glycan_stress_test':rr})
 result={'source_cif':str(R/'inputs/provided/6ARU.cif'),'source_cif_sha256':hashlib.sha256((R/'inputs/provided/6ARU.cif').read_bytes()).hexdigest(),'reference_glycans':records,'bond_preservation_checks':len(tests),'all_bond_preservation_checks_pass':all(tests),'binder_comparisons':binders,'qualification':'Resolved glycan atoms only. Unweighted small torsion grid, not thermal populations or glycan binding free energies. Reference receptor is fixed. Clash-free grid points do not prove access, and a clash-prone grid does not establish complete blockade. Unresolved glycans/branches and induced fit are not modeled.'};(S/'reference/glycan_flexibility.json').write_text(json.dumps(result,indent=2))
 for r in records:print('GLYCAN',r['human_position'],'admissible',r['sterically_admissible_models'],'/',r['sampled_models'],'native',r['unperturbed_model_pass'])
 for r in binders:
  if r['candidate_id'] in ['B00000','N2070600','N2091100','T3070400','H00011','H00020','H00021']:print(r['candidate_id'],[(x['human_glycan_position'],x['models_with_binder_overlap_below2A'],x['admissible_models']) for x in r['glycan_stress_test']])
if __name__=='__main__':main()
