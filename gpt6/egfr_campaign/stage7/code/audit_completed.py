"""Audit completed structures uniformly; retain failures and original scores."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/data/egfr_campaign');S=R/'stage7';sys.path[:0]=[str(R/x/'code') for x in ['stage6','stage5','stage4','stage3','stage2']]+[str(R/'code'),str(R/'stage2/runtime/python')]
from global_audit import model_audit,sha
from mm_refine import app,u

def glycan(path):
 p=app.PDBFile(str(path));x=np.asarray(p.positions.value_in_unit(u.angstrom));bx=x[[a.index for a in p.topology.atoms() if a.residue.chain.id=='B' and a.element!=app.element.hydrogen]];tree=cKDTree(bx);f=R/'stage5/intermediate/glycan_flexibility/admissible_resolved_glycans.npz';data=np.load(f);rows=[]
 for key in data.files:
  distances=[float(tree.query(g)[0].min()) for g in data[key]];n=sum(v<2 for v in distances);rows.append({'root':key,'count':len(distances),'overlaps':n,'fraction':n/len(distances)})
 return {'source':str(f),'source_sha256':sha(f),'by_root':rows,'maximum_overlap_fraction':max(r['fraction'] for r in rows),'pass':max(r['fraction'] for r in rows)<=.5,'qualification':'Unweighted finite glycan conformer stress test, not a probability.'}

def main():
 out=S/'intermediate/audits';out.mkdir(exist_ok=True);rows=[]
 for f in sorted((S/'intermediate/refined').glob('*_proton.json')):
  cid=f.name.split('_')[0];sp='mouse' if 'mouseAF' in f.name else 'human';dpath=S/'intermediate/designs'/f'{cid}.json';d=json.loads(dpath.read_text());p=json.loads(f.read_text());dest=out/f.name
  if dest.exists():a=json.loads(dest.read_text())
  else:
   rec={'key':'stage7:'+cid,'id':cid,'stage':'stage7','species':sp,'model':str(f),'coordinate_source':p['source'],'design_source':str(dpath),'sequence':d['sequence'],'sequence_sha256':hashlib.sha256(d['sequence'].encode()).hexdigest()};a=model_audit(rec);a['glycan_grid']=glycan(Path(p['source']));a['combined_context_and_glycan_pass']=a['combined_geometry_pass'] and a['glycan_grid']['pass'];dest.write_text(json.dumps(a,indent=2))
  narrow=a['independent_priors'];rows.append({'candidate_id':cid,'species':sp,'all_geometry_pass':a['combined_context_and_glycan_pass'],'strict_geometry_pass':a['strict_geometry']['pass'],'peptide_warnings':a['strict_geometry']['peptide_omega_warnings'],'human_alternative_clashes':{k:a['reference_geometry'][k]['heavy_atoms_below_2A'] for k in ['1IVO','1NQL']},'minimum_histidine_contrast':narrow.get('minimum_contrast_kcal'),'maximum_histidine_contrast':narrow.get('maximum_contrast_kcal'),'worst_case':narrow.get('worst_case'),'glycan_fraction':a['glycan_grid']['maximum_overlap_fraction'],'source':str(dest)})
 (S/'reference/structure_summary.json').write_text(json.dumps(rows,indent=2));print('Audited',len(rows),'passing',sum(r['all_geometry_pass'] for r in rows),flush=True)
 for row in rows:print(row['candidate_id'],row['species'],row['all_geometry_pass'],round(row['minimum_histidine_contrast'],4),row['glycan_fraction'],flush=True)
if __name__=='__main__':main()
