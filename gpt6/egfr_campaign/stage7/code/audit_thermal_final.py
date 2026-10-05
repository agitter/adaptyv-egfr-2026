"""Apply the existing independent peptide/chirality/disulfide gate to endpoints.
This does not inspect intermediate heavy-atom geometries: only CA trajectories
and final all-atom structures were saved, so that limitation remains explicit.
"""
from pathlib import Path
import sys,json,hashlib
R=Path('/mnt/data/egfr_campaign');S=R/'stage7'
sys.path[:0]=[str(R/'stage4/code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from audit_refined import audit
from mm_refine import app
from Bio.SeqUtils import seq1

def main():
 result=json.loads((S/'reference/thermal_analysis.json').read_text());rows=[]
 for r in result['replicas']:
  p=Path(r['source_result']);dest=p.with_name(p.stem+'_final.pdb');pp=app.PDBFile(str(dest));chains=list(pp.topology.chains());assert len(chains)==1
  residues=list(pp.topology.residues());sequence=''.join(seq1(x.name,custom_map={'HIE':'H','HID':'H','HIP':'H'}) for x in residues)
  cid=r['candidate_id']
  if cid=='NATIVE_3EAK':d={'candidate_id':cid,'sequence':sequence};assert len(sequence)==127
  else:
   d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());assert sequence==d['sequence']
  a=audit(dest,d,chain=chains[0].id)
  rows.append({'candidate_id':cid,'seed':r['seed'],'added_chirality_restraints':r['added_chirality_restraints'],'result_source':str(p),'final_structure':str(dest),'final_structure_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'independent_audit':a})
 out={'records':rows,'strict_pass_count':sum(r['independent_audit']['pass'] for r in rows),'completed_endpoints':len(rows),'scope':'Final all-atom endpoints only; no all-atom claim for intermediate CA-only saved frames.'}
 (S/'reference/thermal_final_geometry.json').write_text(json.dumps(out,indent=2));print('Independent final geometry',out['strict_pass_count'],'of',len(rows))
 for r in rows:
  a=r['independent_audit'];print(r['candidate_id'],r['seed'],'extra',r['added_chirality_restraints'],'pass',a['pass'],'peptide_warnings',a['peptide_omega_warnings'])
if __name__=='__main__':main()
