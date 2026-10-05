"""Historical positive control: native FcRn/Fc contact histidines, no design seed.
Receptor A+B, interacting Fc chain C from 1I1A; noncontact Fc chain D omitted.
Only the <=5.5 A contact histidines are titrated; this is not a full FcRn pKa fit.
"""
from pathlib import Path
import sys,json
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3'
sys.path[:0]=[str(S3/'code'),str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import *
from mm_binding_probe import complete_target
import proton_ensemble as pe
source=R/'inputs/transfer2/egfr_refinement_fetch/1I1A.pdb';st=PDBParser(QUIET=True).get_structure('fcrn',str(source));out=S3/'intermediate/controls/fcrn';out.mkdir(parents=True,exist_ok=True)
mod=None;mapping={}
for label in ['A','B','C']:
    rows=[]
    for res in st[0][label]:
        if res.id[0]!=' ' or 'CA' not in res:continue
        rows.append({'aa':seq1(res.resname),'id':res.id[1],'atoms':{a.name:a.coord.tolist() for a in res if a.element!='H'}})
    mapping[label]={str(i+1):r['id'] for i,r in enumerate(rows)}
    rows=complete_target(rows);path=out/(label+'.pdb');write_pdb(rows,path,label);p=app.PDBFile(str(path))
    if mod is None:mod=app.Modeller(p.topology,p.positions)
    else:mod.add(p.topology,p.positions)
mod.topology.createDisulfideBonds(mod.positions)
path=out/'native_completed.pdb'
with path.open('w') as f:app.PDBFile.writeFile(mod.topology,mod.positions,f)
original_identify=pe.identify_sites
pe.identify_sites=lambda heavy,target_chains:original_identify(heavy,target_chains,distance_A=5.5)
result=pe.probe(path,out/'native_contact_control',('A','B'))
for h in result['sites']:h['author_residue_id']=mapping[h['chain']][h['pdb_residue_id']]
result['control']='1I1A native FcRn/Fc, contact-His subset; known acidic binding; not a sequence-design seed'
result['control_scope']={'target_chains':['A','B'],'binder_chains':['C'],'titration_distance_A':5.5,'glycans_omitted_from_energy':True,'other_His_fixed':'HIE','not_absolute_affinity_calibration':True}
(out/'native_contact_control_proton.json').write_text(json.dumps(result,indent=2))
