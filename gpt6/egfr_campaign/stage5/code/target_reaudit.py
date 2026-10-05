"""Recheck intended target sequences, crop numbering and proximal glycan sites."""
from pathlib import Path
import sys,json
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage5';O=R/'intermediate'
sys.path[:0]=[str(R/'stage2/code'),str(R/'stage2/runtime/python'),str(R/'code')]
from mm_refine import app,u

def fasta(p):return ''.join(l.strip() for l in p.read_text().splitlines() if not l.startswith('>'))
def main():
 hu=fasta(O/'P00533.fasta');mo=fasta(O/'Q01279.fasta');records=[]
 for name in ['human6ARU','mouseAF','1IVO','1NQL','1YY9']:
  d=json.loads((O/(name+'_aligned.json')).read_text());seq=mo if name=='mouseAF' else hu;diff=[{'human_position':r['human_pos'],'uniprot_position':r['uniprot'],'reference':seq[r['uniprot']-1],'structure':r['aa']} for r in d if r['aa']!=seq[r['uniprot']-1]];crop=[r for r in d if 334<=r['human_pos']<=505];idx=[r['human_pos'] for r in crop];records.append({'target':name,'mismatches_to_uniprot':diff,'crop_count':len(crop),'crop_is_complete_contiguous':idx==list(range(334,506))})
 assert all(r['crop_is_complete_contiguous'] for r in records)
 assert not [d for r in records if r['target'] in ['human6ARU','mouseAF'] for d in r['mismatches_to_uniprot'] if 334<=d['human_position']<=505]
 (S/'reference/target_sequence_reaudit.json').write_text(json.dumps(records,indent=2))
 t={r['human_pos']:r for r in json.loads((O/'human6ARU_aligned.json').read_text())};p=app.PDBFile(str(S/'intermediate/refined/B00000_human6ARU_relaxed.pdb'));xyz=np.asarray(p.positions.value_in_unit(u.angstrom));b=xyz[[a.index for a in p.topology.atoms() if a.residue.chain.id=='B' and a.element!=app.element.hydrogen]];rows=[]
 for hp in sorted(set([540,634]+[i+1 for i in range(24,643) if hu[i]=='N' and hu[i+1]!='P' and hu[i+2] in 'ST'])):
  if hp not in t:continue
  c=np.array(list(t[hp]['atoms'].values()));rows.append({'position':hp,'wildtype_triplet':hu[hp-1:hp+2],'coordinate_residue':t[hp]['aa'],'nearest_binder_to_reference_residue_heavy_A':float(np.linalg.norm(b[:,None]-c[None,:],axis=2).min()),'nearest_binder_to_CA_A':float(np.linalg.norm(b-np.array(t[hp]['atoms']['CA']),axis=1).min())})
 (S/'reference/wildtype_glycan_site_proximity.json').write_text(json.dumps({'candidate':'B00000','source':'supplied human UniProt and6ARU aligned data; stage5 neutral refined binder','sites':rows,'qualification':'Proximity only. Actual forcefield crops match reference sequences exactly. Full6ARU structural sequence differs at540/634, both remote and outside the forcefield crop. Unresolved glycans are not invented.'},indent=2));print('Crops match intended human/mouse sequence exactly; all numbering contiguous172residues.')
if __name__=='__main__':main()
