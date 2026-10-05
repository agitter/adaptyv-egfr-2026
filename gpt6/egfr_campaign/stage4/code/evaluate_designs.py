"""Independent audits; paired extra cysteines are not a format violation."""
from pathlib import Path
import sys,json,re,hashlib,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage4'
sys.path[:0]=[str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from evaluate_pool import structural_audit,contact_audit,sequence_audit
from legacy_geometry import dihedral

def main():
    dest=S/'intermediate/evaluation';dest.mkdir(parents=True,exist_ok=True);summary=[]
    for p in sorted((S/'intermediate/designs').glob('*.json')):
        d=json.loads(p.read_text());seq=d['sequence'];a=sequence_audit(d);a['legacy_no_extra_Cys_filter']=a.pop('format_pass');a['format_pass']=10<=len(seq)<=250 and set(seq)<=set('ACDEFGHIKLMNPQRSTVWY');a['liabilities'].pop('unpaired_cdr_cysteines',None)
        extra=d.get('intended_extra_disulfide_pairs',[]);native=[i+1 for i,c in enumerate(seq) if c=='C' and all(i+1 not in pair for pair in extra)];coverage=len(native)==2 and all(seq[i-1]=='C' and seq[j-1]=='C' for i,j in extra) and len(set(x for pair in extra for x in pair))==2*len(extra)
        a['intended_cysteine_pairing_coverage_pass']=coverage;a['native_pair_positions']=native;a['intended_extra_pairs']=extra
        g=structural_audit(d);identity=seq==''.join(r['aa'] for r in d['structure']) and d['cdr_sequences']==[seq[x:y] for x,y in d['cdr_intervals_zero_based']]
        pro=[]
        for i,row in enumerate(d['structure']):
            if row['aa']=='P':
                at=row['atoms'];pro.append({'position':i+1,'N_CD_distance_A':float(np.linalg.norm(np.array(at['N'])-at['CD']))})
        parent=Path(d['ancestry']['parent_file']);phash=hashlib.sha256(parent.read_bytes()).hexdigest()==d['ancestry']['parent_sha256']
        data={'candidate_id':d['candidate_id'],'sequence_audit':a,'structure_audit':g,'coordinate_sequence_match':identity,'parent_hash_match':phash,'proline_rings':pro,'candidate_ion_pairs':contact_audit(d),'audit_pass':bool(a['format_pass'] and coverage and identity and phash and g['geometry_prefilter_pass'] and not a['liabilities']['glycosylation_sequons'])}
        (dest/p.name).write_text(json.dumps(data,indent=2));summary.append({'candidate_id':d['candidate_id'],'audit_pass':data['audit_pass'],'extra_disulfides':len(extra),'gating_contacts':len(data['candidate_ion_pairs'])})
    (S/'reference/design_audit_summary.json').write_text(json.dumps(summary,indent=2));print('AUDITED',len(summary),'PASS',sum(r['audit_pass'] for r in summary),'CYS',sum(r['extra_disulfides'] for r in summary))
if __name__=='__main__':main()
