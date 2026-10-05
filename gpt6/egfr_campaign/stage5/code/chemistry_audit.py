"""Rule, ancestry and geometry gates for explicit own-design loop edits.
No protected-acid identity assertion: these are deliberately changed and every
change must instead agree with the recorded edit ledger. No scoring inheritance.
"""
from pathlib import Path
import sys,json,re,hashlib,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S=R/'stage5'
sys.path[:0]=[str(S/'code'),str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from build_return_loops import load_parent,sha
from evaluate_pool import sequence_audit,structural_audit,contact_audit
from legacy_geometry import dihedral

def main():
 parent,original,tx=load_parent();seen=set();records=[];out=S/'intermediate/evaluation';out.mkdir(exist_ok=True)
 for p in sorted((S/'intermediate/designs').glob('H*.json')):
  d=json.loads(p.read_text());seq=d['sequence'];sa=sequence_audit(d);st=structural_audit(d);rows=d['structure'];edits=d['contact_redesign']['edits'];expected={e['parent_position']:(e['from'],e['to']) for e in edits};actual={i+1:(a['aa'],b['aa']) for i,(a,b) in enumerate(zip(original,rows)) if a['aa']!=b['aa']};assert len(rows)==len(original)
  identity=seq==''.join(r['aa'] for r in rows) and d['cdr_sequences']==[seq[a:b] for a,b in d['cdr_intervals_zero_based']]
  hashes=sha(d['ancestry']['parent_file'])==d['ancestry']['parent_sha256'] and sha(d['ancestry']['relaxed_backbone_source'])==d['ancestry']['relaxed_backbone_sha256']
  framework=all(a['aa']==b['aa'] and a['atoms']==b['atoms'] for a,b in zip(original,rows) if not a['loop'])
  backbone=all(all(np.allclose(a['atoms'][n],b['atoms'][n],rtol=0,atol=1e-10) for n in ['N','CA','C','O']) for a,b in zip(original,rows))
  loop_edits=all(original[i-1]['loop'] for i in actual);motif=re.findall(r'C([ACDEFGHIKLMNPQRSTVWY]{5,40}?)[WF][GAS]QG',seq)
  omega=[];inversions=[]
  for i,r in enumerate(rows):
   a={n:np.array(v) for n,v in r['atoms'].items()}
   if 'CB' in a and np.linalg.det(np.array([a[n]-a['CA'] for n in ['N','C','CB']]))<=0:inversions.append(i+1)
   if i:
    prev=rows[i-1]['atoms'];v=math.degrees(dihedral(np.array(prev['CA']),np.array(prev['C']),a['N'],a['CA']))
    if abs(180-abs(v))>25:omega.append({'previous_position':i,'omega_deg':v})
  duplicate=seq in seen;seen.add(seq);passed=bool(sa['format_pass'] and identity and hashes and framework and backbone and loop_edits and expected==actual and motif==[d['cdr_sequences'][2]] and not duplicate and not inversions and not omega and st['geometry_prefilter_pass'] and not sa['liabilities']['glycosylation_sequons'] and not sa['liabilities']['hydrophobic_runs'] and d['metrics']['scores_valid'] is False)
  rec={'candidate_id':d['candidate_id'],'sequence_audit':sa,'structure_audit':st,'edits_exactly_match_ledger':expected==actual,'edits_in_newly_generated_loops_only':loop_edits,'sequence_coordinates_match':identity,'parent_hashes_match':hashes,'framework_unchanged':framework,'backbone_unchanged':backbone,'motif_CDR3_matches_design':motif==[d['cdr_sequences'][2]],'duplicate_sequence':duplicate,'absolute_alpha_inversions':inversions,'peptide_omega_warnings':omega,'pass':passed,'contact_audit':contact_audit(d)}
  (out/p.name).write_text(json.dumps(rec,indent=2));records.append({'candidate_id':d['candidate_id'],'pass':passed,'edits':edits,'length':len(seq),'sequence':seq,'glycan_min_A':st['resolved_glycans']['min_heavy_distance_A'],'geometry_pass':st['geometry_prefilter_pass']});print(d['candidate_id'],passed,flush=True)
 (S/'reference/chemistry_design_audit_summary.json').write_text(json.dumps(records,indent=2))
 # Chemically motivated representatives chosen before any human energy outcomes.
 wanted=['H00001','H00005','H00007','H00011','H00016','H00018','H00020','H00021'];selected=[r for r in records if r['candidate_id'] in wanted and r['pass']]
 pre={'rule':'Mechanistic coverage, not cherry-picked MM: longer D54E; alternative H383 arm; acid-removal/paired-arm variants; E106D contact negative control; two continuously optimized E104 proposals. All geometry/rule gates first. Same B00000 matched comparator.','selected':selected,'control':'B00000','selection_before_all_atom':True,'continuous_proposal_extension':['H00020','H00021'],'neutral_pH_nondetection_not_established':True}
 (S/'reference/chemistry_human_selection.json').write_text(json.dumps(pre,indent=2));print('TOTAL',len(records),'PASS',sum(r['pass'] for r in records),'SELECTED',[r['candidate_id'] for r in selected],flush=True)
if __name__=='__main__':main()
