"""Summarize real completed outputs; do not promote hypotheses to binders."""
from pathlib import Path
import json, collections, math
import numpy as np
S=Path('/mnt/data/egfr_campaign/stage2'); R=S.parent

def load(p): return json.loads(p.read_text())
def write(p,d): p.write_text(json.dumps(d,indent=2))

def main():
    designs=[load(p) for p in sorted((S/'intermediate/designs').glob('P*.json')) if not p.stem.endswith(('_summary','_failed'))]
    evaluations={p.stem:load(p) for p in (S/'intermediate/evaluation').glob('*.json')}
    novelty=load(S/'intermediate/novelty/candidate_cdr3_screen_augmented.json')
    print('novelty_schema',type(novelty).__name__)
    good=[e for e in evaluations.values() if e['structure_audit']['geometry_prefilter_pass']]
    robust=[e for e in good if e.get('human_proton_model',{}).get('all_scenarios_acid_on',False)]
    summary={'status':'CHECKPOINT - NOT A SUBMISSION PACKAGE','date_utc':'2026-10-01',
      'design_count':len(designs),'unique_sequences':len({d['sequence'] for d in designs}),
      'unique_cdr3':len({d['cdr_sequences'][2] for d in designs}),
      'poses_with_sequences':len({d['pose_id'] for d in designs}),
      'length_range':[min(len(d['sequence']) for d in designs),max(len(d['sequence']) for d in designs)],
      'geometry_prefilter_pass':len(good),'geometry_prefilter_rejected':len(designs)-len(good),
      'coarse_acid_on_all_12_scenarios':len(robust),'geometry_pass_but_direction_not_robust':len(good)-len(robust),
      'candidates_with_no_glycosylation_sequon_and_robust_coarse_direction':sum(not e['sequence_audit']['liabilities']['glycosylation_sequons'] for e in robust),
      'geometry_failure_counts_overlapping':{n:sum(e['structure_audit'][n]['heavy_atoms_below_2A']>0 for e in evaluations.values()) for n in ['human6ARU','mouseAF']},
      'alternate_target_conformation_zero_clash_counts_among_geometry_pass':{n:sum(e['structure_audit'][n]['heavy_atoms_below_2A']==0 for e in good) for n in ['1IVO','1NQL']},
      'interpretation':'Positive surrogate pH contrast is only a direction check. It does not establish no detectable binding at pH 7.4, affinity, expression, or cross-species binding.',
      'third_phone_a_friend_used':False}
    summary['geometry_failure_counts_overlapping']['resolved_target_glycans']=sum(e['structure_audit']['resolved_glycans']['atoms_below_2A']>0 for e in evaluations.values())
    summary['novelty']=load(S/'reference/novelty_screen_augmented_stats.json')
    summary['novelty_reference']=load(S/'reference/novelty_augmentation.json')
    summary['tests']={n:load(S/'reference'/f'{n}.json')['passed'] for n in ['classical_tests','proton_polynomial_tests']}
    summary['allatom_joint_audit']=load(S/'intermediate/mm_binding/P00024_07_human6ARU_complex_audit.json')
    probe=load(S/'intermediate/mm_binding/P00024_07_human6ARU_joint_probe.json')
    summary['joint_mmgbsa_diagnostic']={'candidate_id':probe['candidate_id'],'neutral_reference_kcal_proxy':probe['neutral_reference']['delta_kcal_proxy'],
        'protonation_shifts':[{'partner':a['partner'],'position':a['residue_position'],'delta_delta_binding_kcal_proxy':a['delta_delta_binding_kcal_proxy']} for a in probe['single_states']],
        'interpretation':'Positive matched-component protonation shift opposes acid-on binding in this diagnostic. This is not a binding free energy or calibrated affinity prediction. Joint relaxation corrects large steric artifacts of the earlier monomer-only diagnostic.'}
    ht=load(R/'intermediate/human6ARU_aligned.json')
    summary['receptor_histidine_opportunities']=[{'human_uniprot_position':a['human_pos'],'human_aa':a['aa'],'mouse_aa':a['mouse_aa'],'sasa_A2':a['sasa']} for a in ht if a['aa']=='H' and 334<=a['human_pos']<=505]
    summary['conformation_warning']='All 211 current human/mouse geometry survivors collide with the ligand-bound 1IVO domain-I arrangement. This is a state-selectivity limitation, not evidence of active-state binding. 154 do not collide with the unliganded 1NQL conformation.'
    table=[]
    for d in designs:
        e=evaluations[d['candidate_id']];m=e.get('human_proton_model',{});q=e['sequence_audit'];g=e['structure_audit']
        table.append({'id':d['candidate_id'],'pose':d['pose_id'],'sequence':d['sequence'],'cdrs':d['cdr_sequences'],'length':len(d['sequence']),
          'geometry_prefilter_pass':g['geometry_prefilter_pass'],'coarse_all_scenarios_acid_on':m.get('all_scenarios_acid_on',False),
          'coarse_worst_contrast_kcal':m.get('worst_contrast_kcal_surrogate'),
          'liabilities':q['liabilities'],'alternate_1NQL_clashing_atoms':g['1NQL']['heavy_atoms_below_2A'],'alternate_1IVO_clashing_atoms':g['1IVO']['heavy_atoms_below_2A'],
          'status':'unvalidated intermediate; do not submit'})
    write(S/'output/intermediate_candidate_inventory.json',table)
    write(S/'reference/stage2_summary.json',summary)
    # Archive-only library. This is not the requested final ranked 100-sequence file.
    with (S/'output/intermediate_366_NOT_FOR_SUBMISSION.fasta').open('w') as f:
        for d in designs:
            f.write('>'+d['candidate_id']+' status=UNVALIDATED_INTERMEDIATE\n')
            for start in range(0,len(d['sequence']),80):f.write(d['sequence'][start:start+80]+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['novelty','allatom_joint_audit','receptor_histidine_opportunities']},indent=2))
if __name__=='__main__':main()
