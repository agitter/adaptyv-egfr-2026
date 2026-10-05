"""Reconstruct checkpoint counts and comparisons from actual retained outputs.
No protein design or scientific scoring is introduced by this reporting audit.
"""
from pathlib import Path
from collections import defaultdict
import json,hashlib,math,datetime,re
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
from Bio.SeqUtils import seq1
R=Path('/mnt/data/egfr_campaign');S=R/'stage3';OUT=S/'output';OUT.mkdir(exist_ok=True)

def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ensemble(d):
    scenarios=d['proton_model']['scenarios']
    key='assumed_unbound_pKa' if 'assumed_unbound_pKa' in scenarios[0] else 'assumed_unbound_histidine_pKa'
    chosen=[s for s in scenarios if s['solute_dielectric']==1 and s['kappa_nm_inverse']==1.25 and s[key]==6.3 and s['assumed_unbound_HIE_fraction']==.5 and s.get('assumed_acid_pKa_offset',0)==0]
    contrasts=[s['contrast_kcal'] for s in scenarios]
    out={'microstates':d['microstate_count'],'scenario_count':len(scenarios),'minimum_contrast_kcal':min(contrasts),'maximum_contrast_kcal':max(contrasts),'all_scenarios_acid_on':min(contrasts)>0,'sites':d['sites']}
    if chosen:
        assert len(chosen)==1
        out['representative_assumptions']={k:v for k,v in chosen[0].items() if k not in ['low_pH','high_pH','conditional_log10_affinity_ratio']}
        out['representative_human_acid_proxy_kcal']=chosen[0]['low_pH']['conditional_binding_energy_kcal']
        out['representative_human_neutral_proxy_kcal']=chosen[0]['high_pH']['conditional_binding_energy_kcal']
    return out

def peptide_audit(path,expected):
    chain=PDBParser(QUIET=True).get_structure('x',path)[0]['B']; rs=list(chain)
    observed=''.join(seq1(r.resname,custom_map={'HIE':'H','HID':'H','HIP':'H','CYX':'C','ASH':'D','GLH':'E'}) for r in rs)
    twists=[];omega=[]
    for a,b in zip(rs,rs[1:]):
        w=float(np.rad2deg(calc_dihedral(a['CA'].get_vector(),a['C'].get_vector(),b['N'].get_vector(),b['CA'].get_vector())))
        omega.append(abs(abs(w)-180))
        if abs(abs(w)-180)>30:twists.append({'before':a.id[1],'after':b.id[1],'after_resname':b.resname,'omega_degrees':w})
    return {'sequence_match':observed==expected,'peptide_bonds':len(omega),'max_trans_omega_deviation_deg':max(omega),'bonds_over_30_degree_deviation':twists,'qualification':'Omega check only, not a complete Ramachandran or folding assessment.'}

def main():
    files=sorted(p for p in (S/'intermediate/designs').glob('*.json') if not p.stem.endswith('_summary'))
    designs={p.stem:load(p) for p in files}; groups=defaultdict(list); rows=[];errors=[]
    for cid,d in designs.items():
        seq=d['sequence'];groups[seq].append(cid)
        checks={'standard_amino_acids':bool(re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',seq)),'length_allowed':10<=len(seq)<=250,'structure_sequence_match':seq==''.join(r['aa'] for r in d['structure']),'cdr_slices_match':[seq[a:b] for a,b in d['cdr_intervals_zero_based']]==d['cdr_sequences'],'all_new_loop_mask_consistent':all(all(d['structure'][i]['loop'] for i in range(a,b)) for a,b in d['cdr_intervals_zero_based'])}
        if not all(checks.values()):errors.append({'candidate_id':cid,'checks':checks})
        rows.append({'candidate_id':cid,'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'length':len(seq),'branch':d['branch'],'pose_id':d['pose_id'],'checks':checks})
    assert not errors,errors
    duplicate_groups=[{'sequence_sha256':hashlib.sha256(s.encode()).hexdigest(),'candidate_ids':v} for s,v in groups.items() if len(v)>1]
    evs=[load(p) for p in (S/'intermediate/evaluation').glob('*.json')]
    mm={};torsions={};audits={}
    for p in sorted((S/'intermediate/refined').glob('*_human6ARU_proton.json')):
        cid=p.name.removesuffix('_human6ARU_proton.json');mm[cid]=ensemble(load(p));ap=p.with_name(cid+'_human6ARU_audit.json');audits[cid]=load(ap)
        torsions[cid]=peptide_audit(p.with_name(cid+'_human6ARU_relaxed.pdb'),designs[cid]['sequence'])
    assert all(a['sequence_match'] for a in torsions.values())
    mouse={}
    for p in sorted((S/'intermediate/refined').glob('*_mouseAF_proton.json')):
        cid=p.name.removesuffix('_mouseAF_proton.json');v=ensemble(load(p))
        for k in ['representative_human_acid_proxy_kcal','representative_human_neutral_proxy_kcal']:
            if k in v:v[k.replace('human','mouse')]=v.pop(k)
        mouse[cid]=v
    variants={}
    for folder,pattern in [('capped','*_proton.json'),('acid_ensemble','*_acid_histidine.json'),('common_sites','*_proton.json'),('refined','*_acid_proton.json')]:
        for p in sorted((S/'intermediate'/folder).glob(pattern)):
            d=load(p)
            if 'proton_model' in d:variants[str(p.relative_to(S))]=ensemble(d)
    unbound={}
    for p in sorted((S/'intermediate/unbound').glob('*_free_audit.json')):
        d=load(p);x={'geometry_pass_after_minimization':d['geometry_audit']['pass'],'minimization_displacements':d['displacements'],'physical_potential_relaxation_kcal':d['potential_relaxation_kcal'],'potential_relaxation_is_not_binding_free_energy':True}
        if 'NVE' in d:
            n=d['NVE'];f=n['frames'];x['short_NVE']={'duration_ps':n['steps']*n['timestep_ps'],'frames':len(f),'all_geometry_checks_pass':all(z['geometry_pass'] for z in f),'last_loop_RMSD_A':f[-1]['loop_CA_RMSD_A'],'last_CDR_RMSD_A':f[-1]['CDR_CA_RMSD_A'],'last_framework_RMSD_A':f[-1]['framework_CA_RMSD_A'],'max_CA_displacement_A':max(z['maximum_CA_displacement_A'] for z in f),'initial_to_1ps_total_energy_change_kJ':f[1]['total_kj']-f[0]['total_kj'],'post_1ps_energy_range_kJ':max(z['total_kj'] for z in f[1:])-min(z['total_kj'] for z in f[1:]),'initial_to_final_kinetic_energy_ratio':f[-1]['kinetic_kj']/f[0]['kinetic_kj'],'interpretation':'Substantial loop motion is a warning, not stability validation. No thermostat/equilibration; not an equilibrated 300 K trajectory. Initial energy transient and short duration limit interpretation.'}
        unbound[d['candidate_id']]=x
    surface={p.stem:load(p) for p in sorted((S/'intermediate/surface').glob('*.json'))}
    for d in surface.values():d.pop('per_binder_residue',None)
    blast={};sel=load(S/'reference/blast_selected_status.json')
    for db in sel['databases']:
        p=S/'intermediate/novelty'/f"{db['database']}_selected_hits.tsv";assert sha(p)==db['result_sha256'];hits=defaultdict(list);queries=[];processed=None
        for line in p.read_text().splitlines():
            if line.startswith('# Query: '):queries.append(line[9:].split()[0])
            elif line.startswith('# BLAST processed '):processed=int(line.split()[3])
            elif line and not line.startswith('#'):
                a=line.split('\t');hits[a[0]].append({'subject_id':a[1],'percent_identity':float(a[2]),'alignment_length':int(a[3]),'query_length':int(a[4]),'query_coverage_percent':float(a[6]),'evalue':float(a[7]),'bit_score':float(a[8])})
        assert set(queries)==set(sel['ids']);assert processed==len(sel['ids']);assert db['exit_code']==0
        blast[db['database']]={'query_records':len(queries),'processed_records':processed,'queries_with_hits':len(hits),'max_reported_identity_percent':max(x['percent_identity'] for v in hits.values() for x in v),'exact_full_query_matches_in_reported_hits':sum(any(x['percent_identity']==100 and x['query_coverage_percent']==100 for x in v) for v in hits.values()),'top_scored_hits':{cid:max(v,key=lambda x:x['bit_score']) for cid,v in hits.items()},'qualification':'BLAST top-10 output, not an exhaustive all-hit maximum. Framework similarity is expected; antibody-specific CDR novelty rules apply.'}
    novelty=load(S/'reference/novelty_screen_augmented_stats.json')
    summary={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':3,'not_final_submission':True,'phone_a_friend_used':2,'total_sequence_structure_records':len(designs),'unique_sequences':len(groups),'duplicate_records_beyond_unique_sequences':len(designs)-len(groups),'duplicate_sequence_groups':duplicate_groups,'unique_CDR3':len({d['cdr_sequences'][2] for d in designs.values()}),'length_range':[min(len(s) for s in groups),max(len(s) for s in groups)],'sequence_structure_checks_pass':True,'format_checks_per_record':5,'initial_geometry_pass':sum(e['structure_audit']['geometry_prefilter_pass'] for e in evs),'no_glycosylation_sequon':sum(not e['sequence_audit']['liabilities']['glycosylation_sequons'] for e in evs),'sequon_flagged_ids':[e['candidate_id'] for e in evs if e['sequence_audit']['liabilities']['glycosylation_sequons']],'native_human_detailed_evaluations':len(mm),'native_human_acid_on_all_54':sum(v['all_scenarios_acid_on'] for v in mm.values()),'native_refinement_geometry_pass':sum(v['pass'] for v in audits.values()),'native_refinement_full_receptor_glycan_pass':sum(v.get('postrefinement_exclusion_pass',False) for v in audits.values()),'native_human_backbones_with_omega_flag':sum(bool(v['bonds_over_30_degree_deviation']) for v in torsions.values()),'mouse_detailed_evaluations':len(mouse),'local_cdr3_screen':novelty,'completed_whole_sequence_blast_records':sel['records'],'completed_whole_sequence_blast_unique_sequences':sel['unique_sequences'],'no_experimental_binding_or_folding_validation':True}
    for name,data in [('design_integrity_audit.json',rows),('detailed_human_comparison.json',mm),('detailed_mouse_comparison.json',mouse),('geometry_variant_comparison.json',variants),('peptide_omega_audit.json',torsions),('unbound_audit_summary.json',unbound),('surface_audit_summary.json',surface),('blast_audit_summary.json',blast),('stage3_summary.json',summary)]:
        (OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False))
    unique=OUT/'stage3_unique_sequences_NOT_FOR_SUBMISSION.fasta'
    with unique.open('w') as f:
        for s,cids in sorted(groups.items(),key=lambda x:x[1][0]):f.write('>'+cids[0]+' NOT_FOR_SUBMISSION aliases='+','.join(cids)+'\n'+s+'\n')
    (OUT/'README.txt').write_text('INTERIM SCIENTIFIC CHECKPOINT, NOT A SUBMISSION PACKAGE.\nThe FASTA is unranked, includes weak/rejected designs, and is not the requested final 100.\nDo not upload it to the competition. Duplicate sequence hypotheses are collapsed.\nAll energy values are conditional model diagnostics, not measured KD or binding free energies.\nSee stage3/REPORT.md and stage3/STATE.md for limitations and next work.\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['duplicate_sequence_groups','local_cdr3_screen']},indent=2));print('wrote',OUT)
if __name__=='__main__':main()
