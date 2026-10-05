"""Independent rule/coordinate/ancestry gates and predeclared selection.
Neither a geometry pass nor a novel CDR is an experimental success claim.
"""
from pathlib import Path
import sys,json,re,hashlib,copy,math
import numpy as np
R=Path('/mnt/data/egfr_campaign');S5=R/'stage5'
sys.path[:0]=[str(S5/'code'),str(R/'stage4/code'),str(R/'stage3/code'),str(R/'stage2/code'),str(R/'code'),str(R/'stage2/runtime/python')]
from build_return_loops import load_parent,PARENT,SOURCE,sha
from evaluate_pool import structural_audit,contact_audit,sequence_audit
from legacy_geometry import dihedral

def main():
 dest=S5/'intermediate/evaluation';dest.mkdir(parents=True,exist_ok=True);summary=[];source_design,original,tx=load_parent();seen=set()
 for p in sorted((S5/'intermediate/designs').glob('[RN]*.json')):
  d=json.loads(p.read_text());s=d['sequence'];a=sequence_audit(d);g=structural_audit(d);identity=s==''.join(r['aa'] for r in d['structure']) and d['cdr_sequences']==[s[x:y] for x,y in d['cdr_intervals_zero_based']]
  hashes=sha(d['ancestry']['parent_file'])==d['ancestry']['parent_sha256'] and sha(d['ancestry']['relaxed_backbone_source'])==d['ancestry']['relaxed_backbone_sha256'];framework=all(r['aa']==original[r['parent_position']-1]['aa'] and all(np.max(np.abs(np.array(r['atoms'][n])-original[r['parent_position']-1]['atoms'][n]))<1e-10 for n in ['N','CA','C']) for r in d['structure'] if not r['loop']);protected={i:next(r for r in d['structure'] if r.get('parent_position')==i) for i in [54,100,103]};protection=all(protected[i]['aa']==original[i-1]['aa'] and protected[i]['atoms']==original[i-1]['atoms'] for i in protected)
  motif=re.findall(r'C([ACDEFGHIKLMNPQRSTVWY]{5,40}?)[WF][GAS]QG',s);cdr_match=motif==[d['cdr_sequences'][2]]
  chain=[];omega=[]
  for x,y in zip(d['structure'],d['structure'][1:]):
   ax={k:np.array(v) for k,v in x['atoms'].items()};ay={k:np.array(v) for k,v in y['atoms'].items()};theta=math.degrees(dihedral(ax['CA'],ax['C'],ay['N'],ay['CA']))
   if abs(180-abs(theta))>25:omega.append({'previous_position':x['position'],'omega_deg':theta})
  duplicate=s in seen;seen.add(s);passed=bool(a['format_pass'] and identity and hashes and framework and protection and cdr_match and not duplicate and g['geometry_prefilter_pass'] and not a['liabilities']['glycosylation_sequons'] and not a['liabilities']['hydrophobic_runs'] and not omega and d['loop_rebuild']['filter_status']=='candidate')
  result={'candidate_id':d['candidate_id'],'sequence_audit':a,'structure_audit':g,'sequence_coordinates_match':identity,'parent_hashes_match':hashes,'framework_unchanged':framework,'protected_contact_atoms_unchanged':protection,'motif_CDR3_matches_design':cdr_match,'duplicate_sequence':duplicate,'peptide_omega_warnings':omega,'pass':passed,'contact_audit':contact_audit(d)};(dest/p.name).write_text(json.dumps(result,indent=2));summary.append({'candidate_id':d['candidate_id'],'branch':d['branch'],'region':d['loop_rebuild']['region'],'new_segment_length':d['loop_rebuild']['new_length'],'length':len(s),'pass':passed,'packing_objective':d['loop_rebuild']['packing_objective'],'backbone_index':d['ancestry']['backbone_index'],'sequence':s,'glycosylation_sequons':a['liabilities']['glycosylation_sequons'],'net_formal_charge':a['liabilities']['net_formal_charge_without_histidine'],'glycan_min_A':g['resolved_glycans']['min_heavy_distance_A'],'1IVO_clashes':g['1IVO']['heavy_atoms_below_2A'],'1NQL_clashes':g['1NQL']['heavy_atoms_below_2A'],'omega_warnings':len(omega)})
 (S5/'reference/design_audit_summary.json').write_text(json.dumps(summary,indent=2));print('TOTAL',len(summary),'PASS',sum(r['pass'] for r in summary))
 for prefix in ['N','R']:
  for region,n in [('H2',7),('H2',9),('H3',6),('H3',8)]:
   pool=[r for r in summary if r['candidate_id'].startswith(prefix) and r['region']==region and r['new_segment_length']==n and r['pass']];print(prefix,region,n,'n',len(pool),'best',[(r['candidate_id'],round(r['packing_objective'],2)) for r in sorted(pool,key=lambda x:x['packing_objective'])[:3]])
 # Distinct geometries in each branch; predeclare selections before all-atom outcomes.
 selected=[]
 for prefix,quota in [('N',2),('R',1)]:
  for region,n in [('H2',7),('H2',9),('H3',6),('H3',8)]:
   pool=sorted([r for r in summary if r['candidate_id'].startswith(prefix) and r['region']==region and r['new_segment_length']==n and r['pass']],key=lambda r:(r['1IVO_clashes']+r['1NQL_clashes'],r['packing_objective']));geoms=set();count=0
   for r in pool:
    if r['backbone_index'] in geoms:continue
    geoms.add(r['backbone_index']);selected.append(r);count+=1
    if count>=quota:break
 # Unchanged input gets exactly the same next joint minimization as mutants.
 parent=copy.deepcopy(source_design);parent.update(candidate_id='B00000',structure=original,seed=None,branch='unchanged_equal_refinement_control',campaign_stage=5,ancestry={'parent_candidate':'C00003','parent_file':str(PARENT),'parent_sha256':sha(PARENT),'relaxed_backbone_source':str(SOURCE),'relaxed_backbone_sha256':sha(SOURCE)},metrics={'scores_valid':False,'note':'Control only; not a new sequence.'});(S5/'intermediate/designs/B00000.json').write_text(json.dumps(parent))
 (S5/'reference/human_selection.json').write_text(json.dumps({'rule':'Two neutral and one charged per loop-length group, distinct geometries, initial geometry gates, alternative-state overlap first then coarse packing. Unchanged control same refinement. Selected before all-atom evaluation.','selected':selected,'control':'B00000'},indent=2))
 print('SELECTED',[r['candidate_id'] for r in selected],'+ B00000 control')
if __name__=='__main__':main()
