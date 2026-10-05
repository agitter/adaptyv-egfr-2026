"""Independent sequence, stereochemistry, target/glycan and proton-model gates."""
from pathlib import Path
import sys,json,time,re
R=Path('/mnt/data/egfr_campaign');S3=R/'stage3'
sys.path[:0]=[str(R/'stage2/code'),str(R/'code')]
from evaluate_pool import sequence_audit,structural_audit,contact_audit,TARGETS
from proton_polynomial import evaluate_design

def coarse_rank_eligible(d):
    """A transformed candidate cannot inherit an unrecomputed parent score."""
    return (not d.get('parent_metrics_not_recomputed')
            and d.get('metrics', {}).get('scores_valid', True) is not False)

def main():
    dest=S3/'intermediate/evaluation';dest.mkdir(parents=True,exist_ok=True);start=time.time();count=0
    for p in sorted((S3/'intermediate/designs').glob('*.json')):
        if p.stem.endswith('_summary'):continue
        out=dest/p.name
        if out.exists():continue
        d=json.loads(p.read_text());a=sequence_audit(d);g=structural_audit(d);c=contact_audit(d)
        e={'candidate_id':d['candidate_id'],'branch':d['branch'],'sequence_audit':a,'structure_audit':g,'candidate_ion_pairs':c}
        if a['format_pass'] and g['geometry_prefilter_pass']:
            try:
                e['human_proton_model']=evaluate_design(d['structure'],TARGETS['human6ARU'],True)
                e['mouse_proton_model']=evaluate_design(d['structure'],TARGETS['mouseAF'],False)
            except Exception as exc:e['proton_model_error']=repr(exc)
        out.write_text(json.dumps(e));count+=1
        if count%50==0:print('evaluated',count,'seconds',round(time.time()-start,1),flush=True)
    files=list(dest.glob('*.json'));es=[json.loads(p.read_text()) for p in files];surv=[e for e in es if e['structure_audit']['geometry_prefilter_pass']]
    clean=[e for e in surv if not e['sequence_audit']['liabilities']['glycosylation_sequons']]
    def rank(e):
        d=json.loads((S3/'intermediate/designs'/(e['candidate_id']+'.json')).read_text())
        p=e.get('human_proton_model',{}).get('worst_contrast_kcal_surrogate',-9)
        # Conservative opportunity score, not a calibrated probability.
        return 8*p-d['metrics']['human_acid_surrogate']*.15-d['metrics']['packing_regularizer']*.15
    rankable=[e for e in clean if coarse_rank_eligible(json.loads((S3/'intermediate/designs'/(e['candidate_id']+'.json')).read_text()))]
    order=sorted(rankable,key=rank,reverse=True);chosen=[];perpose={}
    for e in order:
        pose=e['candidate_id'].split('_')[0]
        if perpose.get(pose,0)>=1:continue
        chosen.append(e['candidate_id']);perpose[pose]=1
        if len(chosen)>=36:break
    (S3/'reference/mm_selection.json').write_text(json.dumps(chosen,indent=2))
    summary={'evaluated':len(es),'geometry_pass':len(surv),'geometry_and_no_sequon':len(clean),'coarse_all_12_acid_on':sum(e.get('human_proton_model',{}).get('all_scenarios_acid_on',False) for e in surv),'MM_selection':chosen,'excluded_from_coarse_ranking_invalid_parent_scores':len(clean)-len(rankable),'ranking_qualification':'Exploratory selection only; not final ranking or binding validation.','seconds_this_run':time.time()-start}
    (S3/'reference/evaluation_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
if __name__=='__main__':main()
