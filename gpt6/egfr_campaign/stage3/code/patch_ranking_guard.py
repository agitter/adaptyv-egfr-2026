"""Prevent unrecomputed parent scores from entering exploratory selection.
This changes no saved candidate, geometry, or energy result.
"""
from pathlib import Path
import hashlib,json,datetime
S=Path('/mnt/data/egfr_campaign/stage3')
p=S/'code/evaluate_new.py'
old=p.read_text()
backup=S/'reference/evaluate_new_before_ranking_guard.py.txt'
if not backup.exists(): backup.write_text(old)
needle='def main():\n'
helper='''def coarse_rank_eligible(d):
    """A transformed candidate cannot inherit an unrecomputed parent score."""
    return (not d.get('parent_metrics_not_recomputed')
            and d.get('metrics', {}).get('scores_valid', True) is not False)

'''
if 'def coarse_rank_eligible' not in old:
    new=old.replace(needle,helper+needle,1)
    new=new.replace('order=sorted(clean,key=rank,reverse=True);chosen=[];perpose={}',
        "rankable=[e for e in clean if coarse_rank_eligible(json.loads((S3/'intermediate/designs'/(e['candidate_id']+'.json')).read_text()))]\n    order=sorted(rankable,key=rank,reverse=True);chosen=[];perpose={}")
    new=new.replace("'MM_selection':chosen,'seconds_this_run'", "'MM_selection':chosen,'excluded_from_coarse_ranking_invalid_parent_scores':len(clean)-len(rankable),'ranking_qualification':'Exploratory selection only; not final ranking or binding validation.','seconds_this_run'")
    p.write_text(new)
for name in ['mm_selection.json','evaluation_summary.json']:
    src=S/'reference'/name; dst=S/'reference'/('before_ranking_guard_'+name)
    if not dst.exists():dst.write_bytes(src.read_bytes())
import sys
sys.path.insert(0,str(S/'code'))
from evaluate_new import coarse_rank_eligible
cases=[({},True),({'parent_metrics_not_recomputed':{'value':3}},False),({'metrics':{'scores_valid':False}},False)]
for cid,expected in [('G00117_C04',True),('A00001',False),('B00013',False),('C00003',False)]:
    d=json.loads((S/'intermediate/designs'/f'{cid}.json').read_text());cases.append((d,expected))
for d,expected in cases:assert coarse_rank_eligible(d)==expected
record={'reason':'Exclude inherited or placeholder score terms from exploratory automatic MM selection. Actual manually selected refinement batches preserved unchanged.',
'original_sha256':hashlib.sha256(backup.read_bytes()).hexdigest(),'updated_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'tests_passed':len(cases),'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(S/'reference/ranking_guard_patch.json').write_text(json.dumps(record,indent=2));print(record)
