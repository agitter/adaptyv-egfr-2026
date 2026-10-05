"""Common low-cost triage for geometry survivors. Not a binding validation.
Reuses the campaign's historical charge/occlusion binding polynomial. This
coarse model previously overpredicted pH selectivity; detailed results override it.
"""
from pathlib import Path
import json,sys,time
R=Path('/mnt/data/egfr_campaign');S=R/'stage6'
sys.path[:0]=[str(R/'stage2/code'),str(R/'code')]
from evaluate_pool import TARGETS
from proton_polynomial import evaluate_design
start=time.time();rows=json.loads((S/'reference/campaign_audit.json').read_text());out=S/'intermediate/coarse_common';out.mkdir(exist_ok=True);records=[];seen=set()
for row in rows:
 if not (row['initial_context_pass'] and row['local_novelty_pass']) or row['sequence'] in seen:continue
 seen.add(row['sequence']);key=row['key'];path=out/(key.replace(':','_')+'.json')
 if path.exists():d=json.loads(path.read_text())
 else:
  design=json.loads(Path(row['design_source']).read_text());d={'key':key,'sequence':row['sequence'],'source':row['design_source'],'human':evaluate_design(design['structure'],TARGETS['human6ARU'],True),'mouse':evaluate_design(design['structure'],TARGETS['mouseAF'],False),'qualification':'Coarse triage only; known false-positive tendency. Cannot override detailed failure.'};path.write_text(json.dumps(d,indent=2))
 records.append({'key':key,'sequence':row['sequence'],'source':str(path),'minimum_contrast':d['human']['worst_contrast_kcal_surrogate'],'all_acid_on':d['human']['all_scenarios_acid_on']})
 if len(records)%25==0:print('Coarse common',len(records),round(time.time()-start,1),flush=True)
(S/'reference/common_coarse_summary.json').write_text(json.dumps({'records':records,'total':len(records),'acid_on':sum(r['all_acid_on'] for r in records),'seconds':time.time()-start},indent=2));print('FINISHED',len(records),round(time.time()-start,1),flush=True)
