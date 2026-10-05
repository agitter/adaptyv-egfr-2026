"""Matched physical follow-up of one existing sequence per alternative docking pose.
Selection is made by the already fixed coarse ranking within each alternative
pose before new detailed results are available. These are not new sequences.
"""
from pathlib import Path
import json,collections,copy,hashlib,subprocess,sys,os
S=Path('/mnt/data/egfr_campaign/stage6');a=json.loads((S/'reference/campaign_audit.json').read_text());co={r['key']:r for r in json.loads((S/'reference/common_coarse_summary.json').read_text())['records']};byp=collections.defaultdict(list)
for r in a:
 c=co.get(r['key'])
 if c and c['all_acid_on'] and r['pose_id']!='G00117':byp[r['pose_id']].append((c['minimum_contrast'],r))
selected=[];jobs=[]
for i,(pose,vals) in enumerate(sorted(byp.items()),1):
 score,r=max(vals,key=lambda x:(x[0],x[1]['key']));src=Path(r['design_source']);d=json.loads(src.read_text());cid=f'K{i:05d}';d.update(candidate_id=cid,campaign_stage=6,reassessment_source={'candidate':r['key'],'file':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'sequence_changed':False,'purpose':'Diverse pose physical reassessment; existing sequence, not a newly generated design.'});p=S/'intermediate/designs'/f'{cid}.json';p.write_text(json.dumps(d));selected.append({'candidate_id':cid,'parent_key':r['key'],'pose':pose,'sequence':d['sequence'],'coarse_selection_contrast':score})
 for species in ['human6ARU','mouseAF']:jobs.append({'id':cid+'_'+species,'command':[sys.executable,str(S/'code/refine_one.py'),cid,species],'timeout':600})
plan={'name':'diversity_followup','concurrency':2,'selection':selected,'jobs':jobs};p=S/'reference/diversity_followup_plan.json';p.write_text(json.dumps(plan,indent=2));env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1');f=open(S/'logs/diversity_followup_manager.log','w');proc=subprocess.Popen([sys.executable,str(S/'code/run_batch.py'),str(p)],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,env=env);(S/'reference/diversity_followup_pid.json').write_text(json.dumps({'pid':proc.pid}));print('Started',len(selected),'existing sequences',len(jobs),'models',proc.pid)
for x in selected:print(x['candidate_id'],x['parent_key'],x['pose'])
